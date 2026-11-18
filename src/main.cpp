// Propagation electrique dans le tissu cardiaque 2D
// Modele d'Aliev-Panfilov, differences finies, Euler explicite.
// Version 1 : sequentielle (reference pour OpenMP puis MPI).

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <hdf5.h>
#include <string>
#include <vector>

// Parametres du modele (valeurs classiques de la litterature)
constexpr double k    = 8.0;
constexpr double a    = 0.15;
constexpr double eps0 = 0.002;
constexpr double mu1  = 0.2;
constexpr double mu2  = 0.3;
constexpr double D    = 1.0;   // coefficient de diffusion

struct Grille {
    int n;                     // points interieurs par direction
    std::vector<double> data;  // (n+2) x (n+2), avec cellules fantomes
    explicit Grille(int n_) : n(n_), data((n_ + 2) * (n_ + 2), 0.0) {}
    double& operator()(int i, int j) { return data[i * (n + 2) + j]; }
    double operator()(int i, int j) const { return data[i * (n + 2) + j]; }
};

// Bords isoles (Neumann, flux nul) : on recopie la valeur voisine
void conditions_bord(Grille& u) {
    const int n = u.n;
    for (int i = 1; i <= n; ++i) {
        u(i, 0) = u(i, 1);
        u(i, n + 1) = u(i, n);
    }
    for (int j = 0; j <= n + 1; ++j) {
        u(0, j) = u(1, j);
        u(n + 1, j) = u(n, j);
    }
}

// Condition initiale S1-S2 : un front plan + une zone refractaire
// La rencontre des deux casse le front et cree une spirale.
void initialiser(Grille& u, Grille& v) {
    const int n = u.n;
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j) {
            u(i, j) = (j <= n / 20) ? 1.0 : 0.0;   // bande excitee a gauche
            v(i, j) = (i > n / 2 && j <= n / 2) ? 1.0 : 0.0; // zone refractaire
        }
}

// Sortie HDF5 : un fichier, un dataset par instantane (/u/0000, /u/0001...)
// + un fichier XDMF qui decrit la grille pour ParaView.
class SortieHDF5 {
public:
    SortieHDF5(const std::string& base, int n, double dx)
        : base_(base), n_(n), dx_(dx), buf_(std::size_t(n) * n) {
        fichier_ = H5Fcreate((base + ".h5").c_str(), H5F_ACC_TRUNC,
                             H5P_DEFAULT, H5P_DEFAULT);
        hid_t g = H5Gcreate2(fichier_, "/u", H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT);
        H5Gclose(g);
    }

    void ecrire(const Grille& u, double t) {
        // On copie l'interieur (sans cellules fantomes) dans un buffer contigu
        for (int i = 1; i <= n_; ++i)
            for (int j = 1; j <= n_; ++j)
                buf_[std::size_t(i - 1) * n_ + (j - 1)] = u(i, j);

        char nom[32];
        std::snprintf(nom, sizeof nom, "/u/%04d", int(temps_.size()));
        hsize_t dims[2] = {hsize_t(n_), hsize_t(n_)};
        hid_t espace = H5Screate_simple(2, dims, nullptr);
        hid_t ds = H5Dcreate2(fichier_, nom, H5T_NATIVE_DOUBLE, espace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buf_.data());
        H5Dclose(ds);
        H5Sclose(espace);
        temps_.push_back(t);
    }

    ~SortieHDF5() {
        H5Fclose(fichier_);
        ecrire_xdmf();
    }

private:
    void ecrire_xdmf() const {
        std::ofstream f(base_ + ".xmf");
        const std::string h5 = base_ + ".h5";
        f << "<?xml version=\"1.0\" ?>\n<Xdmf Version=\"3.0\">\n<Domain>\n"
          << "<Grid Name=\"coeur\" GridType=\"Collection\" "
             "CollectionType=\"Temporal\">\n";
        for (std::size_t s = 0; s < temps_.size(); ++s) {
            char nom[32];
            std::snprintf(nom, sizeof nom, "/u/%04d", int(s));
            f << "  <Grid Name=\"t" << s << "\" GridType=\"Uniform\">\n"
              << "    <Time Value=\"" << temps_[s] << "\"/>\n"
              // Grille 3D d'epaisseur 1 (ordre Z Y X) : ParaView l'affiche dans le plan X-Y
              << "    <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"1 "
              << n_ << " " << n_ << "\"/>\n"
              << "    <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
              << "      <DataItem Dimensions=\"3\" Format=\"XML\">0 0 0</DataItem>\n"
              << "      <DataItem Dimensions=\"3\" Format=\"XML\">1 " << dx_
              << " " << dx_ << "</DataItem>\n"
              << "    </Geometry>\n"
              << "    <Attribute Name=\"u\" AttributeType=\"Scalar\" Center=\"Node\">\n"
              << "      <DataItem Dimensions=\"1 " << n_ << " " << n_
              << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">"
              << h5 << ":" << nom << "</DataItem>\n"
              << "    </Attribute>\n  </Grid>\n";
        }
        f << "</Grid>\n</Domain>\n</Xdmf>\n";
    }

    std::string base_;
    int n_;
    double dx_;
    hid_t fichier_;
    std::vector<double> buf_;
    std::vector<double> temps_;
};

int main(int argc, char** argv) {
    const int    n      = (argc > 1) ? std::atoi(argv[1]) : 256;
    const double T      = (argc > 2) ? std::atof(argv[2]) : 1000.0;
    const int    nsnap  = (argc > 3) ? std::atoi(argv[3]) : 50;
    const double dx     = 0.5;
    const double dt     = 0.04;  // stabilite : dt <= dx^2 / (4 D) = 0.0625
    const long   nsteps = static_cast<long>(T / dt);
    const double alpha  = D * dt / (dx * dx);

    Grille u(n), v(n), u_new(n);
    initialiser(u, v);

    std::printf("Grille %d x %d, %ld pas de temps\n", n, n, nsteps);
    auto t0 = std::chrono::steady_clock::now();
    const long periode = std::max(1L, nsteps / std::max(1, nsnap));
    SortieHDF5 sortie("sortie", n, dx);

    for (long s = 0; s < nsteps; ++s) {
        conditions_bord(u);

        #pragma omp parallel for
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j) {
                const double uc = u(i, j), vc = v(i, j);
                const double lap = u(i - 1, j) + u(i + 1, j) + u(i, j - 1)
                                 + u(i, j + 1) - 4.0 * uc;
                const double eps = eps0 + mu1 * vc / (uc + mu2);
                u_new(i, j) = uc + alpha * lap
                            + dt * (-k * uc * (uc - a) * (uc - 1.0) - uc * vc);
                v(i, j) = vc + dt * eps * (-vc - k * uc * (uc - a - 1.0));
            }
        std::swap(u.data, u_new.data);

        if (nsnap > 0 && s % periode == 0)
            sortie.ecrire(u, s * dt);
    }

    auto t1 = std::chrono::steady_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    double mcups = 1e-6 * n * double(n) * nsteps / sec;
    std::printf("Temps : %.2f s   Debit : %.1f Mcellules/s\n", sec, mcups);
    return 0;
}
