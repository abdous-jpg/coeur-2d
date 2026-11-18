"""Lit sortie.h5 et fabrique spirale.gif (aperçu rapide sans ParaView)."""
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

with h5py.File("sortie.h5", "r") as f:
    noms = sorted(f["u"].keys())
    images = [f["u"][nom][...] for nom in noms]

fig, ax = plt.subplots(figsize=(4, 4))
im = ax.imshow(images[0], cmap="inferno", vmin=0, vmax=1, origin="lower")
ax.axis("off")


def maj(k):
    im.set_data(images[k])
    return [im]


anim = FuncAnimation(fig, maj, frames=len(images), blit=True)
anim.save("spirale.gif", writer=PillowWriter(fps=10))
print(f"{len(images)} images -> spirale.gif")
