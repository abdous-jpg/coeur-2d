"""Lit speedup.csv et trace la courbe de speedup (speedup.png)."""
import csv

import matplotlib.pyplot as plt

with open("speedup.csv") as f:
    lignes = list(csv.DictReader(f))

threads = [int(l["threads"]) for l in lignes]
temps = [float(l["temps_s"]) for l in lignes]
speedup = [temps[0] / t for t in temps]

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(threads, threads, "--", color="gray", label="Idéal")
ax.plot(threads, speedup, "o-", label="Mesuré")
for t, s in zip(threads, speedup):
    ax.annotate(f"{s:.2f}", (t, s), textcoords="offset points", xytext=(6, -12))
ax.set_xlabel("Nombre de threads OpenMP")
ax.set_ylabel("Speedup")
ax.set_xticks(threads)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("speedup.png", dpi=120)
print("Courbe enregistree : speedup.png")
