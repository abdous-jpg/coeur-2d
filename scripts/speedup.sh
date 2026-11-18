#!/bin/bash
# Mesure le temps de calcul pour 1, 2 et 4 threads (sans ecriture disque).
# Usage, depuis run/ : ../scripts/speedup.sh 512 500
N=${1:-512}
T=${2:-500}
export OMP_WAIT_POLICY=passive

echo "threads,temps_s" > speedup.csv
for t in 1 2 4; do
  temps=$(OMP_NUM_THREADS=$t ../build/coeur2d "$N" "$T" 0 | grep Temps | awk '{print $3}')
  echo "$t,$temps" >> speedup.csv
  echo "$t thread(s) : $temps s"
done
