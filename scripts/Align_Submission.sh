#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --mem=15G
#SBATCH --account=abrunet1

cd /home/gareeves/WGA/Scripts/

/home/gareeves/WGA/packages/lastz_D $ONE $TWO --gap=400,30 --hspthresh=3000 --inner=0 --ydrop=9400 --gappedthresh=3000 --masking=50 --scores=/home/gareeves/WGA/Scripts/HoxD55.q > $THREE

