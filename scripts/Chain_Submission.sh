#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --mem=15G
#SBATCH --account=abrunet1

echo hello world

cd /home/gareeves/WGA/Scripts/

/home/gareeves/WGA/packages/axtChain $THREE $ONE $TWO $FOUR -minScore=2000 -linearGap=medium -psl

