#!/bin/bash
#SBATCH -J emu10_baoan_all
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=30
#SBATCH --time=50-00:00:00
#SBATCH --mem=204800
#SBATCH -p bigmem
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
emu -m -p GLs.merged.chr1-22.redup.0.05 -e 10 -t 30 -o GLs.merged.chr1-22.redup.0.05.emu10
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
