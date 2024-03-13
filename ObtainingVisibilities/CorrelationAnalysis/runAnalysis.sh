#!/bin/bash

#SBATCH --account=PAS1977
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rose.1655@osu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6000
#SBATCH --job-name=analysis
#SBATCH --time=00:15:00

### josie 12 mar 2024 - bash script to run analysis over many zipped frames files at once
### note: you must run this ONCE for EACH observation period directory (ex: once for all files in 2024FebDataVersii, again in 2023DecDataVersii, etc)
### list all the files: ls ../2024FebDataVersii | grep gam | grep cas | grep ZippedFrames.root

source /users/PAS1977/jrose/anaconda3/bin/activate
conda activate RootEnv

declare -a zips=(
    "\"gam cas_y2024m02d18h19m31s45_ZippedFrames.root\""
    "\"gam cas_y2024m02d18h20m38s27_ZippedFrames.root\""
    "\"gam cas_y2024m02d18h21m44s15_ZippedFrames.root\""
    "\"gam cas_y2024m02d19h19m07s37_ZippedFrames.root\""
    "\"gam cas_y2024m02d19h20m14s25_ZippedFrames.root\""
    "\"gam cas_y2024m02d20h19m22s00_ZippedFrames.root\""
    "\"gam cas_y2024m02d21h20m26s13_ZippedFrames.root\""
    "\"gam cas_y2024m02d21h20m59s31_ZippedFrames.root\""
    "\"gam cas_y2024m02d21h21m32s24_ZippedFrames.root\""
    )
    
for i in "${zips[@]}"
do
    echo "$i"
    srun root /users/PAS1977/jrose/macros/versiiAnalysisBothFixNorm.C\("$i",1\)
done
