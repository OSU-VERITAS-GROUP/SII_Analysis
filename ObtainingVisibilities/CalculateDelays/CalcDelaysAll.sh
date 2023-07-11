#!/bin/bash

source ~/anaconda3/bin/activate
conda activate my_root_env

index_array=(
UTC20211217/m12d16h23/delays/
UTC20211217/m12d17h01/delays/
UTC20211217/m12d17h03/delays/
UTC20211217/m12d17h05/delays/
UTC20220212/m02d11h19/delays/
UTC20220212/m02d11h20/delays/
UTC20220212/m02d11h23/delays/
UTC20220212/m02d12h01/delays/
UTC20220212/m02d12h03/delays/
UTC20220214/m02d13h19/delays/
UTC20220214/m02d13h21/delays/
UTC20220214/m02d13h23/delays/
UTC20220214/m02d14h01/delays/
UTC20220214/m02d14h03/delays/
UTC20220313/m03d12h22/delays/
UTC20220313/m03d13h01/delays/
UTC20220511/m05d10h21/delays/
)
 
for i in ${index_array[@]}
do
    cd $i 
    echo $i
    #cp /users/PAS1977/jrose/macros/CalcDelaysTelCoords.py .
    python CalcDelaysEtc.py
    
    cd ~/betUMaReproduce/
done