#!/bin/bash

args=("$@")

if [ "${args}" == "1" ];then

    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    
else

    export MKL_NUM_THREADS=12
    export NUMEXPR_NUM_THREADS=12
    export OMP_NUM_THREADS=12
    
fi

python3 dd.py 
