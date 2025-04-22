#!/bin/bash

cd 1/
./runpy.sh 1
python3 generategandq.py

cd ../2/
./runpy.sh 1
python3 generategandq.py

cd ../3/
./runpy.sh 1
python3 generategandq.py

cd ../
