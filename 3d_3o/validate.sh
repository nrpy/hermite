#!/bin/bash

python3 hermite.py
cd AEILocalInterp_output
./convert_to_py.sh .
cd ..
python3 validator.py
