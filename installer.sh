#!/bin/bash
git clone https://github.com/AlexandrovLab/SigProfilerMatrixGenerator.git
mkdir input output 
cd SigProfilerMatrixGenerator 
python3 install_catalogue_ref.py
