#!/bin/bash
# this scripts setups the virtual environment for python, installs all the libraries
#$ -P fraser.prjc
#$ -q short.qc
#$ -N test_DL_data
#$ -e log
#$ -o log
#$ -cwd
module load Python/3.7.2-GCCcore-8.2.0
python3 -m venv env #local python3 install
source env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
#deactivate
