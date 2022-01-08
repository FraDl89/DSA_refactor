#!/bin/bash
#$ -P fraser.prjc
#$ -q short.qc
#$ -N test_DL_data
#$ -e log
#$ -o log
#$ -t 1-5:1
#$ -pe shmem 3
#$ -cwd
#This script runs the algorithms on the cluster
sim=$SGE_TASK_ID
module load Python/3.6.6-foss-2018b
source env/bin/activate
outputdir=Outputs/synthetic_Gamma
mkdir -p $outputdir
python3 1_Generate_data.py $sim $outputdir
python3 2_likelihood_li.py $sim  > $outputdir/li_$sim.tmp
python3 3_likelihood_lr.py $sim > $outputdir/lr_$sim.tmp
python3 4_likelihood_li_lr.py $sim > $outputdir/li_lr_$sim.tmp
deactivate
rm -f $outputdir/Gamma_distr_seed$sim.npy
