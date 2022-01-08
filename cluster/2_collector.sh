#!/bin/bash
outputdir=Outputs/synthetic_Gamma
#count how many files in the output dir
rm -f $outputdir/*.txt
#sim=$(ls $outputdir|grep li_lr*.tmp|wc -l)
sim=1000
echo $sim
for i in $(seq 1 $sim) ; do echo "$(tail -n 1 $outputdir/li_$i.tmp)">>$outputdir/li.txt; done
for i in $(seq 1 $sim) ; do echo "$(tail -n 1 $outputdir/lr_$i.tmp)">>$outputdir/lr.txt; done
for i in $(seq 1 $sim) ; do echo "$(tail -n 1 $outputdir/li_lr_$i.tmp)">>$outputdir/li_lr.txt; done
rm -f outputdir/*.tmp
