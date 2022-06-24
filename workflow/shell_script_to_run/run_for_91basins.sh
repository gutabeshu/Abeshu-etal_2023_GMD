#!/bin/bash

# get a listing of the directories to use
# for runoff
cd /workflow/runoff-watch-file/
# for streamflow
cd /workflow/runoff-flow-file/

for dir in `ls -d xanthos* ` ; do
	dirnum=`echo $dir | cut -c 8-`
	DirNum=$dirnum
	
	sbatch --job-name=Basin${DirNum} --output=Basin${DirNum}.out --export=DirNum=$dirnum /workflow/run_for_91basins.sbatch

done

