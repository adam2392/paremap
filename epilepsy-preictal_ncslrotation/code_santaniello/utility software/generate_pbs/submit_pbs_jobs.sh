#!/bin/bash
#-------------------------------------------------------------------------------
# bash code for submitting multiple pbs jobs in a local folder
#-------------------------------------------------------------------------------

# get the list of pbs jobs to be submitted. The jobs are in the local folder
ls *.pbs >list_pbs.txt

# read the list and submit each job...
while read LINE
do
  qsub $LINE
  echo -e "job $LINE has been submitted\n"
done <list_pbs.txt
#-------------------------------------------------------------------------------

