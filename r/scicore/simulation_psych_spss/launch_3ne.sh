#!/bin/bash

#SBATCH --job-name=ah_450                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem=8G
#Total memory reserved: 2GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=out/out.o     #These are the STDOUT and STDERR files
#SBATCH --error=err/err.o

#You selected an array of jobs from 1 to 108 (108 population models)


#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
#################################
ml R/4.0.3-foss-2018b

#export your required environment variables below
#################################################
cp $HOME/simulation_psych_spss/output/model_recovery_results.RDS $TMPDIR/

#add your command lines below
#############################
time Rscript 3_analysis_NE.R path=\"$TMPDIR\"
