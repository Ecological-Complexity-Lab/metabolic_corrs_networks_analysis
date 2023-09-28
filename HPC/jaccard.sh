#! /bin/bash

# these are variables to be used in the job queueing by the HPC:
#$ -q shai.q
#$ -cwd
#$ -N jaccard_rnd_vs_true_pos_and_neg
#$ -l h_vmem=2G

# running the desired R script
/gpfs0/shai/projects/R4/R-4.2.0/bin/Rscript jaccard.R
