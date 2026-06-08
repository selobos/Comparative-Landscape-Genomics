#!/bin/bash
    
# ==========================================
# BayPass Core Model Runs
# Run BayPass under the core model mode to generate covariance matrix
# ==========================================

#Set working directory 
cd 3_GEA/3_1_BayPass/3_1_1_BM/3_1_1_1_Core

#Load BayPass
module load baypass

# Verify BayPass is available
g_baypass

# Run 5 independent replicates
for i in 1 2 3 4 5
do
seed=$((1000 + RANDOM % 9999))
echo "$seed"
nohup g_baypass -npop 30 -gfile BM_Baypass -seed $seed -outprefix bm_core.$i -nval 5000 -burnin 50000 -nthreads 48 &
done

#In R: calculate the means from the initial BayPass core model runs, and produce a covariance matrix using the "BM_Get_median_mat_omega.R" script.
