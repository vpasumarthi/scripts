#!/bin/bash

n_traj=100
sleep_time=2

for ((i=1; i<=$n_traj; i++))
do
    cd "traj"$i
    sbatch slurmscript | awk -vORS=, '{ print $4 >> "../job_ids.txt" }'
    sleep $sleep_time
    cd ..
done

./generate_slurm_msd.py
sbatch slurmscript_msd
