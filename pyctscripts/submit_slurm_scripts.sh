#!/bin/bash

n_traj=100
sleep_time=2

doping_tag="undoped"
start_carriers=4
end_carriers=9
t_final="1.00E-02"
time_interval="1.00E-08"
n_traj_key="1.00E+02"

for ((i=1; i<=$n_traj; i++))
do
    cd "traj"$i
    sbatch slurmscript | awk -vORS=, '{ print $4 >> "../job_ids.txt" }'
    sleep $sleep_time
    cd ..
done

./generate_slurm_msd.py
sbatch slurmscript_msd
