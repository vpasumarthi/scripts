#!/bin/bash

n_traj=100
sleep_time=2

doping_tag="undoped"
start_carriers=4
end_carriers=9
t_final="1.00E-02"
time_interval="1.00E-08"
n_traj_key="1.00E+02"
execute_pre_prod=1

for ((j=$start_carriers; j<=$end_carriers; j++))
do
	cd "SimulationFiles/"$j"electrons,0holes/"$t_final"SEC,"$time_interval"TimeInterval,"$n_traj_key"Traj/no_field_"$doping_tag
	if [ $execute_pre_prod -eq 1 ]
	then
		./pre_prod.py
	fi
	for ((i=1; i<=$n_traj; i++))
	do
	    cd "traj"$i
	    sbatch slurmscript | awk -vORS=, '{ print $4 >> "../job_ids.txt" }'
	    sleep $sleep_time
	    cd ..
	done

	./generate_slurm_msd.py
	sbatch slurmscript_msd
	cd ../../../../
done
