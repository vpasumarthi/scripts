#!/usr/bin/env python

from pathlib import Path

import numpy as np
import yaml


class SimulationFiles(object):
    """Class definition to generate simulation files"""

    # constants
    MIN2SEC = 60
    HR2SEC = 60 * MIN2SEC

    def __init__(self, sim_param_file_name, var_species_type_index,
                 var_species_count_list, kmc_prec):
        # Load simulation parameters
        with open(sim_param_file_name, 'r') as stream:
            try:
                params = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        for key, value in params.items():
            setattr(self, key, value)
        self.var_species_type_index = var_species_type_index
        self.non_var_species_type_index = int(not var_species_type_index)
        self.var_species_count_list = var_species_count_list
        self.n_runs = len(var_species_count_list)
        self.kmc_prec = kmc_prec

    def generate_files(self, gen_params_files, gen_run_files, gen_msd_files,
                       gen_slurm_files):
        if gen_params_files:
            self.simulation_parameter_files()
        if gen_run_files:
            self.run_files()
        if gen_msd_files:
            self.msd_files()
        if gen_slurm_files:
            self.slurm_files(self.kmc_prec)
        return None

    def dst_path(self, species_count_list):
        # determine destination path
        child_dir1 = 'SimulationFiles'
        child_dir2 = ('ion_charge_type=' + self.system['ion_charge_type']
                      + ';species_charge_type='
                      + self.system['species_charge_type'])
        child_dir3 = (
                str(species_count_list[0])
                + ('electron' if species_count_list[0] == 1 else 'electrons')
                + ',' + str(species_count_list[1])
                + ('hole' if species_count_list[1] == 1 else 'holes'))
        child_dir4 = str(self.system['temp']) + 'K'
        child_dir5 = (('%1.2E' % self.run['t_final']) + 'SEC,'
                      + ('%1.2E' % self.run['time_interval']) + 'TimeInterval,'
                      + ('%1.2E' % self.run['n_traj']) + 'Traj')
        electric = self.system['external_field']['electric']
        if electric['active']:
            ld_tag = 'ld_' if electric['ld'] else ''
            self.field_tag = ('ef_' + ld_tag +
                              str(electric['dir']).replace(' ', '') + '_'
                              + ('%1.2E' % electric['mag']))
        else:
            self.field_tag = 'no_field'
        work_dir = self.field_tag
        system_directory_path = Path.cwd()
        work_dir_path = (system_directory_path / child_dir1 / child_dir2
                         / child_dir3 / child_dir4 / child_dir5 / work_dir)
        work_dir_depth = (len(work_dir_path.parts)
                          - len(system_directory_path.parts))
        return (work_dir_path, work_dir_depth)

    def simulation_parameter_files(self):
        species_count_list = [0] * len(self.system['species_count'])
        for i_run in range(self.n_runs):
            species_count_list[self.non_var_species_type_index] = (
                self.system['species_count'][self.non_var_species_type_index])
            species_count_list[self.var_species_type_index] = (
                                            self.var_species_count_list[i_run])

            (work_dir_path, work_dir_depth) = self.dst_path(species_count_list)
            self.system['work_dir_depth'] = work_dir_depth
            self.system['species_count'] = species_count_list
            Path.mkdir(work_dir_path, parents=True, exist_ok=True)
            dst_file_path = work_dir_path.joinpath(self.system['dst_file_name']
                                                   )
            # generate simulation parameter file
            with dst_file_path.open('w') as dst_file:
                dst_file.write('# System parameters:\n')
                yaml.dump(self.system, dst_file)
                dst_file.write('\n')
                dst_file.write('# Run parameters:\n')
                yaml.dump(self.run, dst_file, default_flow_style=False)
                dst_file.write('\n')
                dst_file.write('# MSD parameters:\n')
                yaml.dump(self.msd, dst_file, default_flow_style=False)
        return None

    def run_files(self):
        species_count_list = [0] * len(self.system['species_count'])
        for i_run in range(self.n_runs):
            species_count_list[self.non_var_species_type_index] = (
                self.system['species_count'][self.non_var_species_type_index])
            species_count_list[self.var_species_type_index] = (
                                            self.var_species_count_list[i_run])
            (work_dir_path, _) = self.dst_path(species_count_list)
            dst_file_path = work_dir_path.joinpath(self.run['dst_file_name'])

            # generate simulation parameter file
            with dst_file_path.open('w') as dst_file:
                dst_file.write(
                    "#!/usr/bin/env python\n\n"
                    "from pathlib import Path\n\n"
                    "from PyCT.material_run import material_run\n\n"
                    "dst_path = Path.cwd()\n"
                    "material_run(dst_path)\n")
            dst_file_path.chmod(0o755)
        return None

    def msd_files(self):
        species_count_list = [0] * len(self.system['species_count'])
        for i_run in range(self.n_runs):
            species_count_list[self.non_var_species_type_index] = (
                self.system['species_count'][self.non_var_species_type_index])
            species_count_list[self.var_species_type_index] = (
                                            self.var_species_count_list[i_run])
            (work_dir_path, _) = self.dst_path(species_count_list)
            dst_file_path = work_dir_path.joinpath(self.msd['dst_file_name'])

            # generate simulation parameter file
            with dst_file_path.open('w') as dst_file:
                dst_file.write(
                    "#!/usr/bin/env python\n\n"
                    "from pathlib import Path\n\n"
                    "from PyCT.material_msd import material_msd\n\n"
                    "dst_path = Path.cwd()\n"
                    "material_msd(dst_path)\n")
            dst_file_path.chmod(0o755)
        return None

    def run_time(self, species_count_list, kmc_prec):
        k_total = np.dot(self.k_total_per_species, species_count_list)
        time_step = 1 / k_total
        kmc_steps = int(np.ceil(self.run['t_final'] / time_step / kmc_prec)
                        * kmc_prec)
        num_states_per_step = np.dot(self.num_states_per_species,
                                     species_count_list)
        total_states_per_traj = num_states_per_step * kmc_steps
        est_run_time = int(self.time_per_state[self.var_species_type_index]
                           * self.run['n_traj'] * total_states_per_traj
                           + (self.add_on_time_limit * self.HR2SEC))
        return est_run_time

    def slurm_files(self, kmc_prec):
        # keywords
        job_name_key = ('--job-name="' + self.system['material'] + '-'
                        + 'x'.join(str(element)
                                   for element in self.system['system_size']))
        output_key = ('--output=' + self.system['material'] + '-'
                      + 'x'.join(str(element)
                                 for element in self.system['system_size']))

        charge_comb = (self.system['ion_charge_type'][0]
                       + self.system['species_charge_type'][0])
        species_count_list = [0] * len(self.system['species_count'])
        species_tag = 'e' if self.var_species_type_index == 0 else 'h'
        for i_run in range(self.n_runs):
            # estimate simulation run time in sec
            species_count_list[self.non_var_species_type_index] = (
                self.system['species_count'][self.non_var_species_type_index])
            species_count_list[self.var_species_type_index] = (
                                            self.var_species_count_list[i_run])
            (work_dir_path, _) = self.dst_path(species_count_list)
            Path.mkdir(work_dir_path, parents=True, exist_ok=True)
            dst_file_path = work_dir_path.joinpath(self.slurm['dst_file_name'])

            # generate slurm file
            with dst_file_path.open('w') as dst_file:
                dst_file.write('#!/bin/sh\n')
                dst_file.write('#SBATCH ' + job_name_key + '_' + charge_comb
                               + '_' + self.field_tag + '_' + species_tag
                               + str(self.var_species_count_list[i_run])
                               + '"\n')
                dst_file.write('#SBATCH ' + output_key + '_' + charge_comb
                               + '_' + self.field_tag + '_' + species_tag
                               + str(self.var_species_count_list[i_run])
                               + '.out\n')
                dst_file.write(
                    f"#SBATCH --partition={self.slurm['partition_value']}\n")
                if self.slurm['partition_value'] == 'mdupuis2':
                    dst_file.write('#SBATCH --clusters=chemistry\n')
                num_days = num_hours = num_mins = num_sec = 0
                if self.slurm['partition_value'] == 'debug':
                    num_hours = 1
                elif self.slurm['partition_value'] == 'mdupuis2':
                    num_days = self.slurm['md_slurm_job_max_time_limit']
                else:
                    est_run_time = self.run_time(species_count_list, kmc_prec)
                    if est_run_time > self.gc_slurm_job_max_time_limit:
                        num_hours = self.gc_slurm_job_max_time_limit
                    else:
                        num_hours = est_run_time // self.HR2SEC
                        num_mins = (est_run_time // self.MIN2SEC) % self.MIN2SEC
                time_limit = f'{num_days:02d}-{num_hours:02d}:{num_mins:02d}:{num_sec:02d}'
                dst_file.write(f'#SBATCH --time={time_limit}\n')
                dst_file.write(f"#SBATCH --nodes={self.slurm['num_nodes']}\n")
                dst_file.write(f"#SBATCH --tasks-per-node={self.slurm['num_tasks_per_node']}\n")
                if self.slurm['exclusive']:
                    dst_file.write('#SBATCH --exclusive\n')
                if self.slurm['mem']:
                    dst_file.write(f"#SBATCH --mem={self.slurm['mem']}\n\n")
                if self.slurm['email']:
                    dst_file.write(f"#SBATCH --mail-user={self.slurm['email']}\n")
                    dst_file.write("#SBATCH --mail-type=END\n")
                dst_file.write(
                    "#SBATCH --constraint=IB\n"
                    "\n# Job description:\n"
                    "# run KMC simulation followed by performing MSD analysis"
                    " of the output trajectories\n\n"
                    "echo \"SLURM_JOBID=\"$SLURM_JOBID\n"
                    "echo \"SLURM_JOB_NODELIST\"=$SLURM_JOB_NODELIST\n"
                    "echo \"SLURM_NNODES\"=$SLURM_NNODES\n"
                    "echo \"SLURMTMPDIR=\"$SLURMTMPDIR\n\n"
                    "echo \"working directory = \"$SLURM_SUBMIT_DIR\n\n"
                    "module load python\n"
                    "source activate py36\n"
                    "module list\n"
                    "ulimit -s unlimited\n\n"
                    "# The initial srun will trigger the SLURM prologue on"
                    " the compute nodes.\n"
                    "NPROCS=`srun --nodes=${SLURM_NNODES}"
                    " bash -c 'hostname' |wc -l`\n"
                    "echo NPROCS=$NPROCS\n"
                    "echo \"Launch mymodel with srun\"\n\n"
                    "#The PMI library is necessary for srun\n"
                    "export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so\n")
                if self.slurm['submit_run']:
                    dst_file.write("srun Run.py\n")
                if self.slurm['submit_m_s_d']:
                    dst_file.write("srun MSD.py\n")
                dst_file.write("\necho \"All Done!\"\n")
        return None
