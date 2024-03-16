#!/bin/bash
#SBATCH --job-name=matlab        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=8G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=96:01:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=kschen@princeton.edu

cd ~/github/Chemotaxis_function

module purge
module load matlab/R2020b

matlab -singleCompThread -nodisplay -nosplash -r "scan_staPAW"
