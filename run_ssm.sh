#!/bin/bash
#SBATCH --job-name=matlab        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=4               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=16G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=96:01:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=kschen@princeton.edu


# Check if the number of command-line arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input1 input2"
    exit 1
fi

# Extract command-line arguments
input1="$1"
input2="$2"

cd ~/github/Chemotaxis_function

module purge
module load matlab/R2020b

# Call the MATLAB script with input variables
matlab -nodisplay -nosplash -r "run_staPAW_func('$input1', '$input2'); exit;"

