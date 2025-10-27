#!/bin/bash
#SBATCH --partition=gpu                  # Define the partition on which the job shall run
#SBATCH --job-name=2jpz_POH                  # Job name
#SBATCH --ntasks=1                       # Total number of tasks
#SBATCH --nodes=1                        # Number of nodes on which to run
#SBATCH --cpus-per-task=20               # Threads per task
#SBATCH --gres=gpu:2                     # number of GPUs per node (gres=gpu:N)
#SBATCH --mem=50000                      # Memory per job in MB
#SBATCH --time=0-999:00:00               # Time limit ( dd-hh:mm:ss ) ← e.g. Maximum run time: 1 hour
#SBATCH --time=0-999:00:00               # Time limit ( dd-hh:mm:ss ) ← e.g. Maximum run time: 1 hour
#SBATCH --output=slurm_simulation_%j.out       # Standard output log
#SBATCH --error=slurm_simulation_%j.out        # Standard error log

  echo "Start at `date`"
  echo "Running on hosts: $SLURM_NODELIST"
  echo "Running on $SLURM_NNODES nodes."
  echo "Running $SLURM_NTASKS tasks."
  echo "Using ${SLURM_CPUS_PER_TASK} cpus per task"
  echo "Job id is $SLURM_JOBID"
  echo "Current working directory is `pwd`"

mdpFolder="/home/user/workfolder/g4dna_protocol/smsl/opt/mdp"
mdFolder="."

## load program
export PATH=$HOME/opt/gcc-10.4.0/bin:$PATH
export LD_LIBRARY_PATH=$HOME/opt/gcc-10.4.0/lib64
source $HOME/opt/gromacs-2023/bin/GMXRC

cd ${mdFolder}
## HEAT
gmx grompp -f ${mdpFolder}/heat.mdp -c min2.pdb -r 2jpz_POH_solv.pdb -p topol.top -o heat.tpr
gmx mdrun -deffnm heat -c heat.pdb -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0011 -ntomp 6
## NVT
gmx grompp -f ${mdpFolder}/nvt1.mdp -c heat.pdb -r 2jpz_POH_solv.pdb -t heat.cpt -p topol.top -o nvt1.tpr
gmx mdrun -deffnm nvt1 -c nvt1.pdb -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0011 -ntomp 6
## NPT1
gmx grompp -f ${mdpFolder}/npt1.mdp -c nvt1.pdb -r 2jpz_POH_solv.pdb -t nvt1.cpt -p topol.top -o npt1.tpr
gmx mdrun -deffnm npt1 -c npt1.pdb -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0011 -ntomp 6
## NPT2
gmx grompp -f ${mdpFolder}/npt2.mdp -c npt1.pdb -r 2jpz_POH_solv.pdb -t npt1.cpt -p topol.top -o npt2.tpr
gmx mdrun -deffnm npt2 -c npt2.pdb -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0011 -ntomp 6
