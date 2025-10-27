#!/bin/bash
#SBATCH --partition=xeon                 # Define the partition on which the job shall run
#SBATCH --job-name=spOrca                # Job name
#SBATCH --ntasks=8                       # Total number of tasks
#SBATCH --nodes=1                        # Number of nodes on which to run
#SBATCH --cpus-per-task=1                # Threads per task
#SBATCH --mem=8000                       # Memory per job in MB
#SBATCH --output=spOrca.%J.out
#SBATCH --error=spOrca.%J.out
#SBATCH --array=1-360

scratchFile="/scratch/chu_02/allen/amber_gaff2"
confFilename=$(basename $(realpath .))
currentConfFile=$(pwd)
scratchConfFile="${scratchFile}/${confFilename}"
orca="$HOME/opt/orca-5.0.4/orca"

module load openmpi/gcc/64

mkdir -p ${scratchConfFile}
cp -r conf1/Conf${SLURM_ARRAY_TASK_ID}.orca ${scratchConfFile}

${orca} ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.orca > ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.oout

cp ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.oout ${currentConfFile}/conf1
echo -e "# of primitive gaussian functions\n" >> ${currentConfFile} ## add mark to let 3.ExtractEnergies.mdgx can reconginze
