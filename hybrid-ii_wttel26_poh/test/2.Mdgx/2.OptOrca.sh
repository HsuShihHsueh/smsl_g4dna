scp -r conf1/Conf{1..140}.orca           scale:~/workfolder/gaff/tmpyp4_amber_gaff2/2.Mdgx/conf1
scp -r conf1/Conf{141..280}.orca  multiphysics:~/workfolder/gaff/tmpyp4_amber_gaff2/2.Mdgx/conf1
scp -r conf1/Conf{281..360}.orca        chu_02:~/user/allen/workfolder/gaff/tmpyp4_amber_gaff2/2.Mdgx/conf1


qsub run_orca.qsub

scp        scale:/home/allen/workfolder/gaff/tmpyp4_amber_gaff2/2.Mdgx/conf1/Conf*.oout conf1/
scp multiphysics:/home/allen/workfolder/gaff/tmpyp4_amber_gaff2/2.Mdgx/conf1/Conf*.oout conf1/
scp       chu_02:/home/chu_02/user/allen/workfolder/gaff/tmpyp4_amber_gaff2/2.Mdgx/conf1/Conf*.oout conf1/

bash add_mark.sh

## small script ##

run_orca.qsub
#!/bin/bash
### Declare this job
#PBS -N spOrca
### Using ncpus threads
#PBS -l select=1:ncpus=16
### Combine PBS standard output and error files
#PBS -j oe
### Create an array
#PBS -J 1-140

cd ${PBS_O_WORKDIR}
scratchFile="/scratch/allen/tmpyp4_amber_gaff2"
confFilename="`basename \`realpath .\``"
currentConfFile="${PBS_O_WORKDIR}"
scratchConfFile="${scratchFile}/${confFilename}"


export PATH="$PATH:$HOME/opt/openmpi-4.1.6/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/opt/openmpi-4.1.6/lib"


mkdir -p ${scratchConfFile}
cp -r ${currentConfFile}/Conf${PBS_ARRAY_INDEX}.orca ${scratchConfFile}

sed -i 's/PAL2/PAL16/g' ${scratchConfFile}/Conf${PBS_ARRAY_INDEX}.orca
sed -i 's/MaxCore 512/MaxCore 2000/g' ${scratchConfFile}/Conf${PBS_ARRAY_INDEX}.orca
~/opt/orca-5.0.4/orca ${scratchConfFile}/Conf${PBS_ARRAY_INDEX}.orca > ${scratchConfFile}/Conf${PBS_ARRAY_INDEX}.oout

cp ${scratchConfFile}/Conf${PBS_ARRAY_INDEX}.oout ${currentConfFile}


run_orca.slurm
#!/bin/bash
#SBATCH --partition=xeon                 # Define the partition on which the job shall run
#SBATCH --job-name=spOrca               # Job name
#SBATCH --ntasks=8                       # Total number of tasks
#SBATCH --nodes=1                        # Number of nodes on which to run
#SBATCH --cpus-per-task=1                # Threads per task
#SBATCH --mem=360000                     # Memory per job in MB
#SBATCH --output=spOrca.%J.out
#SBATCH --error=spOrca.%J.out
#SBATCH --array=281-360

scratchFile="/scratch/chu_02/allen/tmpyp4_amber_gaff2"
confFilename=$(basename $(realpath .))
currentConfFile=$(pwd)
scratchConfFile="${scratchFile}/${confFilename}"

module load openmpi/gcc/64

mkdir -p ${scratchConfFile}
cp -r Conf${SLURM_ARRAY_TASK_ID}.orca ${scratchConfFile}

sed -i 's/PAL16/PAL8/g' ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.orca
sed -i 's/MaxCore 512/MaxCore 7000/g' ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.orca
~/opt/orca-5.0.4/orca ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.orca > ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.oout

cp ${scratchConfFile}/Conf${SLURM_ARRAY_TASK_ID}.oout ${currentConfFile}


add_mark.sh
#!/bin/bash
for filename in conf1/Conf*.oout; do
    echo -e "# of primitive gaussian functions\n" >> ${filename}
done


