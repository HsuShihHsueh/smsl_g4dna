#!/bin/bash
#SBATCH --partition=gpu                  # Define the partition on which the job shall run
#SBATCH --job-name=hybrid-ii          # Job name
#SBATCH --ntasks=1                       # Total number of tasks
#SBATCH --nodes=1                        # Number of nodes on which to run
#SBATCH --cpus-per-task=20               # Threads per task
#SBATCH --gres=gpu:2                     # number of GPUs per node (gres=gpu:N)
#SBATCH --mem=50000                      # Memory per job in MB
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

  # Unload all loaded module and reset everything to original state.
  nstart=0
  nend=9 
  nstep=$((500*500000)) #500ns 500000
  nstep_ns=$((2*${nstep}/1000000)) 
 
  user="allen"
  jobname="2jpz"
  data_dir=`pwd`
  scratch_dir="/scratch/chu_02/${user}/mdrun/hybrid-ii/2.simulation"
  mdp_folder="/home/user/workfolder/g4dna_protocol/smsl/opt/mdp"
  echo "data dircetory: ${data_dir}"
  echo -n "scratch_directory: ${scratch_dir} ..."
  echo "mdp dircetory: ${mdp_folder}"


  if [ -d "${scratch_dir}" ]; then
      echo "exists!"
  else
      echo "mkdir -p ${scratch_dir}"
      mkdir -p ${scratch_dir}
  fi
  
  echo "moving data to scratch ..."
  cp -r ${data_dir} ${scratch_dir}/..

  ## load program
  export PATH=$HOME/opt/gcc-10.4.0/bin:$PATH
  export LD_LIBRARY_PATH=$HOME/opt/gcc-10.4.0/lib64
  source $HOME/opt/gromacs-2023/bin/GMXRC

  cd ${scratch_dir}
  prename="npt2" 
  for nrun in $(seq ${nstart} ${nend})
  do
      ninit=$((${nrun}*${nstep}))
      curname=md${nrun}
      
      echo "#################### nrun: ${nrun} ##########################"
      echo "curname  is: ${curname}"
      echo "prename is: ${prename}" 

      echo "; set initial step for trajectory continuation" > ${curname}.mdp
      echo "init-step             = ${ninit}" >> ${curname}.mdp
      echo "nsteps                = ${nstep} ; 2 * ${nstep} = ${nstep_ns} ns" >> ${curname}.mdp
      cat ${mdp_folder}/template.mdp  >> ${curname}.mdp

      gmx grompp -f ${curname}.mdp -po ${curname}.out.mdp -c ${prename}.pdb -t ${prename}.cpt -p topol.top -o ${curname}.tpr
      gmx mdrun -deffnm ${curname} -c ${curname}.pdb -ntmpi 4 -nb gpu -pme gpu -npme 1 -gputasks 0011 -ntomp 6

      echo "moving ${curname} from scratch to data ..."
      cp ${scratch_dir}/${curname}* ${data_dir}/
      prename=${curname}

      ssh frontend << EOF
        mail -s "${curname}: job \"${SLURM_JOB_NAME}\" is completed" \
        shhsu.bt11@nycu.edu.tw <<< "${curname}: job \"${SLURM_JOB_NAME}\"(${SLURM_JOBID}) is completed at `date`"
EOF

  done   
  echo "End at `date`"

  ssh frontend << EOF
    mail -s "job \"${SLURM_JOB_NAME}\"(${SLURM_JOBID}) is completed" \
    shhsu.bt11@nycu.edu.tw <<< "job \"${SLURM_JOB_NAME}\"(${SLURM_JOBID}) is completed at `date`"
EOF

