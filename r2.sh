#!/bin/bash
#SBATCH --job-name=OpenFOAM_simulation
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=23:00:00
#SBATCH --output=run_%j.out

# Load modules and OpenFOAM environment
module load EasyBuild
module load Workspace_Home
module load OpenFOAM/v2306-foss-2022b
source /storage/homefs/am20n186/Software/epyc2.9/easybuild/software/OpenFOAM/v2306-foss-2022b/OpenFOAM-v2306/etc/bashrc
# Set OpenMP threads (not used here, but good practice)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Preprocessing
# Run snappyHexMesh in parallel
srun -n 64 pimpleFoam -parallel | tee log.pimpleFoam

