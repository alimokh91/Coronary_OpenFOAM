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
# NOTE: confirm the exact v2506 module name on your cluster with `module avail OpenFOAM`
#       (the toolchain suffix, e.g. -foss-2022b / -foss-2023a, may differ for v2506)
module load OpenFOAM/v2506-foss-2022b
source /storage/homefs/am20n186/Software/epyc2.9/easybuild/software/OpenFOAM/v2506-foss-2022b/OpenFOAM-v2506/etc/bashrc
# Set OpenMP threads (not used here, but good practice)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# IMPORTANT (v2506 upgrade): rebuild the custom BC library against v2506 once on the
# cluster before running, e.g.:
#   ( cd lumpedParameterNetworkBC/lumpedParameterNetworkBC && wclean libso && wmake libso )
# Also ensure the mesh is generated/decomposed (blockMesh, snappyHexMesh, decomposePar)
# into 64 subdomains before this solver step.

# Run the solver in parallel
srun -n 64 pimpleFoam -parallel | tee log.pimpleFoam

