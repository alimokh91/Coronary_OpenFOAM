#!/bin/bash
#------------------------------------------------------------------------------
# Local run script for the Coronary case on OpenFOAM v2506 (this WSL machine).
# Meshes (blockMesh + snappyHexMesh in parallel), then runs pimpleFoam in parallel.
# Usage:
#   ./Alllocal.sh          # full pipeline: mesh + solve
#   ./Alllocal.sh mesh     # only generate the mesh
#   ./Alllocal.sh solve    # only run the solver (mesh must already exist)
#------------------------------------------------------------------------------
cd "$(dirname "$0")"

# Load local OpenFOAM v2506 environment (before `set -e`: its bashrc can emit a
# harmless pop_var_context warning that would otherwise abort the script)
source "$HOME/openfoam/OpenFOAM-v2506/etc/bashrc"

set -e

NP=24                                   # MPI ranks for local run (28-core box)
DECDICT="-decomposeParDict system/decomposeParDict.local"
STAGE="${1:-all}"

do_mesh()
{
    echo "==> blockMesh"
    blockMesh > log.blockMesh 2>&1

    echo "==> (re)build custom BC library against v2506"
    ( cd lumpedParameterNetworkBC/lumpedParameterNetworkBC && wclean libso && wmake libso ) \
        > log.wmake.lib 2>&1

    # The 0/ fields reference patches (inlet, outlet, artery, Lc*, Rc*) that only
    # exist AFTER snappyHexMesh. Move them aside so decomposePar/snappy don't try to
    # read fields against the bare blockMesh (a/b/c patches), then restore.
    echo "==> hide 0/ during meshing"
    rm -rf .0_hold && mv 0 .0_hold

    echo "==> decomposePar (mesh only, for parallel snappyHexMesh)"
    decomposePar $DECDICT -force > log.decomposePar.mesh 2>&1

    echo "==> snappyHexMesh -parallel -overwrite on $NP ranks"
    mpirun -np $NP snappyHexMesh $DECDICT -parallel -overwrite > log.snappyHexMesh 2>&1

    echo "==> reconstructParMesh"
    reconstructParMesh -constant > log.reconstructParMesh 2>&1

    echo "==> restore 0/"
    mv .0_hold 0

    echo "==> checkMesh"
    checkMesh -constant > log.checkMesh 2>&1 || true

    # clean processor dirs so the solver can re-decompose cleanly
    rm -rf processor*
}

do_solve()
{
    echo "==> decomposePar (for solver)"
    decomposePar $DECDICT -force > log.decomposePar 2>&1

    echo "==> pimpleFoam -parallel on $NP ranks"
    mpirun -np $NP pimpleFoam $DECDICT -parallel | tee log.pimpleFoam
}

case "$STAGE" in
    mesh)  do_mesh ;;
    solve) do_solve ;;
    all)   do_mesh; do_solve ;;
    *) echo "unknown stage: $STAGE (use: mesh | solve | all)"; exit 1 ;;
esac

echo "==> done ($STAGE)"
