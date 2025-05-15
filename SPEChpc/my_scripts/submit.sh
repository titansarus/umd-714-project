#!/bin/bash
#SBATCH --job-name=spechpc_tiny_sweep
#SBATCH --output=output_spechpc_%j.out
#SBATCH --error=error_spechpc_%j.err
#SBATCH -A cmsc714-class
#SBATCH --time=01:30:00
#SBATCH --nodes=4
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=64
#SBATCH --exclusive

#export OMP_PROCESSOR_BIND=true

# Load necessary modules or environment
module purge
module load intel/2021.4.0
module load intel-mpi
module load intel-vtune
#module load openmpi
#module load hpctoolkit

cd /home/jumeike/scratch.cmsc714/final-project/

source shrc

#cd benchspec/HPC/605.lbm_s/run/run_base_ref_gnu_mpi.0010
#cd benchspec/HPC/613.soma_s/run/run_base_ref_gnu_mpi.0010
cd benchspec/HPC/621.miniswp_s/run/run_base_ref_gnu_mpi.0011

# Rank sizes to test
for RANKS in 256 # 32 # 16 8 4 2 1
do
  # Calculate number of nodes assuming 32 tasks per node
  NODES=$(( (RANKS + 63) / 64 ))

  echo "Running SPEChpc with $RANKS ranks on $NODES nodes..."

  #srun --nodes=$NODES --ntasks=$RANKS runhpc --config=Project-try1 --ranks=$RANKS tiny
  #srun --nodes=$NODES runhpc --config=Project-try1 --ranks=$RANKS tiny
  #runhpc --config=Project-try1 --ranks=$RANKS tiny
  #mpirun -np $RANKS sph_exa -n 256 -s 80 -w -1 > sph_exa.out
  #mpirun  -np 256 vtune -collect hotspots -trace-mpi -r 605-vtune-hotspot-result-openmpi  ../run_base_ref_gnu_mpi.0010/lbm_base.gnu_open_mpi_debug 0<&- > lbm_vtune_hotspot_open_mpi.out 2>> lbm_vtune_hotspot_open_mpi.err 
  
  #mpirun  -np 256 vtune -collect threading -trace-mpi -r 605-vtune-threading-result-openmpi  ../run_base_ref_gnu_mpi.0010/lbm_base.gnu_open_mpi_mpi_debug 0<&- > lbm_vtune_threading_open_mpi.out 2>> lbm_vtune_threading_open_mpi.err 
  
  #mpirun  -np 256 vtune -collect hotspots -trace-mpi -r 613-vtune-hotspot-result-openmpi ../run_base_ref_gnu_mpi.0010/soma_base.gnu_open_mpi_debug -r 42 -t 400 --npoly=25000000 --gen-state-file 0<&- > soma_vtune_hotspot_open_mpi.out 2>> soma_vtune_hotspot_open_mpi.err
  
  #mpirun  -np 256 vtune -collect threading -trace-mpi -r 613-vtune-threading-result-openmpi ../run_base_ref_gnu_mpi.0010/soma_base.gnu_open_mpi_debug -r 42 -t 400 --npoly=25000000 --gen-state-file 0<&- > soma_vtune_threading_open_mpi.out 2>> soma_vtune_threading_open_mpi.err
  
  #mpirun  -np 256 hpcrun -e CYCLES -e INSTRUCTIONS -e CPU-CLOCK -e TASK-CLOCK e PERF_COUNT_HW_STALLED_CYCLES_FRONTEND -e STALLED-CYCLES-FRONTEND -e IDLE-CYCLES-FRONTEND -e PERF_COUNT_HW_STALLED_CYCLES_BACKEND -e STALLED-CYCLES-BACKEND -e IDLE-CYCLES-BACKEND -e PAGE-FAULTS -e CONTEXT-SWITCHES -e CPU-MIGRATIONS -e CACHE_MISSES -e L1-DCACHE-STORE-MISSES -e L1-DCACHE-LOAD-MISSES -e L1-DCACHE-PREFETCH-MISSES -e L1-ICACHE-LOAD-MISSES -e L1-ICACHE-PREFETCH-MISSES -e LLC-LOAD-MISSES -e LLC-STORE-MISSES -e LLC-PREFETCH-MISSES -e DTLB-LOAD-MISSES -e  DTLB-PREFETCH-MISSES -e ITLB-LOAD-MISSES ../run_base_ref_gnu_mpi.0011/sweep_base.gnu_open_mpi --niterations 80 --ncell_x 128 --ncell_y 64 --ncell_z 64 --ne 64 --na 32 --nblock_z 8 --nthread_e 1 0<&- > sweep_prof_mpi_1.out 2>> sweep_prof_mpi_1.err
  
  #mpirun -np 256 vtune -collect hotspots -trace-mpi -r 621-vtune-hotspot-result-openmpi ../run_base_ref_gnu_mpi.0011/sweep_base.gnu_open_mpi_debug --niterations 80 --ncell_x 128 --ncell_y 64 --ncell_z 64 --ne 64 --na 32 --nblock_z 8 --nthread_e 1 0<&- > sweep_vtune_hotspot_open_mpi.out 2>> sweep_vtune_hotspot_open_mpi.err
  
  #mpirun -np 256 vtune -collect threading -trace-mpi -r 621-vtune-hotspot-result-openmpi ../run_base_ref_gnu_mpi.0011/sweep_base.gnu_open_mpi_debug --niterations 80 --ncell_x 128 --ncell_y 64 --ncell_z 64 --ne 64 --na 32 --nblock_z 8 --nthread_e 1 0<&- > sweep_vtune_threading_open_mpi.out 2>> sweep_vtune_threading_open_mpi.err
  echo "Completed run with $RANKS ranks"
  echo "--------------------------------"
done

