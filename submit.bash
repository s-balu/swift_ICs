#!/bin/bash
#PBS -N 2LPTic_L20N256
#PBS -P iv23
#PBS -j oe

#PBS -l ncpus=4
#PBS -l mem=48gb

#PBS -l walltime=00:05:00

#PBS -l storage=gdata/iv23+scratch/iv23

module load openmpi/4.1.2
module load fftw3/3.3.8 hdf5/1.12.1p gsl/2.6

ulimit -c unlimited
cd /scratch/iv23/balu/ICs 

# Run 2LPTic
#mpirunscript="mpirun --mca pml ob1 -x UCX_TLS=ud_x,shm,self -np $PBS_NCPUS"
mpirunscript="mpirun -np $PBS_NCPUS"
executable=./2LPTic
param=L20_N256.param 
#lfs setstripe -c 25 /g/data/iv23/ICs/L35_N2400_5.hdf5
echo ${mpirunscript} ${executable} ${param}
${mpirunscript} ${executable} ${param} > output.log

