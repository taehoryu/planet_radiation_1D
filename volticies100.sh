#!/bin/bash
#PBS -l nodes=2:ppn=24,walltime=1:00:00
#PBS -N HJ_rotation
#PBS -q debug

#module load shared
#module load  mpiexec/0.84_432
#module load gcc/5.1.0
#module load anaconda/2
#module load openmpi/open64/64/1.8.5
# module load mvapich2/mlnx/gcc/64/2.0

#activate my_env
module load shared
#module load  mpiexec/0.84_432                                                                                       
module load gcc/7.1.0
module load anaconda/3
#module load openmpi/open64/64/1.8.5                                                                                 
# module load mvapich2/mlnx/gcc/64/2.0                                                                               
module load mvapich2/gcc/64/2.2rc1

module load maui/3.3.1


cd /gpfs/scratch/taryu/project7_HJ/Castro/Exec/science/planet_rotation/
#make realclean
#make clean
#make -j12
mpiexec -n 48 ./Castro3d.gnu.MPI.ex inputs_3d   1>stdout.txt 2>stderr.txt

echo End Job
