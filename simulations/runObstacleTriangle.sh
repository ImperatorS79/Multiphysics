#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=01:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=out.txt

# Load the modules & set the exports
module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16
export OMP_CANCELLATION=true

# Generate the .msh
cd $HOME/Multiphysics
cd ./geometry/obstacle/
gmsh -2 -order 3 triangle.geo -o triangle.msh
cd ../../

# Run the simulation
clear
srun ./build/bin/main ./geometry/obstacle/triangle.msh ./params/obstacleTriangle.dat

# Get back to the initial repository
cd ./simulations
