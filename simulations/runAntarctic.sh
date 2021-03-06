#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=03:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2048 # megabytes 
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
cd ../geometry/antarctic/
gmsh -2 -order 1 ant.geo -o ant.msh
cd ../../

# Run the simulation
clear
srun ./build/bin/main ./geometry/antarctic/ant.msh ./params/antarctic.dat ./simulations/resultsAntarctic.msh

# Get back to the initial repository
cd ./simulations
