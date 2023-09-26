#!/bin/bash

#SBATCH --job-name=dns
#SBATCH --output=dns.out

#SBATCH --partition=shortq7
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=2

#SBATCH --time=06:00:00
#SBATCH --exclusive

module load suite-sparse-5.3.0-gcc-8.3.0-no2so7t
module load openblas/0.3.7
module load gsl-2.5-gcc-8.3.0-4sokk4l

#make 3-config options=../dns.cfg THORNLIST=../dns.th
make -j12 dns
