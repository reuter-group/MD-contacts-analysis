#!/bin/bash

### Example for submitting a large md-analysis job on HPC

#SBATCH --account=nn4700k
#SBATCH --job-name=md-ana
#SBATCH --nodes=1

#SBATCH --time=02:00:00  
#SBATCH --mem-per-cpu=2G  # can change to 4G, depending on number of used cpus
#SBATCH --cpus-per-task=32  # can change to 64, or 128

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors


# choose  working folder:
cd /cluster/projects/nn4700k/md_analysis

# turn on the environment
. /cluster/projects/nn4700k/miniconda3/etc/profile.d/conda.sh
conda activate md-ana-py37

## example test
./cation_pi_ana.py test_data/phosphod_bound_to_pure_popc.old.psf test_data/100frames.dcd test_data/cation_pi_candidates.txt 2  > cation_pi_test_results.txt
