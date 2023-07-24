#!/bin/bash
#SBATCH --job-name=mus-hal2maf
#SBATCH --output=slurm-logs/%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=greggwct@gmail.com
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=400g
#SBATCH --time=48:00:00

hal="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/mus-t-haplotype-6.hal"
ref="mm10"
maf="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/00-maf-rmdups/mus-t-haplotype"
cactus_path="/n/holylfs05/LABS/informatics/Users/gthomas/turtles/cactus/cactus_v2.2.0.sif"
outdir="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/00-maf-rmdups/"
tmpdir="/n/holylfs05/LABS/informatics/Users/gthomas/tmp/"

cd $outdir

singularity exec --nv --cleanenv --bind $tmp:/tmp $cactus_path hal2mafMP.py --numProc 46 --splitBySequence --noDupes --refGenome $ref $hal $maf