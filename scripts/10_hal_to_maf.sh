#!/bin/bash
#SBATCH --job-name=mus-hal2maf
#SBATCH --output=slurm_logs/%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=greggwct@gmail.com
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=400g
#SBATCH --time=48:00:00

PROJECT_ROOT="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype"

hal="$PROJECT_ROOT/analysis/02-mus-t-windows-new-tree/mus-t-haplotype-6.hal"
ref="mm10"
maf="$PROJECT_ROOT/analysis/02-mus-t-windows-new-tree/00-maf-rmdups/mus-t-haplotype"
cactus_path="/n/holylfs05/LABS/informatics/Users/gthomas/turtles/cactus/cactus_v2.2.0.sif"
outdir="$PROJECT_ROOT/analysis/02-mus-t-windows-new-tree/00-maf-rmdups/"
tmpdir="/n/holylfs05/LABS/informatics/Users/gthomas/tmp/"

cd $outdir

singularity exec --nv --cleanenv --bind $tmp:/tmp $cactus_path hal2mafMP.py --numProc 46 --splitBySequence --noDupes --refGenome $ref $hal $maf
