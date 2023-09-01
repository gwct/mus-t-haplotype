#!/bin/bash
#SBATCH --job-name=mus_maf2fasta
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=greggwct@gmail.com
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g
#SBATCH --time=24:00:00

PROJECT_ROOT="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype"

wd="$PROJECT_ROOT/scripts/phast_scripts/"
maf="$PROJECT_ROOT/analysis/02-mus-t-windows-new-tree/00-maf-rmdups/mus-t-haplotype_chr17.sorted.maf"
bed="$PROJECT_ROOT/analysis/02-mus-t-windows-new-tree/01-bed-windows/mm10-17-masked_10kb-windows-named.bed"
outdir="$PROJECT_ROOT/analysis/02-mus-t-windows-new-tree/02-fasta/chr17/"
windowsize="10kb"
ref="mm10"

cd $wd
mkdir $outdir

python maf2fasta.py --maf $maf --bed $bed --ref_species $ref --out_folder $outdir

