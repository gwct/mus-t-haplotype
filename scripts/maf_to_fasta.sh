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

wd="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/windows/phast_scripts/"
maf="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/00-maf-rmdups/mus-t-haplotype_chr17.sorted.maf"
bed="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/windows/mm10-17-masked_10kb-windows-named.bed"
outdir="/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/windows/fasta/chr17/"
windowsize="10kb"
ref="mm10"

cd $wd

python maf2fasta.py --maf $maf --bed $bed --ref_species $ref --out_folder $outdir

