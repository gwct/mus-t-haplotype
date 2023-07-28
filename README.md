## Mus t-haplotype phylogenies

This repository contains all scripts and data files associated with the Mus t-haplotype project

### Workflow outline

Briefly, the workflow was:

1. Download genomes
2. Extract t-haplotype homologous chromosomes 
3. Mask chromosomes
4. Align chromosomes with Cactus
5. Convert to MAF and FASTA alignments
6. Partition alignments in windows across the chromosome (10kb)
7. Filter alignments
8. Make window trees with IQ-tree
9. Count topologies

For more detailed, step-by-step instructions see `notes.txt`.

### Files

| Folder | Description | 
| ------ | ----------- |
| `analysis/` | All raw input, intermediate, and raw output files from all steps of the workflow -- **NOTE: This folder and all its sub-folders are excluded from the github repo due to size limitations** |
| `analysis/00-genomes/` | The raw and masked FASTA sequences for the selected genomes |
| `analysis/01-cactus-all-mask/` | The output sequences and HAL file from Cactus |
| `analysis/02-mus-t-windows/` | The input, intermediate, and raw output files for the windowed phylogeny analyses |
| `analysis/02-mus-t-windows/00-maf-rmdups/` | The Cactus alignment in MAF format (converted from HAL), with mm10 as the reference |
| `analysis/02-mus-t-windows/01-bed-windows/` | Bed files that partition the input chromosome (chr17) into 10kb windows |
| `analysis/02-mus-t-windows/02-fasta/` | FASTA files for the alignments of each 10kb window, extracted from the MAF |
| `analysis/02-mus-t-windows/03-fasta-no-anc-filter/` | FASTA files for the alignments of each 10kb window with the ancestral sequences removed and gappy sites filtered |
| `analysis/02-mus-t-windows/04-iqtree/` | Raw IQ-tree and ASTRAL output files, including window trees, a concatenated species tree, and an ASTRAL species tree |
| `analysis/02-mus-t-windows/mus-t-haplotype-6.hal` | The final HAL output file from Cactus, copied over from `analysis/01-cactus-all-mask/` |
| `data/` | Processed summary files from the analyses |
| `data/cactus-stats.csv` | Output from halStats, which summarizes the Cactus alignment |
| `data/cactus-summarize.csv` | Output from halSummarizeMutations, which counts inferred mutations in the Cactus alignment |
| `data/genomes.txt` | The tree and paths to the genomes for input to Cactus |
| `data/mm10-10kb-spec-counts.tsv` | The counts of each species in the final 10kb window alignments |
| `data/mm10-10kb-topo-counts.csv` | The counts and inferred topologies for the 10kb window alignments |
| `data/mm10-10kb-window-stats.tsv` | Alignment stats for filtering of the 10kb window alignments |
| `data/mus-t-haplotype-cactus.tre` | The tree used as input for Cactus in Newick format |   
| `data/samples.csv` | Information on the input genomes |
| `docs/` | Scripts for analyzing the data files and generating HTML reports for the web |
| `figs/` | Figures |
| `scripts/` | All scripts and config files used in the workflow |
| `notes.txt` | Detailed, step-by-step instructions for the workflow |
| `README.md` | This file! |