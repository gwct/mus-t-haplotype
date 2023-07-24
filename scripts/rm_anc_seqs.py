#!/usr/bin/python
############################################################
# For fasta files from cactus hal > maf files, 07.2023
# Removes the ancestral sequences from a cactus alignment in
# FASTA files that have been converted from the MAF file
############################################################

import sys, os, argparse, core, seqparse as SEQ

############################################################

aln_dir = "/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/windows/fasta/chr17/";
outdir = "/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/windows/fasta-no-anc/chr17/";
if not os.path.isdir(outdir):
    os.makedirs(outdir);

print("start");
sw = 0;
aln_list = [ f for f in os.listdir(aln_dir) if f.endswith(".fa") ];
for aln in aln_list:
    aln_file = os.path.join(aln_dir, aln);
    out_file = os.path.join(outdir, aln.replace(".fa", ".noanc.fa"));

    cur_seqs = SEQ.fastaGetDict(aln_file);
    out_seqs = { header : cur_seqs[header] for header in cur_seqs if "Anc" not in header };

    with open(out_file, "w") as outfile:
        for header in out_seqs:
            outfile.write(header + "\n");
            outfile.write(out_seqs[header] + "\n");
            sw += 1;
print(sw);