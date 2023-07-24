#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re
from itertools import product

#############################################################################
# Example cmd

# snakemake -p -s make_window_trees.smk --configfile windows-config.yaml --profile profiles/slurm_profile/ --keep-going --dryrun

#############################################################################
# Functions

#############################################################################
# Input and output info

work_dir = os.path.dirname(__file__);
#print("working dir: " + work_dir);

winsize_kb = config["winsize_kb"];
#print("winsize_kb: " + winsize_kb);

#repeat_cutoff = config["repeat_cutoff"];
#print("repeat_cutoff: " + repeat_cutoff);

#missing_cutoff = config["missing_cutoff"];
#print("missing_cutoff: " + missing_cutoff);

OUTGROUP = config["outgroup"];
INPUT_CHROMES = config["chromes"];

#print("reading window ids...");

ALNDIR = config["aln_base_dir"];
TREE_OUTDIR = config["tree_out_dir"];
SUB_DIR = winsize_kb# + "-" + repeat_cutoff + "-" + missing_cutoff;
# Dirs

######################

#WINDOW_INPUT = os.path.join(config["data_dir"], SUB_DIR + "-windows-filter.tsv");
WINDOW_INPUT = "/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/scripts/get-windows-default.tsv";
#print("window file: " + WINDOW_INPUT);

CHROMES, WINDOWS, first = [], [], True;
for line in open(WINDOW_INPUT):
    if line[0] == "#":
        continue;
    if first:
        first = False;
        continue;
    # Skip header lines in the window file

    line = line.strip().split("\t");
    #print(line[0]);
    #print(line[0], float(line[46]), int(line[28]) >= 4 and float(line[46]) < 0.5);
    if line[1] in INPUT_CHROMES and int(line[32]) == 6 and int(line[34]) >= 4 and line[42] == "0":
    # 32: 6 sequences are present in the alignment
    # 34: At least 4 unique sequences are present in the alignment
    # 42: There are no sequences that are comprised of all missing or gap characters
    #if line[0] in INPUT_CHROMES and all( filters == "PASS" for filters in [line[5], line[7], line[8]] ):
        CHROMES.append(line[1]);
        WINDOWS.append(line[0].split(".")[0]);
        #WINDOWS.append(line[3].replace(":", "-"));
## Read the windows to align given the input chromes and the filters

#print(len(WINDOWS), "window ids read");

######################

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        gene_trees_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}.noanc.filter.treefile.rooted"), zip, chrome=CHROMES, window=WINDOWS),
        astral_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile.rooted"), chrome=CHROMES),
        astral_cf_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree.rooted"), chrome=CHROMES),
        concat_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile.rooted"), chrome=CHROMES),
        concat_cf_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree.rooted"), chrome=CHROMES)


        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}.noanc.filter.treefile"), zip, chrome=CHROMES, window=WINDOWS),
        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile"), chrome=CHROMES),
        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "astral.treefile"), chrome=CHROMES),

        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile"), chrome=CHROMES)

#############################################################################
# Pipeline rules

rule iqtree:
    input:
        aln_file = os.path.join(ALNDIR, "{chrome}-filter", "{window}.noanc.filter.fa")
    output:
        tree_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}.noanc.filter.treefile"),
        tree_file_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}.noanc.filter.treefile.rooted"),
        aln_pass_file = os.path.join(ALNDIR, "{chrome}-filter-used", "{window}.noanc.filter.fa")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "logs", "{window}-iqtree.log")
    params:
        bootstrap = "1000",
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}.noanc.filter"),
        outgroup = OUTGROUP
    resources:
        cpus=1
    shell:
        """
        iqtree -s {input.aln_file} --prefix {params.prefix} -B {params.bootstrap} -T {resources.cpus} &> {log}

        nw_reroot {output.tree_file} {params.outgroup} > {output.tree_file_rooted}

        cp {input.aln_file} {output.aln_pass_file}
        """
# Run each locus through iqtree

######################

rule combine_gt:
    input:
        tree_file = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}.noanc.filter.treefile"), zip, chrome=CHROMES, window=WINDOWS)
    output:
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile")
    params:
        loci_dir = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci")
    resources:
        cpus=12
    shell:
        """
        find {params.loci_dir} -name *.treefile -exec cat {{}} \\; > {output.gene_trees_file}
        """
# Combine the gene trees into a single file

######################

rule astral:
    input:
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile")
    output:
        astral_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile"),
        astral_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "astral.log")
    params:
        outgroup = OUTGROUP
    resources:
        mem="24g"
    shell:
        """
        astral -i {input.gene_trees_file} -o {output.astral_tree} 2> {log}

        nw_reroot {output.astral_tree} {params.outgroup} > {output.astral_tree_rooted}
        """
# Infer species tree with ASTRAL

######################

rule astral_cf:
    input:
        aln_file = expand(os.path.join(ALNDIR, "{{chrome}}-filter-used", "{{window}}.noanc.filter.fa"), zip, chrome=CHROMES, window=WINDOWS),
        aln_pass_dir = os.path.join(ALNDIR, "{chrome}-filter-used"),    
        astral_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile"),
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile"),
    output:
        concord_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree"),
        concord_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral-iqtree-cf.log")
    params:
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral"),
        outgroup = OUTGROUP
    resources:
        cpus=1,
        mem="24g"
    shell:
        """
        iqtree -t {input.astral_tree} --gcf {input.gene_trees_file} -p {input.aln_pass_dir} --scf 100 -T {resources.cpus} --prefix {params.prefix} &> {log}

        nw_reroot {output.concord_tree} {params.outgroup} > {output.concord_tree_rooted}
        """
# Calculate concordance factors on the astral tree

######################

rule concat:
    input:
        aln_file = expand(os.path.join(ALNDIR, "{{chrome}}-filter-used", "{{window}}.noanc.filter.fa"), zip, chrome=CHROMES, window=WINDOWS),
        aln_pass_dir = os.path.join(ALNDIR, "{chrome}-filter-used")
    output:
        concat_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile"),
        concat_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat-iqtree.log")
    params:
        bootstrap = "1000",
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat"),
        outgroup = OUTGROUP
    resources:
        cpus=12,
        time="12:00:00",
        mem="24g"
    shell:
        """
        iqtree -redo -s {input.aln_pass_dir} --prefix {params.prefix} -B {params.bootstrap} -T {resources.cpus} &> {log}

        nw_reroot {output.concat_tree} {params.outgroup} > {output.concat_tree_rooted}
        """
# Run each locus through iqtree with concatenation

######################

rule concat_cf:
    input:
        aln_file = expand(os.path.join(ALNDIR, "{{chrome}}-filter-used", "{{window}}.noanc.filter.fa"), zip, chrome=CHROMES, window=WINDOWS),
        aln_pass_dir = os.path.join(ALNDIR, "{chrome}-filter-used"),    
        concat_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile"),
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile"),
    output:
        concord_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree"),
        concord_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat-iqtree-cf.log")
    params:
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat"),
        outgroup = OUTGROUP
    resources:
        cpus=1,
        mem="24g"
    shell:
        """
        iqtree -redo -t {input.concat_tree} --gcf {input.gene_trees_file} -p {input.aln_pass_dir} --scf 100 -T {resources.cpus} --prefix {params.prefix} &> {log}
        # -redo needed because the checkpoint from the main tree search reports that the run has already been done

        nw_reroot {output.concord_tree} {params.outgroup} > {output.concord_tree_rooted}
        """
# Calculate concordance factors on the concatenated tree

######################

