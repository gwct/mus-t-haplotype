#############################################################################
# Samples quartets of tips to calculate D-statistics
#
# Gregg Thomas, November 2022
#############################################################################

import sys, os, random
from itertools import combinations
import core
import treec as tp
import seqparse as seq
import multiprocessing as mp

#############################################################################

def getSubAln(spec_list, aln, aln_skip_chars):

    orig_aln_len = len(aln[spec_list[0]]);
    # Get the alignment length from the first sequence

    sub_aln = { spec : list(aln[spec]) for spec in spec_list };
    # Get the sub-alignment for the current list of species

    rm_pos = [ j for j in range(orig_aln_len) if any(sub_aln[spec][j] in aln_skip_chars for spec in spec_list) ];
    # Get a list of positions that contain missing data in any of the sample species

    final_sub_aln = { spec : "" for spec in spec_list };
    for spec in spec_list:
        for j in range(orig_aln_len):
            if j not in rm_pos:
                final_sub_aln[spec] += sub_aln[spec][j];
    # Remove sites with missing data in any sample from all species here by concatenating the bases together at all positions

    sub_aln_len = len(final_sub_aln[spec_list[0]]);
    # Adjusted alignment length after removing columns with missing data/gaps

    return final_sub_aln, sub_aln_len;

####################

def getHamming(seq1, seq2):
# Calculates hamming distance between two strings (sequences) 
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

####################

def pairReps(tip_pair_info):
# This function takes a set of alignments and set of species and calculates 
# dxy for each gene in the alignment for the given tip pair as well as dxy
# relative to the outgroup species
   
    tip_pair, outgroup, alns, aln_skip_chars = tip_pair_info;
    # Unpack the arguments passed to the function

    outlines = [];
    # Each gene for each quartet will result in site counts and we save the output
    # in outlines

    for locus in alns:
        #print(locus);
        aln = alns[locus];
        # Lookup the alignment dict for the current locus

        if not all(spec in list(aln.keys()) for spec in tip_pair):
            continue;
        # Skip any alignments that don't have both tips        

        pair_sub_aln, pair_sub_aln_len = getSubAln(tip_pair, aln, aln_skip_chars);
        # Subset the alignment to contain only the current tip pair

        pair_diffs = getHamming(pair_sub_aln[tip_pair[0]], pair_sub_aln[tip_pair[1]]);
        # Calculate pairwise diffs between the two tips

        if outgroup in aln:
            cur_trio_list = tip_pair + [ outgroup ];
            # Subset the alignment to contain only the current tip pair and the outgroups

            og_sub_aln, og_sub_aln_len = getSubAln(cur_trio_list, aln, aln_skip_chars);
            # Subset the alignment to contain only the current tip pair and the outgroups

            pair_trio_diffs = getHamming(og_sub_aln[tip_pair[0]], og_sub_aln[tip_pair[1]]);
            # Calculate pairwise diffs between the two tips in the trio given the outgroup must also be present

            p1_out_diffs = getHamming(og_sub_aln[tip_pair[0]], og_sub_aln[outgroup]);
            p2_out_diffs = getHamming(og_sub_aln[tip_pair[1]], og_sub_aln[outgroup]);
            # Calculate pairwise diffs between each tip and the outgroup
        else:
            pair_trio_diffs = "NA";
            p1_out_diffs = "NA";
            p2_out_diffs = "NA";
            # If the outgroup is not present in the alignment, set the pairwise diffs values to NA
        
        outline = [ locus, str(pair_sub_aln_len), str(pair_diffs), str(og_sub_aln_len), str(pair_trio_diffs), str(p1_out_diffs), str(p2_out_diffs) ];
        outlines.append(outline);
        # Append the current outline to the list of outlines for all genes        

    return [outlines, tip_pair];

#############################################################################

num_procs = 12;
proc_pool = mp.Pool(processes=num_procs);

#treefile = "../analysis/02-mus-t-windows-new-tree/04-iqtree-no-pahari/chr17/10kb/concat/chr17-concat.cf.tree.rooted";
alndir = "../analysis/02-mus-t-windows-new-tree/03-fasta-no-anc-filter-no-pahari/chr17-filter/";
outdir = "../data/rnd-new-tree/";
sisdir = os.path.join(outdir, "sister");
nonsisdir = os.path.join(outdir, "non-sister");

for d in [outdir, sisdir, nonsisdir]:
    if not os.path.isdir(d):
        os.makedirs(d);

outgroup = "caroli";
aln_skip_chars = ["-", "N", "n"];

####################

pair_outfile = os.path.join(outdir, "pair-counts.csv");
gene_outfile = os.path.join(outdir, "gene-counts.csv");

####################

print(core.getDateTime() + " | reading tree...");
#tree_str = open(treefile, "r").read();
tree_str = "(pahari:1.0,(caroli:1.0,(spretus:1.0,(spicilegus:1.0,(mm10:1.0,mus-t:1.0)Anc4:1.0)Anc3:1.0)Anc2:1.0)Anc1:1.0)Anc0;";
tree = tp.Tree(tree_str);
#tree.showAttrib("length", "anc", "sis", "type");
# Read the tree

tree_in_tips = [ tip for tip in tree.tips if tip != outgroup ];
# Get the tips that aren't outgroups in the tree

tip_pairs = list(combinations(tree_in_tips, 2));
tip_pairs = [ sorted(tip_pair) for tip_pair in tip_pairs ];
# Get all combinations of tip pairs in the tree

sister_pairs = tree.getSisPairs();
# Get all sister pairs in the tree

print(core.getDateTime() + " | " + str(len(sister_pairs)));
print(core.getDateTime() + " | " + str(len(tip_pairs)));
# Get all tip pairs that are direct sisters in the tree and remove them from the 
# tip_pairs list

####################

print(core.getDateTime() + " | reading alns...");
alns = {};
alns_read = 0;

aln_files = os.listdir(alndir);
for aln_file in aln_files:
    if not aln_file.endswith(".fa"):
        continue;
    # Skip files without the .fa exentsion

    aln = os.path.join(alndir, aln_file);
    # Get the full path to the aln file

    cur_aln = seq.fastaReadSeqs(aln);
    # Read the alignment and save to the alns dict

    # if not any([ og in cur_aln for og in outgroups ]):
    #     continue;
    # Skip the alignment if it doesn't contain any outgroups

    alns[aln_file] = cur_aln;
    # Save the alignment 

    alns_read += 1;
    # if alns_read > 10:
    #     break;
    # For testing: stop after a certain number of alignments have been read

print(core.getDateTime() + " | " + str(len(alns)));
# Print total alignments read

####################

with proc_pool as pool:
    print(core.getDateTime() + " | starting counts...");

    cols = [ "gene", "pair.aln.len", "p1.p2.diffs", "trio.aln.len", "p1.p2.trio.diffs", "p1.out.diffs", "p2.out.diffs" ];
    tip_count = 0;

    # for tip_pair in tip_pairs:
    #     print(core.getDateTime(), " | ", tip_count, " ", tip_pair);
    #     result, tmp = pairReps((tip_pair, tree, alns, reps, aln_skip_chars));

    #     if tip_pair in sister_pairs:
    #         pair_out = os.path.join(sisdir, "-".join(tip_pair) + ".csv");
    #     else:
    #         pair_out = os.path.join(nonsisdir, "-".join(tip_pair) + ".csv");

    #     with open(pair_out, "w") as outfile:
    #         outfile.write(",".join(cols));

    #         for outline in result:
    #             print(outline);
    #             outfile.write(",".join(outline) + "\n");

    #     tip_count += 1;
    # Serial version

    for result in pool.imap_unordered(pairReps, ((tip_pair, outgroup, alns, aln_skip_chars) for tip_pair in tip_pairs)):
        outlines_result, cur_pair = result
        print(core.getDateTime(), " | ", tip_count, " ", cur_pair);

        if cur_pair in sister_pairs:
            pair_out = os.path.join(sisdir, "-".join(cur_pair) + ".csv");
        else:
            pair_out = os.path.join(nonsisdir, "-".join(cur_pair) + ".csv");

        with open(pair_out, "w") as outfile:
            for outline_result in outlines_result:
                outfile.write(",".join(outline_result) + "\n");
        
        tip_count += 1;
    # Parallel version
    ## End tip pair loop

####################


##########
## LAST OUTPUT
##########
# ┌─[ gthomas@holybioinf:scripts > conda:base > git:mus-t-haplotype:main ]
# └─[ 11:39:42 ] $: time -p python 20_tip_intro.py
# 06.04.2024 | 11:40:35 | reading tree...
# 06.04.2024 | 11:40:35 | 2
# 06.04.2024 | 11:40:35 | 6
# 06.04.2024 | 11:40:35 | reading alns...
# 06.04.2024 | 11:40:36 | 9499
# 06.04.2024 | 11:40:36 | starting counts...
# 06.04.2024 | 12:31:40  |  0   ['mm10', 'mus-t']
# 06.04.2024 | 12:33:39  |  1   ['mus-t', 'spretus']
# 06.04.2024 | 12:33:58  |  2   ['mm10', 'spretus']
# 06.04.2024 | 12:34:33  |  3   ['mus-t', 'spicilegus']
# 06.04.2024 | 12:34:47  |  4   ['mm10', 'spicilegus']
# 06.04.2024 | 12:35:47  |  5   ['spicilegus', 'spretus']
# real 3311.89
# user 19207.35
# sys 2.83

##########
## NEW TREE OUTPUT
##########
# time -p python 20_tip_intro.py
# 06.11.2024 | 15:04:35 | reading tree...
# 06.11.2024 | 15:04:35 | 1
# 06.11.2024 | 15:04:35 | 10
# 06.11.2024 | 15:04:35 | reading alns...
# 06.11.2024 | 15:06:31 | 9499
# 06.11.2024 | 15:06:31 | starting counts...
# 06.11.2024 | 15:06:32  |  0   ['pahari', 'spretus']
# 06.11.2024 | 15:06:32  |  1   ['pahari', 'spicilegus']
# 06.11.2024 | 15:06:33  |  2   ['mm10', 'pahari']
# 06.11.2024 | 15:06:33  |  3   ['mus-t', 'pahari']
# 06.11.2024 | 16:25:34  |  4   ['mm10', 'mus-t']
# 06.11.2024 | 16:26:58  |  5   ['mm10', 'spicilegus']
# 06.11.2024 | 16:27:15  |  6   ['mm10', 'spretus']
# 06.11.2024 | 16:27:23  |  7   ['mus-t', 'spicilegus']
# 06.11.2024 | 16:27:39  |  8   ['mus-t', 'spretus']
# 06.11.2024 | 16:28:15  |  9   ['spicilegus', 'spretus']
# real 5020.13
# user 28978.60
# sys 6.01


##########
## STASH
##########

######################################

    # if all([ og in cur_aln for og in outgroups ]):
    # ## If both outgroups are in the alignment, include only sites where they share identical alleles

    #     aln_len = len(cur_aln[list(cur_aln.keys())[0]]);
    #     # Get the alignment length from the first sequence

    #     cur_aln = { spec : list(cur_aln[spec]) for spec in cur_aln }; 
    #     # Convert each alignment to a list

    #     for j in range(aln_len):
    #     ## Loop over every site in the current alignment

    #         if len(set([ cur_aln[og][j] for og in outgroups ])) == 1:
    #         # If the set containing all outgroup alleles has length 1, then they are identical

    #             for spec in cur_aln:
    #                 cur_aln[spec][j] == "";
    #             # Loop over every species in the alignment and replace alleles with empty strings,
    #             # which will be removed when the lists are joined
    #     ## End site loop

    #     for spec in cur_aln:
    #         cur_aln[spec] = "".join(cur_aln[spec]);
    #     # Join every species sequence list back to a string, removing empty strings in the process
    # ## End outgroup block



    ######################################
    ##########

    # pair_lca = tree.LCA(tip_pair);
    # # The LCA of the pair of tips

    # lca_clade = tree.getClade(pair_lca);
    # # The clade descending from the LCA

    # ##########

    # p2_possible = [ tip for tip in lca_clade if tip not in tip_pair ];
    # # Possible p2 tips include those descending from the LCA but NOT one of the tips in the pair

    # p2_chosen = random.choices(p2_possible, k=reps);
    # # Randomly choose one of the possible p2 tips
    # ## P2
    
    # ##########

    # outgroups = ['Lophuromys_woosnami_LSUMZ37793', 'Lophiomys_imhausi_UM5152'];

    # #out_possible = [ tip for tip in tree.tips if tip not in lca_clade ];
    # # Possible outgroup tips include any not descending from the LCA

    # #out_chosen = random.choices(out_possible, k=reps);
    # # Randomly choose one of the possible outgroup tips
    # ## Outgroup

    # ##########

    # for i in range(len(p2_chosen)):
    # ## Loop over every selected p2 and outgroup pair

    #     p1 = tip_pair[0];
    #     p2 = p2_chosen[i];
    #     p3 = tip_pair[1];
    #     #out = out_chosen[i];
    #     # Explicitly define p1, p2, p3, and out for clarity

    #     ##########

    #     lca_d1 = tree.desc[pair_lca][0];
    #     lca_d1_clade = tree.getClade(lca_d1);
    #     lca_d2 = tree.desc[pair_lca][1];
    #     lca_d2_clade = tree.getClade(lca_d2)
    #     # Get the two descendant nodes from the p1-p3 LCA

    #     p1p2_sis = False;
    #     # Set the p1-p2 sister flag to False

    #     if p1 in lca_d1_clade and p2 in lca_d1_clade:
    #         p1p2_sis = True;
    #     elif p1 in lca_d2_clade and p2 in lca_d2_clade:
    #         p1p2_sis = True;
    #     # Check if both p1 and p2 descend from the same branch from the LCA of p1 and p3, if so
    #     # set p1-p2 sister flag to True

    #     ## To determine the site patterns to use for calculating D below, we need to know whether
    #     ## p1 is sister to p2 or p3 is sister to p2
    #     ## (((p1,p2),p3),out): d = (baba - abba) / (baba + abba)
    #     ## (((p3,p2),p1),out): d = (baba - aabb) / (baba + aabb)

    #     ##########

    #     quartet = [p1, p2, p3];
    #     # Set the tips as the current quartet

    #     #print("\t", quartet);
    #     #print("\t", p1p2_sis);

    #     ##########

    #     outline = [ p1, p2, p3, str(p1p2_sis) ];
    #     # Initialize the output line for the current quartet

    #     ##########

    #     invariant_sites, variant_sites, decisive_sites = 0, 0, 0;
    #     aabb, baba, abba = 0, 0, 0;
    #     # Initialize site counts

    #     ##########

    #     for locus in alns:
    #     ## Loop over every alignment

    #         quartet += [ og for og in outgroups if og in alns[locus] ];

    #         assert len(quartet) >= 4, "no quartet";

    #         aln = alns[locus];
    #         # Lookup the alignment dict for the current locus

    #         aln_len = len(aln[list(aln.keys())[0]]);
    #         # Get the alignment length from the first sequence

    #         if not all(spec in list(aln.keys()) for spec in quartet):
    #             continue;
    #         # Skip any alignments that don't have all 4 species from the current quartet

    #         sub_aln = { spec : aln[spec] for spec in quartet };
    #         # Get the sub-alignment for the current quartet

    #         #sites = zip(*sub_aln);
    #         # Transpose the alignment to loop over as sites

    #         ##########

    #         for j in range(aln_len):
    #         #for site in sites:

    #             site = [ sub_aln[spec][j] for spec in quartet ];

    #             if site[-1] != site[-2]:
    #                 continue;
    #             # If the outgroups don't share the same allele, skip

    #             if any(char in site for char in aln_skip_chars):
    #                 continue;
    #             # We only care about sites with full information for the quartet, so skip others here

    #             site = site[:-1];
    #             # Remove one of the outgroups to reduce back to a quartet

    #             uniq_site = set(site);
    #             # The alleles at the current site

    #             if len(uniq_site) == 1:
    #                 invariant_sites += 1;
    #             # If the number of unique alleles in the site is 1, it is invariant

    #             elif len(uniq_site) >= 2:
    #             # If the number of unique alleles in the site is greater than 1, it could be informative

    #                 variant_sites += 1;
    #                 # It is definitely variant

    #                 if len(uniq_site) == 2:
    #                     if all(site.count(allele) == 2 for allele in uniq_site):
    #                         decisive_sites += 1;
    #                         # If the number of uniqe alleles at the site is 2 and they are both present in 2 species, the site is decisive

    #                         if site[0] == site[1] and site[2] == site[3]:
    #                             aabb += 1;
    #                         elif site[0] == site[2] and site[1] == site[3]:
    #                             baba += 1;
    #                         elif site[0] == site[3] and site[1] == site[2]:
    #                             abba += 1;
    #                         # Count the various site patterns for decisive sites
    #         ## End site loop

    #         ##########

    #         if p1p2_sis:
    #             d_num = baba - abba;
    #             d_den = baba + abba;
    #         # Calculate the numerator and denominator for D when p1 and p2 are sister

    #         else:
    #             d_num = baba - aabb;
    #             d_den = baba + aabb;
    #         # Calculate the numerator and denominator for D when p3 and p2 are sister

    #         if d_den:
    #             d = d_num / d_den;
    #         else:
    #             d = "NA";
    #         # Calculate d, accounting for loci with 0 counts by setting them as NA

    #         ##########
    #     ## End aln loop

    #     outline += [ str(invariant_sites), str(variant_sites), str(decisive_sites), str(aabb), str(baba), str(abba), str(d) ];
    #     outlines.append(outline);
    #     # Compile the rest of the output and append the outline for the current quartet-locus pair to the list of outlines

    #     ##########        
    # ## End quartet loop
    ######################################