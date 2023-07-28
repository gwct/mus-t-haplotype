#!/usr/bin/python
############################################################
# Counts topologies from window-based trees
############################################################

import sys, os, re, argparse, lib.core as CORE, lib.treeparse as TP
from collections import defaultdict

############################################################
# Functions

def treeDictRmBlength(tdict):
    for node in tdict:
        tdict[node][0] = 'NA';
    return tdict;

############################################################
# Globals
#chrome_sizes = {"1" : 195471971, "2" : 182113224, "3" : 160039680, "4" : 156508116, "5" : 151834684, "6" : 149736546, "7" : 145441459,
#                "8" : 129401213, "9" : 124595110, "10" : 130694993, "11" : 122082543, "12" : 120129022, "13" : 120421639, "14" : 124902244,
#                "15" : 104043685, "16" : 98207768, "17" : 94987271, "18" : 90702639, "19" : 61431566, "X" : 171031299 }; #, "Y" : 91744698

#chrome_sizes = { "19" : 61431566 }; #, "Y" : 91744698

chrome_sizes = { "17" : 94987271 }
# Mouse chromosome sizes from /mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.chrome.sizes"

############################################################

parser = argparse.ArgumentParser(description="Rodent window topology counter");
parser.add_argument("-i", dest="window_file", help="A file that has information about the alignment windows in a certain format...", default=False);
parser.add_argument("-w", dest="window_size", help="The size of the sliding window in kb.", default=False);
parser.add_argument("-t", dest="base_tree_dir", help="A directory with subdirectories for chromosomes and window sizes that contains the output from iqtree and astral.", default=False);
parser.add_argument("-o", dest="outfile", help="A csv file to write the topology counts to", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
args = parser.parse_args();
# Options

if not args.window_size:
    sys.exit(" * Error 1: Window size in kb (-w) must be defined.");
if int(args.window_size) < 1:
    sys.exit(" * Error 2: Window size in kb (-w) must be a positive integer.");
else:
    wsize = float(args.window_size) * 1000;
    wsize_str = str(wsize);
# Parse the window size

if not args.window_file or not os.path.isfile(args.window_file):
    sys.exit(" * Error 3: A valid windows .tsv file must be provided with -i.");
else:
    windowfile = args.window_file;
#windowfile = "get-windows-default.tsv";
# Parse the input window file

if not args.base_tree_dir or not os.path.isdir(args.base_tree_dir):
    sys.exit(" * Error 3: A valid directory must be provided with -t.");
else:
    base_treedir = args.base_tree_dir;
#base_treedir = "/n/holylfs05/LABS/informatics/Users/gthomas/mus-t-haplotype/cactus-all-mask/windows/tree/";
# Check the tree directory

in_prefix = args.window_size + "kb";
out_prefix = args.window_size + "kb";

if not args.outfile:
    outfilename = out_prefix + "-topo-counts.csv";
else:
    outfilename = args.outfile;
# Get the output file

if os.path.isfile(outfilename) and not args.overwrite:
    sys.exit( " * Error: Output file (" + outfilename + ") already exists! Explicity specify --overwrite to overwrite it.");
# Make sure --overwrite is set if the output file already exists
# File names

###################

headers = ["window", "chr", "chr.len", "start", "end", 
            "filter", 
            "gcf", "scf", 
            "topo.num.overall", "topo.count.overall", "topo.rank.overall", "topo.color", 
            "topo.num.chrome", "topo.count.chrome", "topo.rank.chrome", 
            "topo", "concat.chrome.topo.match", "astral.chrome.topo.match", "unparsed.tree"];
# Column headers for the output file.

###################

with open(outfilename, "w") as outfile:
    CORE.runTime("# Window-based topology counter", outfile);
    CORE.PWS("# Window size:          " + wsize_str, outfile);
    CORE.PWS("# Window file:          " + windowfile, outfile);
    CORE.PWS("# Tree dir:             " + base_treedir, outfile);
    CORE.PWS("# Output file:          " + outfilename, outfile);
    CORE.PWS("# ----------------");
    # Input info for log file

    CORE.PWS("# " + CORE.getDateTime() + " Reading windows...", outfile);
    windows, marker_win_trees, first = {}, {}, True;
    for line in open(windowfile):
        if line[0] == "#":
            continue;
        if first:
            first = False;
            continue;
        line = line.strip().split("\t");
        #print(line);

        window, chrome, start, end = line[0].replace(".noanc", ""), line[1], line[2], line[3].replace(".noanc", "");

        filter_str = "FILTER";
        if int(line[32]) == 6 and int(line[34]) >= 4 and line[42] == "0":
            filter_str = "PASS";

        chr_len = str(chrome_sizes[chrome.replace("chr", "")]);

        windows[window] = { 'chr' : chrome, 'chr.len' : chr_len, 'start' : str(start), 'end' : str(end), 
                            'filter' : filter_str,
                            'gcf' : "NA", 'scf' : "NA",
                            'topo.num.overall' : "NA", 'topo.count.overall' : "NA", 'topo.rank.overall' : "NA", 'topo.color' : "NA",
                            'topo.num.chrome' : "NA", 'topo.count.chrome' : "NA", 'topo.rank.chrome' : "NA", 
                            'topo' : "NA", 'concat.chrome.topo.match' : "NA", 'astral.chrome.topo.match' : "NA", 'unparsed.tree' : "NA", 'clade.set' : "NA" };
        # Initialize info for this tree window.
    CORE.PWS("# Windows read:         " + str(len(windows)), outfile);
    CORE.PWS("# ----------------");
    # Read the tree windows.

    topo_colors, num_colors, colors = {}, 0, ["#db6d00", "#004949", "#006ddb", "#920000", "#490092", "#6cb6ff", "#24ff24", "#fdb4da", "#ffff6d", "#009292"];
    # Variables for assigning colors to each topology.

    uniq_topos_all, topo_counts_all = [], defaultdict(int);
    # Variables for overall counting.

    ##########

    for chromosome in chrome_sizes:
    # Count trees for each chromosome.

        chrstr = "chr" + chromosome
        chr_end = chrome_sizes[chromosome];
        chr_end_str = str(chr_end);
        # Chromosome variables.

        indir = os.path.join(base_treedir, chrstr, in_prefix);
        # Tree directory for this chromosomes.
        
        locidir = os.path.join(indir, "loci");
        assert os.path.isdir(locidir), "\nCANNOT FIND TREE DIRECTORY: " + locidir;

        concattreefile = os.path.join(indir, "concat", chrstr + "-concat.treefile.rooted");
        assert os.path.isfile(concattreefile), "\nCANNOT FIND CONCATENATED TREE: " + concattreefile;

        concordfile = os.path.join(indir, "astral", chrstr + "-astral.cf.stat");
        assert os.path.isfile(concordfile), "\nCANNOT FIND CONCORDANCE FILE: " + concordfile;

        astraltreefile = os.path.join(indir, "astral", chrstr + "-astral.treefile.rooted");
        assert os.path.isfile(concattreefile), "\nCANNOT FIND ASTRAL TREE: " + astraltreefile;
        # Tree files for this chromosome.

        CORE.PWS("# Chromosome:           " + chromosome, outfile);
        CORE.PWS("# Chromosome length:    " + chr_end_str, outfile);
        CORE.PWS("# Locus tree directory: " + locidir, outfile);
        CORE.PWS("# Concat tree file:     " + concattreefile, outfile);
        CORE.PWS("# ASTRAL tree file:     " + astraltreefile, outfile);
        CORE.PWS("# Concordance file:     " + concordfile, outfile);
        CORE.PWS("# ----------------", outfile);
        # Log info for this chromosome.

        CORE.PWS("# " + CORE.getDateTime() + " Reading concatenated chromosome tree...", outfile);
        concattree = open(concattreefile, "r").read().strip();
        CORE.PWS("# Concat chromosome tree:      \"" + concattree + "\"", outfile);
        tinfo, concat_tree, root = TP.treeParse(concattree);
        concat_tree_topo = set([frozenset(TP.getClade(node, tinfo)) for node in tinfo if tinfo[node][2] != 'tip']);
        CORE.PWS("# Concat chromosome topology:  \"" + concat_tree + "\"", outfile);
        # Read the concatenated tree for this chromosome.

        CORE.PWS("# " + CORE.getDateTime() + " Reading ASTRAL chromosome tree...", outfile);
        astraltree = open(astraltreefile, "r").read().strip();
        CORE.PWS("# ASTRAL chromosome tree:      \"" + astraltree + "\"", outfile);
        tinfo, astral_tree, root = TP.treeParse(astraltree);
        astral_tree_topo = set([frozenset(TP.getClade(node, tinfo)) for node in tinfo if tinfo[node][2] != 'tip']);
        CORE.PWS("# ASTRAL chromosome topology:  \"" + astral_tree + "\"", outfile);
        # Read the ASTRAL tree for this chromosome.

        astral_concat_match = "FALSE";
        if concat_tree_topo == astral_tree_topo:
            astral_concat_match = "TRUE";
        CORE.PWS("# ASTRAL-concat match:   " + astral_concat_match, outfile);
        CORE.PWS("# ----------------");
        # Check if the conactenated and ASTRAL topologies are the same for this chromosome.

        CORE.PWS("# " + CORE.getDateTime() + " Calculating average CFs...", outfile);
        gcfs, scfs, first = [], [], True;
        for line in open(concordfile):
            if line[0] == "#":
                continue;
            if first:
                first = False;
                continue;
            line = line.strip().split("\t");
            if line[1] != "NA":
                gcfs.append(float(line[1]));
            if line[10] != "NA":
                scfs.append(float(line[10]));

        avg_gcf = str(CORE.mean(gcfs))
        avg_scf = str(CORE.mean(scfs));
        CORE.PWS("# Average gCF:          " + avg_gcf, outfile);
        CORE.PWS("# Average sCF:          " + avg_scf, outfile);
        CORE.PWS("# ----------------");
        # Using the IQ-tree stat file to get average CFs for this chromosome.

        CORE.PWS("# " + CORE.getDateTime() + " Reading windows and counting topologies...", outfile);
        uniq_topos, topo_counts = [], defaultdict(int);
        # Lists of unique clade sets for each window and the unique lists of clade sets to count topologies.

        num_trees, chrome_to_all = 0, {};
        # The number of trees read and the conversion between chromosome topology number and overall topology number for the coloring.

        filtered_windows, no_tree_file = 0,0;

        for window in windows:
            if windows[window]['chr'] != chrstr:
                continue;
            # Get the tree for each window on this chromosome.

            windows[window]['gcf'] = avg_gcf;
            windows[window]['scf'] = avg_scf;
            # Add concordance factors.

            if windows[window]['filter'] == "PASS":
            # If the window was previously filtered then skip it.

                #window_str = window.replace(":", "-");
                window_name = "chr" + chromosome + ":" + windows[window]['start'] + "-" + windows[window]['end'];
                window_dir = os.path.join(locidir, window_name);
                cur_window_files = os.listdir(window_dir);
                treefile_list = [ f for f in cur_window_files if f.endswith(".filter.treefile.rooted") ];
                assert len(treefile_list) == 1, "\nMORE THAN ONE TREE FILE FOUND: " + os.path.join(locidir, window_name);
                treefile = os.path.join(window_dir, treefile_list[0]);
                #treefile = os.path.join(locidir, window_name, window_name + ".noanc.filter.treefile.rooted");
                # Get window name and tree file.
                
                if not os.path.isfile(treefile):
                    print("# CANNOT FIND TREE FILE, SETTING STATUS TO FILTER: " + treefile);
                    no_tree_file += 1;
                    windows[window]['missing-filter'] = "FILTER";
                    continue;
                # If the file isn't found because one of the alignment or tree steps failed, just set it to filter here.

                tree = open(treefile, "r").read().strip();
                # Read the tree from the file.

                tinfo, tree_parsed, root = TP.treeParse(tree);
                #tree_parsed = re.sub("<[\d]+>", "", tree_parsed);
                # Parse the tree.

                clade_topo = set([frozenset(TP.getClade(node, tinfo)) for node in tinfo if tinfo[node][2] != 'tip']);
                # Get the topology as the set of all clades present in the tree.

                if clade_topo not in uniq_topos_all:
                    uniq_topos_all.append(clade_topo);
                # If this topology hasn't been seen before at all, add it to the list of unique topologies.

                if clade_topo not in uniq_topos:
                    uniq_topos.append(clade_topo);
                # If this topology hasn't been seen before on this chromosome, add it to the list of unique topologies.

                concat_match = "FALSE";
                if clade_topo == concat_tree_topo:
                    concat_match = "TRUE";
                astral_match = "FALSE";
                if clade_topo == astral_tree_topo:
                    astral_match = "TRUE";
                # Check if the current topology matches the chromosome level topologies.

                topo_num_all = uniq_topos_all.index(clade_topo)+1;
                topo_counts_all[topo_num_all] += 1;
                topo_num = uniq_topos.index(clade_topo)+1;
                topo_counts[topo_num] += 1;
                chrome_to_all[topo_num] = topo_num_all;
                # The topology ID number is then the index in the unique topology list.

                windows[window]['topo.num.overall'] = topo_num_all;
                windows[window]['topo.num.chrome'] = topo_num;
                windows[window]['topo'] = tree_parsed + ";";
                windows[window]['concat.chrome.topo.match'] = concat_match;
                windows[window]['astral.chrome.topo.match'] = astral_match;
                windows[window]['unparsed.tree'] = tree;
                windows[window]['clade.set'] = clade_topo;
                # Save output info for this window.

                num_trees += 1;
                # Increment the number of trees read.
            # If this window has a tree.
            else:
                filtered_windows += 1;
        CORE.PWS("# Total trees read:   " + str(num_trees), outfile);
        CORE.PWS("# No tree file:       " + str(no_tree_file), outfile);
        CORE.PWS("# Filtered windows:   " + str(filtered_windows), outfile);
        CORE.PWS("# Unique topologies:  " + str(len(uniq_topos)), outfile);
        CORE.PWS("# ----------------", outfile);

        CORE.PWS("# " + CORE.getDateTime() + " Counting and ranking topologies for this chromosome...", outfile);
        topo_counts_sorted = sorted(topo_counts.items(), key=lambda x: x[1], reverse=True);
        # Sort the counts
        topo_ranks, j = {}, 1;
        for i in topo_counts_sorted:
            topo_ranks[i[0]] = j;
            j += 1;
        # Rank the sorted counts

        for topo in topo_ranks:
            if topo_ranks[topo] <= 3:
                topo_num_all = chrome_to_all[topo];
                if topo_num_all not in topo_colors:
                    topo_colors[topo_num_all] = colors[num_colors];
                    num_colors += 1;
        # Assign a color to the top 3 topologies if they haven't been assigned already.

        for window in windows:
            if windows[window]['chr'] != chrstr:
                continue;
            if windows[window]['filter'] == "PASS":
                chrome_topo_count = topo_counts[windows[window]['topo.num.chrome']];
                chrome_topo_rank = topo_ranks[windows[window]['topo.num.chrome']];
                #print(window, chrome_topo_count, chrome_topo_rank);
                windows[window]['topo.count.chrome'] = str(chrome_topo_count);
                windows[window]['topo.rank.chrome'] = str(chrome_topo_rank);
            # If this window had a tree made, count the number of times that topology occurred.
        CORE.PWS("# ----------------", outfile);
    # END CHROMSOME LOOP
    ####################

    CORE.PWS("# " + CORE.getDateTime() + " Counting and ranking overall topologies...", outfile);
    topo_counts_sorted = sorted(topo_counts_all.items(), key=lambda x: x[1], reverse=True);
    # Sort the counts
    topo_ranks, j = {}, 1;
    for i in topo_counts_sorted:
        topo_ranks[i[0]] = j;
        j += 1;
    # Rank the sorted counts

    for window in windows:
        if windows[window]['filter'] == "PASS":
            #print(window, windows[window]);
            #print(topo_counts_all);
            topo_count = topo_counts_all[windows[window]['topo.num.overall']];
            topo_rank = topo_ranks[windows[window]['topo.num.overall']];
            windows[window]['topo.count.overall'] = str(topo_count);
            windows[window]['topo.rank.overall'] = str(topo_rank);
    # If this window had a tree made, count the number of times that topology occurred.
    CORE.PWS("# ----------------", outfile);         

    print(topo_colors)
    for topo in topo_ranks:
        if topo not in topo_colors:
            topo_colors[topo] = "#999999";
    # Fill in the remaining colors with grey

    CORE.PWS("# " + CORE.getDateTime() + " Writing output: " + outfilename, outfile);
    outfile.write(",".join(headers) + "\n");

    for window in windows:
        w = windows[window];

        w['topo.color'] = "NA";
        if w['topo.num.overall'] != "NA":
            w['topo.color'] = topo_colors[w['topo.num.overall']];

        outline = [window];
        for header in headers[1:]:
            outline.append('"' + str(w[header]) + '"');
        outfile.write(",".join(outline) + "\n");
        # Write the current output (minus the clade topology if a tree is present for window);
    CORE.PWS("# ----------------");
    print("# " + CORE.getDateTime() + " Done!");
