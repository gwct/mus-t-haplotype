#!/usr/bin/python
############################################################
# For turtles genomes, 01.2023
# Build a table of windows and assess filters (e.g. missing
# data, repeats) from a whole genome alignment of multiple
# species
############################################################

import sys, os, argparse, core, seqparse as SEQ
import multiprocessing as mp

############################################################

def isPosFloat(x, x_min=0.0, x_max=1.0):
    try:
        x = float(x);
    except:
        return False;

    if x < x_min or x > x_max:
        return False;
    else:
        return True;

#########################

def processScaff(scaff_item):

    scaff_id, base_indir, cur_scaff_num, num_scaffs, base_outdir = scaff_item;
    scaff_outlines = [];

    scaff_indir = os.path.join(base_indir, scaff_id);

    if base_outdir:
        scaff_outdir = os.path.join(base_outdir, scaff_id +"-filter");
        if not os.path.isdir(scaff_outdir):
            os.makedirs(scaff_outdir);

    #scaff_dir_base = os.path.basename(scaff_dir);
    #print("# " + core.getDateTime() + " | " + str(cur_scaff_num) + " / " + str(num_scaffs) + " - " + scaff_dir_base);
    # Status update

    scaff_aln_files = os.listdir(scaff_indir);
    # Get the list of all alignments for the current scaffold

    num_alns = len(scaff_aln_files);
    cur_aln_num = 1;
    # Counting for status updates

    for scaff_aln_file in scaff_aln_files:
        print("# " + core.getDateTime() + " | " + scaff_id + " > " + str(cur_aln_num) + " / " + str(num_alns) + " - " + scaff_aln_file);
        # Status update

        file_base = os.path.splitext(scaff_aln_file)[0];
        file_list = file_base.split(":");
        scaff = file_list[0]
        start, end = file_list[1].split("-");
        # Parse the window info out of the alignment file name

        cur_outline = [file_base, scaff, start, end];
        # Add the window info to the output list

        cur_aln_raw = SEQ.fastaReadSeqs(os.path.join(scaff_indir, scaff_aln_file));
        cur_aln = { header : cur_aln_raw[header] for header in cur_aln_raw if "Anc" not in header };
        num_seqs = len(cur_aln);
        aln_len = len(cur_aln[list(cur_aln.keys())[0]]);
        # Read the alignment and remove ancestral sequences

        seq_stats = countSeqStats(cur_aln);
        #avg_no_missing_len, num_missing, num_repeat, num_gappy = countNonMissingLength(cur_aln, args.missing_thresh, args.gap_thresh);
        # Count aln stats by sequence

        site_stats = siteCount(cur_aln, aln_len);
        #invariant_sites, informative_sites, sites_w_missing, perc_sites_w_missing, high_missing_sites, perc_sites_high_missing, all_missing_sites, sites_w_gap, perc_sites_w_gap, high_gap_sites, perc_high_gap_sites, all_gap_sites, mg_sites, percent_mg_sites = siteCount(cur_aln, aln_len, args.missing_thresh, args.gap_thresh);
        # Count aln stats by site

        filtered_aln, gappy_cols = windowSiteFilter(cur_aln, aln_len);

        num_seqs_f = len(filtered_aln);
        filtered_aln_len = len(cur_aln[list(filtered_aln.keys())[0]]);
        
        cur_outline = cur_outline + seq_stats + site_stats;

        ###################

        seq_stats_f = countSeqStats(filtered_aln);
        # Count aln stats by sequence
        #invariant_sites_f, informative_sites_f, sites_w_missing_f, perc_sites_w_missing_f, high_missing_sites_f, perc_sites_high_missing_f, all_missing_sites_f, sites_w_gap_f, perc_sites_w_gap_f, high_gap_sites_f, perc_high_gap_sites_f, all_gap_sites_f, mg_sites_f, percent_mg_sites_f = siteCount(filtered_aln, filtered_aln_len, args.missing_thresh, args.gap_thresh);
        # Count aln stats by site

        site_stats_f = siteCount(cur_aln, filtered_aln_len);
        cur_outline = cur_outline + seq_stats_f + site_stats_f + ["FALSE"];
        
        ###################

        if base_outdir:
            cur_aln_outfile = os.path.join(scaff_outdir, scaff_aln_file.replace(".fa", ".filter.fa"));
            with open(cur_aln_outfile, "w") as cur_outfile:
                for header in filtered_aln:
                    cur_outfile.write(">" + header + "\n");
                    cur_outfile.write(filtered_aln[header] + "\n");
            cur_outline[-1] = "TRUE";

        cur_aln_num += 1;
        # Increment aln counter

        cur_outline_str = [ str(col) for col in cur_outline ];
        scaff_outlines.append(cur_outline_str);          
        # Add stats to output list, convert all entries to strings, and write to file
    ## End aln loop

    return scaff_outlines;

#########################

def countSeqStats(seqs):
# This function goes through every sequence in an alignment and calculates
# the average length of each sequence excluding gaps/missing data and counts
# how many seqs above filter thresholds.

    # repeat_filter = 0.5;

    gap_chars = "-";
    missing_chars = "NnXx";
    repeat_chars = "atcgnx";
    # Strings of missing and repeat/softmask characters

    uniq, ident = countUniqIdentSeqs(seqs);

    len_sum, num_seqs = 0, 0;
    above_missing, all_missing = 0, 0;
    above_gap, all_gap = 0, 0;
    above_mg, all_mg = 0, 0;
    
    for seq in seqs:
        num_seqs += 1;
        full_len = len(seqs[seq]);

        num_missing, num_gap = 0, 0;
        for site in seqs[seq]:
            if site in missing_chars:
                num_missing += 1;
            if site in gap_chars:
                num_gap += 1;
        # Count missing and gap chars for current seq

        num_mg = num_missing + num_gap;
        # Count total missing + gap chars

        no_missing_len = full_len - num_missing;
        len_sum += no_missing_len;
        # Get len of sequence without missing chars

        if (num_missing / full_len) > args.missing_seq_thresh:
            above_missing += 1;
        if num_missing == full_len:
            all_missing += 1;
        # Check missing

        # no_repeat_len = len( [ site for site in seqs[seq] if site not in repeat_chars ] );

        # if 1 - (no_repeat_len / full_len) > repeat_filter:
        #     num_repeat_seqs += 1;
        # If the number of repeat chars in the sequence (calculated as 1 - the fraction of non-repeat length to full length)
        # is above some threshold, mark it as a repeat seq

        if (num_gap / full_len) > args.gap_seq_thresh:
            above_gap += 1;
        if num_gap == full_len:
            all_gap += 1;
        # Check gap

        if (num_mg / full_len) > args.mg_seq_thresh:
            above_mg += 1;
        if num_mg == full_len:
            all_mg += 1;
        # Check missing + gap
    ## End seq loop

    no_missing_len = len_sum / full_len;
    # Get avg length of seq in aln without missing chars

    return [ num_seqs, full_len, uniq, ident, above_missing, all_missing, no_missing_len, above_gap, all_gap, above_mg, all_mg ]

#########################

def countUniqIdentSeqs(seqs):
# This function goes through every sequence in an alignment and counts how 
# many sequences are unique or identical.

    uniq_seqs, ident_seqs, found = 0, 0, [];
    seq_list_raw = list(seqs.values());
    seq_list = [ seq.replace("-", "") for seq in seq_list_raw ];
    for seq in seq_list:
        if seq_list.count(seq) == 1:
            uniq_seqs += 1;
        if seq_list.count(seq) != 1 and seq not in found:
            ident_seqs += 1;
            found.append(seq);

    return uniq_seqs, ident_seqs;

#########################

def siteCount(seqs, aln_len):
# This function goes through every site in an alignment to check for invariant sites, gappy sites, 
# and stop codons.

    # [ "invariant.sites", "informative.sites" ] 
    # missing_headers = [ "sites.w.missing", "perc.sites.w.missing", "sites.high.missing", "perc.sites.high.missing", "sites.all.missing" ];
    # gap_headers = [ "sites.w.gap", "perc.sites.w.gap", "sites.high.gap", "perc.sites.high.gap", "sites.all.gap" ];
    # mg_headers = [ "sites.missing.gap", "perc.sites.missing.gap", "sites.high.missing.gap", "perc.sites.high.missing.gap", "sites.all.missing.gap" ];
    
    gap_chars = "-";
    missing_chars = "NnXx";

    invariant, informative = 0, 0;
    missing, perc_missing, high_missing, perc_high_missing, all_missing = 0, 0.0, 0, 0.0, 0;
    gap, perc_gap, high_gap, perc_high_gap, all_gap = 0, 0.0, 0, 0.0, 0;
    mg, perc_mg, high_mg, perc_high_mg, all_mg = 0, 0.0, 0, 0.0, 0;
    # Counts

    seq_list = list(seqs.values());
    num_spec = len(seq_list);

    for i in range(aln_len):
    # Loop over every site

        site = [];
        for j in range(num_spec):
            site.append(seq_list[j][i]);
        site_len = len(site);
        # Get the current (ith) site from every sequence (j)

        if site.count(site[0]) == site_len:
            invariant += 1;
        # If all the nts at the site match the first one, it is invariant

        allele_counts = { allele : site.count(allele) for allele in site if allele not in missing_chars+gap_chars };
        # Count the occurrence of each allele in the site

        if len(allele_counts) > 1:
            multi_allele_counts = [ allele for allele in allele_counts if allele_counts[allele] >= 2 ];
            # Count the number of alleles present in at least 2 species

            if len(multi_allele_counts) >= 2:
                informative += 1;
            # If 2 or more alleles are present in 2 or more species, this site is informative

        ##########

        num_missing, num_gap = 0, 0;
        for allele in site:
            if allele in missing_chars:
                num_missing += 1;
            if allele in gap_chars:
                num_gap += 1;
        # Count missing and gap chars for current seq

        #####

        if num_missing > 1:
            missing += 1;
            # Increment by one if there is at least one gap

            if (num_missing / site_len) > args.missing_site_thresh:
                high_missing += 1;
            # Count if the number of gaps at this site is above some threshold

            if num_missing == site_len:
                all_missing += 1;
            # Check if the site is all missing
        # Missing
        #####

        if num_gap > 1:
            gap += 1;
            # Increment by one if there is at least one gap

            if (num_gap / site_len) > args.gap_site_thresh:
                high_gap += 1;
            # Count if the number of gaps at this site is above some threshold

            if num_gap == site_len:
                all_gap += 1;
            # Check if the site is all gaps
        # Gaps
        #####

        if num_missing > 1 or num_gap > 1:
            mg += 1;

            if (num_missing + num_gap) / site_len > args.mg_site_thresh:
                high_mg += 1;

            if (num_missing + num_gap) == site_len:
                all_mg += 1;
        # Gaps OR missing
        ##########
    ## End site loop

    perc_missing = missing / aln_len;
    perc_high_missing = high_missing / aln_len;
    perc_all_missing = all_missing / aln_len;

    perc_gap = gap / aln_len;
    perc_high_gap = high_gap / aln_len;
    perc_all_gap = all_gap / aln_len;

    perc_mg = mg / aln_len;
    perc_high_mg = high_mg / aln_len;
    perc_all_mg = all_mg / aln_len;
    # Calculate the percentage of sites that have at least one missing character or one gap

    ##########

    return [ invariant, informative, missing, perc_missing, high_missing, perc_high_missing, all_missing, gap, perc_gap, high_gap, perc_high_gap, all_gap, mg, perc_mg, high_mg, perc_high_mg, all_mg ];

#########################

def windowSiteFilter(seqs, aln_len):
# This function takes a sliding window along a codon alignment and filters windows
# where 2 or more codons have 2 or more gaps.

    window_size = 5;
    gappy_window = 3;

    num_seqs = len(seqs);

    exclude_sites = [];
    # The indices of the sites to filter

    i = 0;
    while i < (aln_len - window_size):
    # Loop over every site in the alignment

        seq_exclude = 0;
        # The number of sequences at the current window that are too gappy

        for seq in seqs:
        # Go over every sequence in the alignment

            cur_window = seqs[seq][i:i+window_size]
            # For every sequence get the current window of <window_size> nucleotides

            if cur_window.count("-") >= gappy_window:
                seq_exclude += 1;
            # If there are <gappy_window> or more positions that are gaps in the window in the current sequence, add to the
            # list of exclude sequences for this window
        ## End seq loop

        if seq_exclude / num_seqs > args.gap_window_thresh:
            exclude_sites += [pos for pos in range(i, i+window_size)];
        # If the total number of sequences excluded for being gappy at these sites is over some threshold,
        # add all <window_size> sites to the list of excluded sites for this alignment

        i += 1;
    ## End site loop

    exclude_sites_sorted = sorted(list(set(exclude_sites)), reverse=True);
    # Reverse sort the excluded sites to remove with del()

    for seq in seqs:
        seqs[seq] = list(seqs[seq]);
        for site in exclude_sites_sorted:
            del(seqs[seq][site]);
        seqs[seq] = "".join(seqs[seq]);
    #     seqs_filtered = { seq : i for seq,i in seqs.items() }
    #    seqs[seq] = [ seqs[seq][i] for i in range(aln_len) if i not in exclude_sites ];
    # From every sequence, remove the sites from windows determined to be gappy in too many sequences above

    return seqs, len(exclude_sites_sorted);
    #return codon_seqs, len(list(set(exclude_sites)));

############################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Get rodent chromosome window coordinates");
    parser.add_argument("-i", dest="indir", help="The input directory, containing subdirectories for each reference scaffold that contain window alignments.", default=False);
    parser.add_argument("-w", dest="window_size", help="The size of the sliding window in kb.", type=int, default=False);
    parser.add_argument("-gseq", dest="gap_seq_thresh", help="The proportion of a sequence that must be gap characters for it to be considered 'high gap'. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-gsite", dest="gap_site_thresh", help="The proportion of an alignment column that must be gap characters for it to be considered 'high gap'. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-g", dest="gap_window_thresh", help="In sliding windows of 5 columns, if the proportion of sequences that have 3 or more gaps is greater than this threshold, remove all 5 columns of the alignment from every sequence. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-mseq", dest="missing_seq_thresh", help="The proportion of a sequence that must be missing characters for it to be considered 'high missing'. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-msite", dest="missing_site_thresh", help="The proportion of an alignment column that must be missing characters for it to be considered 'high missing'. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-mgseq", dest="mg_seq_thresh", help="The proportion of a sequence that must be missing OR gap characters for it to be considered 'high mg'. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-mgsite", dest="mg_site_thresh", help="The proportion of an alignment column that must be missing OR gap characters for it to be considered 'high mg'. Default: 0.5", type=float, default=0.5);
    parser.add_argument("-o", dest="outfile", help="A file to output the csv values and log info to.", default="get-windows-default.tsv");
    parser.add_argument("-d", dest="outdir", help="If provided, a directory to write the filtered sequences to.", default=False);
    parser.add_argument("-p", dest="procs", help="The number of processes to use. Default: 1", type=int, default=1);
    parser.add_argument("--overwrite", dest="overwrite", help="A file to output the csv values and log info to.", action="store_true", default=False);
    # parser.add_argument("-p", dest="procs", help="The number of processes the script should use. Default: 1.", type=int, default=1);
    args = parser.parse_args();

    #########################

    if not os.path.isdir(args.indir):
        sys.exit(" * Error 1: Input directory (-i) not found.");
    else:
        indir = os.path.abspath(args.indir);

    #####

    if not args.window_size:
        sys.exit(" * Error 2: Window size in kb (-w) must be defined.");
    if args.window_size < 1:
        sys.exit(" * Error 3: Window size in kb (-w) must be a positive integer.");
    else:
        wsize_str = str(args.window_size);

    #####

    for opt in [args.gap_seq_thresh, args.gap_site_thresh, args.gap_window_thresh, args.missing_seq_thresh, args.missing_site_thresh, args.mg_seq_thresh, args.mg_site_thresh ]:
        if not isPosFloat(opt):
            sys.exit(" * Error 4: Thresholds must be a value between 0 and 1.");

    #####

    outfilename = os.path.abspath(args.outfile);

    if os.path.isfile(outfilename) and not args.overwrite:
        sys.exit(" * Error 6: Output file (-o) already exists. Please move the current file or specify to --overwrite it.");

    #####

    if args.outdir:
        if not os.path.isdir(args.outdir):
            os.makedirs(args.outdir);
        elif not args.overwrite:
            sys.exit(" * Error 7: Output directory (-d) already exists. Please move the current directory or specify to --overwrite it.");


    proc_pool = mp.Pool(processes=args.procs);

    #########################

    meta_headers = ["window", "scaffold", "start", "end"]
    seq_headers = [ "total.seqs", "aln.len", "uniq.seqs", "ident.seqs", "seqs.above.missing", "seqs.all.missing", "avg.seq.len.wo.missing", "seqs.above.gappy", "seqs.all.gappy", "seqs.above.missing.gappy", "seqs.all.missing.gappy" ];

    missing_headers = [ "sites.w.missing", "perc.sites.w.missing", "sites.high.missing", "perc.sites.high.missing", "sites.all.missing" ];
    gap_headers = [ "sites.w.gap", "perc.sites.w.gap", "sites.high.gap", "perc.sites.high.gap", "sites.all.gap" ];
    mg_headers = [ "sites.missing.gap", "perc.sites.missing.gap", "sites.high.missing.gap", "perc.sites.high.missing.gap", "sites.all.missing.gap" ];
    
    site_headers = [ "invariant.sites", "informative.sites" ] + missing_headers + gap_headers + mg_headers;

    seq_filter_headers = [ col + ".filter" for col in seq_headers ];
    site_filter_headers = [ col + ".filter" for col in site_headers ];

    headers = meta_headers + seq_headers + site_headers + seq_filter_headers + site_filter_headers + ["written"];

    print("\t".join(headers));
    sys.exit();

    # headers = [ "window", "scaffold", "start", "end", "total.seqs", "aln.len", "seqs.above.missing", "avg.seq.len.wo.missing", "seqs.above.gappy", "gappy.cols", "uniq.seqs", "ident.seqs", "invariant.sites", "informative.sites", "sites.w.missing", 
    # "percent.sites.w.missing", "sites.high.missing", "percent.sites.high.missing", "sites.all.missing", "sites.w.gap", "percent.sites.w.gap", "sites.high.gap", "percent.sites.high.gap", "sites.all.gap", "sites.missing.gap", "percent.sites.missing.gap" ];
    # headers += ["total.seqs.filter", "aln.len.filter", "seqs.above.missing.filter", "avg.seq.len.wo.missing.filter", "seqs.above.gappy.filter", "uniq.seqs.filter", "ident.seqs.filter", "invariant.sites.filter", "informative.sites.filter", "sites.w.missing.filter", "percent.sites.w.missing.filter", "sites.high.missing.filter", "percent.sites.high.missing.filter", "sites.all.missing.filter", "sites.w.gap.filter", "percent.sites.w.gap.filter", "sites.high.gap.filter", "percent.sites.high.gap.filter", "sites.all.gap.filter",  "sites.missing.gap.filter", "percent.sites.missing.gap.filter" ];
    # headers.append("written");

    #########################

    pad = 40;
    header_pad = 40;
    with open(outfilename, "w") as outfile, proc_pool as pool:
        core.runTime(writeout=outfile);
        core.PWS("# Window size (-w):\t" + wsize_str + "kb", outfile);
        core.PWS("# Input directory:\t" + indir, outfile);
        core.PWS("# Missing threshold (seq):\t" + str(args.missing_seq_thresh), outfile);
        core.PWS("# Missing threshold (site):\t" + str(args.missing_site_thresh), outfile);
        core.PWS("# Gap threshold (seq):\t" + str(args.gap_seq_thresh), outfile);
        core.PWS("# Gap threshold (site):\t" + str(args.gap_site_thresh), outfile);
        core.PWS("# Gap threshold (window):\t" + str(args.gap_window_thresh), outfile);
        core.PWS("# Missing+gap threshold (seq):\t" + str(args.mg_seq_thresh), outfile);
        core.PWS("# Missing+gap threshold (site):\t" + str(args.mg_site_thresh), outfile);
        core.PWS("# Output file:\t" + outfilename, outfile);
        if args.outdir:
            core.PWS("# Seq outdir:\t" + args.outdir, outfile);
        core.PWS("# ----------------", outfile);
        # Run time info for logging.
        ###################

        # outfile.write("# HEADER INFO:\n");
        # outfile.write("# window:\t" + "Unique window ID (scaffold:start-end)\n");
        # outfile.write("# scaffold:\t" + "Scaffold ID of window\n");
        # outfile.write("# start:\t" + "The start coordinate of the window (0-based)\n");
        # outfile.write("# end:\t" + "The end coordinate of the window\n");
        # outfile.write("# total.seqs:\t" + "The number of sequences in the alignment (excluding ancestral seqs)\n");
        # outfile.write("# aln.len:\t" + "The length of the alignment, including indels\n");
        # outfile.write("# seqs.above.missing:\t" + "The number of sequences in the alignment with a percentage of sites above the specified missing data threshold (-m)\n");
        # outfile.write("# seqs.all.missing:\t" + "The number of sequences in the alignment with a all sites missing.\n");
        # outfile.write("# avg.seq.len.wo.missing:\t" + "The average length of all sequences in the alignment while excluding gaps or missing data (X, N)\n");
        # outfile.write("# seqs.above.gappy:\t" + "The number of sequences in the alignment with a percentage of sites above the specified gap threshold (-g)\n");
        # outfile.write("# seqs.all.gappy:\t" + "The number of sequences in the alignment that are all gaps\n");
        # outfile.write("# seqs.all.missing.gappy:\t" + "The number of sequences in the alignment that are all missing OR gaps\n");
        # outfile.write("# gappy.cols:\t" + "The number of alignment columns that exist in 5 bp windows that have 3 or more gaps in at least -g fraction of sequences in the alignment\n");
        # outfile.write("# uniq.seqs:\t" + "The number of unique sequences in the alignment\n");
        # outfile.write("# ident.seqs:\t" + "The number of identical sequences in the alignment\n");
        # outfile.write("# invariant.sites:\t" + "The number of sites in the alignment with only 1 allele\n");
        # outfile.write("# informative.sites:\t" + "The number of sites in the alignment that have at least 2 alleles that are present in 2 species\n");
        # outfile.write("# sites.w.missing:\t" + "The number of sites in the alignment with at least one missing data character (X, N)\n");
        # outfile.write("# percent.sites.with.missing:\t" + "The percent of sites in the alignment with at least one missing data character (X, N)\n");
        # outfile.write("# sites.high.missing:\t" + "The number of sites in the alignment that are made up of a high percentage of missing data characters (over -m)\n");
        # outfile.write("# percent.sites.high.missing:\t" + "The percent of sites in the alignment that are made up of a high percentage of missing data characters (over -m)\n");
        # outfile.write("# sites.all.missing:\t" + "The number of sites in the alignment that are made up of all missing data characters\n");
        # outfile.write("# sites.w.gap:\t" + "The number of sites in the alignment with at least one gap\n");
        # outfile.write("# percent.sites.with.gap:\t" + "The percent of sites in the alignment with at least one gap\n");
        # outfile.write("# sites.high.gap:\t" + "The number of sites in the alignment that are made up of a high percentage of gaps (over -g)\n");
        # outfile.write("# percent.sites.high.gap:\t" + "The percent of sites in the alignment that are made up of a high percentage of gaps (over -g)\n");
        # outfile.write("# sites.all.gap:\t" + "The number of sites in the alignment that are made up of all gaps\n");
        # outfile.write("# sites.missing.gap:\t" + "The number of sites in the alignment that are missing OR gaps\n");
        # outfile.write("# percent.sites.missing.gap:\t" + "The percent of sites in the alignment that are missing OR gaps\n");
        # outfile.write("# .filter:\t" + "Columns with the .filter suffix have the same definition as those defined above, but are counted after the window filter (gappy.cols filter)\n");
        # outfile.write("# written:\t" + "Whether or not the filtered sequence was written to a file\n");
        # outfile.write("# ----------------\n");

        outfile.write("\t".join(headers) + "\n");

        # Header info
        ###################

        scaff_ids = [ d for d in os.listdir(indir) ];

        #scaff_dirs = [ os.path.join(indir, d) for d in os.listdir(indir) ];
        # Get the list of all scaffold directories to loop over
        
        num_scaffs = len(scaff_ids);
        cur_scaff_num = 1;
        # Counting for status updates

        for result in pool.imap_unordered(processScaff, ((scaff_id, indir, cur_scaff_num, num_scaffs, args.outdir) for scaff_id in scaff_ids)):
            for outline in result:
                outfile.write("\t".join(outline) + "\n");
            cur_scaff_num += 1;
            # Increment scaffold counter
        ## End scaff loop
    ## Close output file


############################################################