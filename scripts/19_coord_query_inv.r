############################################################
# For penn genomes, 02.2021
# Given a set of coordinates, calculates RF from a given 
# distance around those coordinates and compares to random 
# coordinates (or a given set of null coordinates).
# Gregg Thomas
############################################################

library(tidyverse)
library(ape)
library(phangorn)
library(doParallel)
library(here)

#############################################################
# Functions

isFiltered <- function(w){
  filter_flag = FALSE
  
  if(w$filter != "PASS"){
    filter_flag = TRUE
  }  

  return(filter_flag)
}

######################

checkBounds <- function(direction, adj_start, c_len){
  bounds_flag = TRUE
  if(direction=="backward" && adj_start < 0){
    bounds_flag = FALSE
  }
  
  if(direction=="forward" && adj_start > c_len){
    bounds_flag = FALSE
  }
  
  return(bounds_flag)
}

######################

# getRandWindow <- function(window_df){
#   rand_window = data.frame("window"="", "chr"="", "start"="", "tree"="")
  
#   random_window = as.character(sample(window_df$window, 1))
#   rand_window$window = random_window
  
#   window_data = window_df[window_df$window==random_window,]
#   rand_window$chr = window_data$chr
#   rand_window$start = window_data$start
#   rand_window$tree = as.character(window_data$unparsed.tree)
  
#   return(rand_window)
# }

######################

getDists <- function(direction, adj, target_window, target_tree, target_filt, bed.id, int_start, st){

  result = data.frame("chr"=target_window$chr, "start"=target_window$start, "end"=target_window$end, 
                      "bed.id"=bed.id, "int.start"=int_start, 
                      "window"=NA, "adj"=adj, "filt"="FILTER",
                      "rf"=NA, "wrf"=NA, "spr"=NA, "kf"=NA, "path"=NA,
                      "rf.st"=NA, "wrf.st"=NA, "spr.st"=NA, "kf.st"=NA, "path.st"=NA)
  
  adj_start = target_window$start + (adj * window_size)
  
  within_bounds = checkBounds(direction, adj_start, target_window$chr.len)
  # Will be TRUE if window is beyond bounds of chromosome length.
  
  if(within_bounds){
    query_window = subset(windows, chr==target_window$chr & start==adj_start)
    result$window = query_window$window
    # Get the next window based on the current adjaceny step and window size
    
    win_filtered = isFiltered(query_window)
    # Will be TRUE if any of the filters of the query window are unpassed.

    if(!win_filtered){
      result$filt = "PASS"

      query_tree = read.tree(text=as.character(query_window$unparsed.tree))
      # Get the tree for the adjacent window
      
      if(target_filt == "PASS"){
        result$rf = RF.dist(target_tree, query_tree)
        result$wrf = wRF.dist(target_tree, query_tree)
        spr = SPR.dist(target_tree, query_tree)
        names(spr) = NULL
        result$spr = spr
        result$kf = KF.dist(target_tree, query_tree)
        result$path = path.dist(target_tree, query_tree)
      }
      # Calculate the wRF for the two trees

      result$rf.st = RF.dist(query_tree, st)
      result$wrf.st = wRF.dist(query_tree, st)
      spr = SPR.dist(query_tree, st)
      names(spr) = NULL
      result$spr.st = spr
      result$kf.st = KF.dist(query_tree, st)
      result$path.st = path.dist(query_tree, st)
      # Calculate the wRF for the query window to the species tree 
    }
  }
  
  return(result)
}

######################

parseCoords <- function(feature, st){

  cur_windows = subset(windows, chr == feature$chr & start >= feature$int.start & end <= feature$int.end)

  feature_dists = data.frame("chr"=c(), "start"=c(), "end"=c(), "bed.id"=c(),
                "window"=c(), "filt"=c(),
                "rf.st"=c(), "wrf.st"=c(), "spr.st"=c(), "kf.st"=c(), "path.st"=c())
  # Initialize the df for all windows

  for(i in 1:nrow(cur_windows)){
    cur_window = cur_windows[i,]

    cur_dists = data.frame("chr"=cur_window$chr, "start"=cur_window$start, "end"=cur_window$end, "bed.id"=feature$bed.id, 
                "window"=cur_window$window, "filt"="FILTER",
                "rf.st"=NA, "wrf.st"=NA, "spr.st"=NA, "kf.st"=NA, "path.st"=NA)
    # Initialize the df for this window

    if(!isFiltered(cur_window)){
        cur_tree = read.tree(text=as.character(cur_window$unparsed.tree))
        # Get the tree for the current window

        cur_dists$filt = "PASS"
        # Add that the current window passed filtering

        cur_dists$rf.st = RF.dist(cur_tree, st)
        cur_dists$wrf.st = wRF.dist(cur_tree, st)
        spr = SPR.dist(cur_tree, st)
        names(spr) = NULL
        cur_dists$spr.st = spr
        cur_dists$kf.st = KF.dist(cur_tree, st)
        cur_dists$path.st = path.dist(cur_tree, st)
        # Calculate distances from the target window to the species tree
    }

    feature_dists = rbind(feature_dists, cur_dists)

  }

  return(feature_dists)
}

############################################################
# Options

server = T
# Set if running on server

window_size_kb = 10
# Window size in kb

window_size = window_size_kb * 1000
# Window size in bp

num_cores = 10
# Number of cores

gen_random = T
# Set to generate the random windows

do_random = T
# Set to do the random windows

read_data = T
# Set depending on whether or not data has been read into environment

feature = "inv-bps"
# One of inv-bps

new_tree = F
# Whether or not to use the new tree topology

flank_interval_mb = 5
flank_interval = flank_interval_mb * 1000000
# The distance to calculate outward from both sides of a given coordinate in Mb

windows_per_interval = flank_interval / window_size

test_run = FALSE

serial = TRUE

if(!server){
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
  num_cores = 1
}

if(serial){
  num_cores = 1
}

############################################################
# Input

cat("#", as.character(Sys.time()), "| Feature:", feature, "\n")
if(feature=="inv-bps"){
  infile = here("data", "inversions.bed")
  logfile = "logs/coord-query-inv-bps-"
  cols = c("chr", "start", "end", "bed.id")
}

if(!new_tree){
  logfile = paste(logfile, flank_interval_mb, "mb-no-pahari.log", sep="")
}else{
  logfile = paste(logfile, flank_interval_mb, "mb-new-tree-no-pahari.log", sep="")
}

if(read_data){
  cat("#", as.character(Sys.time()), "| Reading input data\n")

  if(!new_tree){
    window_file = here("data", paste0("mm10-", window_size_kb, "kb-topo-counts-no-pahari.csv"))
  }else{
    window_file = here("data", paste0("mm10-", window_size_kb, "kb-topo-counts-new-tree-no-pahari.csv"))
  }

  windows = read_csv(window_file, comment="#")
  names(windows) = make.names(names(windows))
  chromes = unique(select(windows, chr, chr.len))
  # The window data with trees.
  
  coords = read.table(infile)
  names(coords) = cols
  coords$chr = factor(coords$chr)

  coords = coords[coords$chr %in% windows$chr,]
  coords$int.start = coords$start - flank_interval
  coords$int.end = coords$end + flank_interval
  # The coordinates to test in BED format.
}

############################################################
# Main

split_coords = split(coords, f=coords$chr)
# Split by chromosome.

for(chr_coords in split_coords){
  chrome = as.character(chr_coords[1,]$chr)
  chrome_len = chromes[chromes$chr==chrome,]$chr.len
  
  cat("#", as.character(Sys.time()), "| Starting", chrome, "with", nrow(chr_coords), "features\n")

  cat("#", as.character(Sys.time()), "| Reading species tree\n")

  if(!new_tree){
    species_tree_file = here("analysis", "02-mus-t-windows", "04-iqtree-no-pahari", "chr17", paste0(window_size_kb, "kb"), "concat", paste0(chrome, "-concat.treefile.rooted"))
  }else{
    species_tree_file = here("analysis", "02-mus-t-windows-new-tree", "04-iqtree-no-pahari", "chr17", paste0(window_size_kb, "kb"), "concat", paste0(chrome, "-concat.treefile.rooted"))
  }

  species_tree = read.tree(file=species_tree_file)
  # Read the species tree

  if(!new_tree){
    outfile = here("data", "coord-query-no-pahari", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".dists.inv"))
    #outfile_random = here("data", "coord-query-no-pahari", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".random"))
  }else{
    outfile = here("data", "coord-query-new-tree-no-pahari", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".dists.inv"))
    #outfile_random = here("data", "coord-query-new-tree-no-pahari", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".random"))
  }
  
  ##########

#   if(gen_random){
#     cat("#", as.character(Sys.time()), "|", chrome, "| Selecting random windows\n")

#     null_coords = data.frame("chr"=c(), "start"=c(), "end"=c(), "bed.id"=c())
#     for(j in 1:nrow(chr_coords)){
#       rand_pos = sample(1:chrome_len, 1)
#       null_coords = rbind(null_coords, data.frame("chr"=chrome, "start"=rand_pos, "end"=rand_pos, "bed.id"=NA))
#     }
#     null_coords$chr = as.character(null_coords$chr)
#     null_coords$int.start = null_coords$start
#   }
  # Select the null/random coordinates to compare against.

  ##########

  cat("#", as.character(Sys.time()), "|", chrome, "| Calculating tree dists\n")

  if(serial){
    #print("serial")
    results = NULL
    for(i in 1:nrow(chr_coords)){
      feature = chr_coords[i,]
      #print(feature)
      feature_dists = parseCoords(feature, species_tree)
      results = rbind(results, feature_dists)
    }
  }else{
    #print("parallel")
    cl <- parallel::makeCluster(num_cores, setup_strategy="sequential", outfile=logfile)
    registerDoParallel(cl)
    results <- foreach(i=1:nrow(chr_coords), .packages=c("ape", "phangorn")) %dopar% {
      feature = chr_coords[i,]
      cat("#", "FEATURE", i, do.call(paste, c(feature, sep="    ")), "\n")
      parseCoords(feature, species_tree)
    } 
    results = do.call(rbind, results)
  }

  ##########

  cat("#", as.character(Sys.time()), "|", chrome, "| Writing tree dists to file:", outfile, "\n")
  write.csv(results, file=outfile, row.names=F)
  
  ##########

#   if(do_random){
#     cat("#", as.character(Sys.time()), "|", chrome, "| Calculating null dists\n")
#     null_results <- foreach(i=1:nrow(null_coords), .packages=c("ape", "phangorn")) %dopar% {
#       feature = null_coords[i,]
#       cat("#", "FEATURE", i, do.call(paste, c(feature, sep="    ")), "\n")
#       parseCoords(feature, species_tree)
#     } 
#     null_results = do.call(rbind, null_results)
#     cat("#", as.character(Sys.time()), "|", chrome, "| Writing null dists to file:", outfile_random, "\n")
#     write.csv(null_results, file=outfile_random, row.names=F)
#   }
  if(!test_run && !serial){
    stopCluster(cl)
  }
}

############################################################




