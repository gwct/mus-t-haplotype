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

getDists <- function(direction, adj, target_window, target_tree, bed.id, int_start, st){
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
      
      result$rf = RF.dist(target_tree, query_tree)
      result$wrf = wRF.dist(target_tree, query_tree)
      spr = SPR.dist(target_tree, query_tree)
      names(spr) = NULL
      result$spr = spr
      result$kf = KF.dist(target_tree, query_tree)
      result$path = path.dist(target_tree, query_tree)
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

  cur_window = subset(windows, chr == feature$chr & start <= feature$int.start & end > feature$int.start)

  feature_dists = data.frame("chr"=cur_window$chr, "start"=cur_window$start, "end"=cur_window$end, 
                "bed.id"=feature$bed.id, "int.start"=feature$int.start, 
                "window"=cur_window$window, "adj"=0, "filt"="FILTER",
                "rf"=NA, "wrf"=NA, "spr"=NA, "kf"=NA, "path"=NA,
                "rf.st"=NA, "wrf.st"=NA, "spr.st"=NA, "kf.st"=NA, "path.st"=NA)
  # Initialize the df for this feature with the intersecting target window
  
  # print(paste("    START WINDOW", feature$gid, do.call(paste, c(cur_window, sep="    "))))
  # print(paste("    START FILTER:", isFiltered(cur_window)))
  
  if(!isFiltered(cur_window)){
    cur_tree = read.tree(text=as.character(cur_window$unparsed.tree))
    # Get the tree for the current window

    feature_dists$filt = "PASS"
    # Add that the target window passed filtering

    feature_dists$rf.st = RF.dist(cur_tree, st)
    feature_dists$wrf.st = wRF.dist(cur_tree, st)
    spr = SPR.dist(cur_tree, st)
    names(spr) = NULL
    feature_dists$spr.st = spr
    feature_dists$kf.st = KF.dist(cur_tree, st)
    feature_dists$path.st = path.dist(cur_tree, st)
    # Calculate distances from the target window to the species tree

    cur_adj = 1
    while(cur_adj <= windows_per_interval){
      #print(paste(feature$gid, "-", cur_adj, sep=""))
      
      for_df = getDists("forward", cur_adj, cur_window, cur_tree, feature$bed.id, feature$int.start, st)
      feature_dists = rbind(feature_dists, for_df)
      
      back_df = getDists("backward", (-1*cur_adj), cur_window, cur_tree, feature$bed.id, feature$int.start, st)
      feature_dists = rbind(feature_dists, back_df)
      cur_adj = cur_adj + 1
    }
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

flank_interval_mb = 5
flank_interval = flank_interval_mb * 1000000
# The distance to calculate outward from both sides of a given coordinate in Mb

windows_per_interval = flank_interval / window_size

test_run = FALSE

serial = FALSE

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
  infile = here("data", "inv-breaks.bed")
  logfile = "logs/coord-query-inv-bps-"
  cols = c("chr", "start", "end", "bed.id")
}
logfile = paste(logfile, flank_interval_mb, "mb.log", sep="")

if(read_data){
  cat("#", as.character(Sys.time()), "| Reading input data\n")
  window_file = here("data", paste0("mm10-", window_size_kb, "kb-topo-counts-no-pahari.csv"))
  windows = read_csv(window_file, comment="#")
  names(windows) = make.names(names(windows))
  chromes = unique(select(windows, chr, chr.len))
  # The window data with trees.
  
  coords = read.table(infile)
  names(coords) = cols
  coords$chr = factor(coords$chr)

  coords = coords[coords$chr %in% windows$chr,]
  coords$int.start = round((coords$start + coords$end) / 2)
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
  species_tree_file = here("analysis", "02-mus-t-windows", "04-iqtree", "chr17", paste0(window_size_kb, "kb"), "concat", paste0(chrome, "-concat-no-pahari.treefile.rooted"))
  species_tree = read.tree(file=species_tree_file)
  # Read the species tree

  outfile = here("data", "coord-query", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".dists"))
  outfile_random = here("data", "coord-query", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".random"))
  
  ##########

  if(gen_random){
    cat("#", as.character(Sys.time()), "|", chrome, "| Selecting random windows\n")

    null_coords = data.frame("chr"=c(), "start"=c(), "end"=c(), "bed.id"=c())
    for(j in 1:nrow(chr_coords)){
      rand_pos = sample(1:chrome_len, 1)
      null_coords = rbind(null_coords, data.frame("chr"=chrome, "start"=rand_pos, "end"=rand_pos, "bed.id"=NA))
    }
    null_coords$chr = as.character(null_coords$chr)
    null_coords$int.start = null_coords$start
  }
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

  if(do_random){
    cat("#", as.character(Sys.time()), "|", chrome, "| Calculating null dists\n")
    null_results <- foreach(i=1:nrow(null_coords), .packages=c("ape", "phangorn")) %dopar% {
      feature = null_coords[i,]
      cat("#", "FEATURE", i, do.call(paste, c(feature, sep="    ")), "\n")
      parseCoords(feature, species_tree)
    } 
    null_results = do.call(rbind, null_results)
    cat("#", as.character(Sys.time()), "|", chrome, "| Writing null dists to file:", outfile_random, "\n")
    write.csv(null_results, file=outfile_random, row.names=F)
  }
  if(!test_run && !serial){
    stopCluster(cl)
  }
}

############################################################




