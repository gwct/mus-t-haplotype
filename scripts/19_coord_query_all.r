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

parseCoords <- function(feature, st){






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

new_tree = T
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

# if(!new_tree){
#   logfile = paste(logfile, flank_interval_mb, "mb-no-pahari-all.log", sep="")
# }else{
#   logfile = paste(logfile, flank_interval_mb, "mb-new-tree-no-pahari-all.log", sep="")
# }

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
}

############################################################
# Main

split_windows = split(windows, f=windows$chr)
# Split by chromosome.

for(chr_windows in split_windows){
  chrome = as.character(chr_windows[1,]$chr)
  chrome_len = chromes[chromes$chr==chrome,]$chr.len
  
  cat("#", as.character(Sys.time()), "| Starting", chrome, "with", nrow(chr_windows), "features\n")

  cat("#", as.character(Sys.time()), "| Reading species tree\n")

  if(!new_tree){
    species_tree_file = here("analysis", "02-mus-t-windows", "04-iqtree-no-pahari", "chr17", paste0(window_size_kb, "kb"), "concat", paste0(chrome, "-concat.treefile.rooted"))
  }else{
    species_tree_file = here("analysis", "02-mus-t-windows-new-tree", "04-iqtree-no-pahari", "chr17", paste0(window_size_kb, "kb"), "concat", paste0(chrome, "-concat.treefile.rooted"))
  }

  species_tree = read.tree(file=species_tree_file)
  # Read the species tree

  if(!new_tree){
    outfile = here("data", "coord-query-no-pahari", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".dists.all"))
  }else{
    outfile = here("data", "coord-query-new-tree-no-pahari", paste0(feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".dists.all"))
  }
  
  ##########

  cat("#", as.character(Sys.time()), "|", chrome, "| Calculating tree dists\n")

  window_dists = data.frame("chr"=c(), "start"=c(), "end"=c(), "window"=c(), "filt"=c(),
                "rf.st"=c(), "wrf.st"=c(), "spr.st"=c(), "kf.st"=c(), "path.st"=c())
  # Initialize the df for all windows

  for(i in 1:nrow(chr_windows)){
    cur_window = chr_windows[i,]

    cur_dists = data.frame("chr"=cur_window$chr, "start"=cur_window$start, "end"=cur_window$end, "window"=cur_window$window, "filt"="FILTER",
                "rf.st"=NA, "wrf.st"=NA, "spr.st"=NA, "kf.st"=NA, "path.st"=NA)
    # Initialize the df for this window

    if(!isFiltered(cur_window)){
        cur_tree = read.tree(text=as.character(cur_window$unparsed.tree))
        # Get the tree for the current window

        cur_dists$filt = "PASS"
        # Add that the current window passed filtering

        cur_dists$rf.st = RF.dist(cur_tree, species_tree)
        cur_dists$wrf.st = wRF.dist(cur_tree, species_tree)
        spr = SPR.dist(cur_tree, species_tree)
        names(spr) = NULL
        cur_dists$spr.st = spr
        cur_dists$kf.st = KF.dist(cur_tree, species_tree)
        cur_dists$path.st = path.dist(cur_tree, species_tree)
        # Calculate distances from the target window to the species tree
    }

    window_dists = rbind(window_dists, cur_dists)
  }

  ##########

  window_dists = window_dists %>% arrange(start)

  cat("#", as.character(Sys.time()), "|", chrome, "| Writing tree dists to file:", outfile, "\n")
  write.csv(window_dists, file=outfile, row.names=F)
}

############################################################




