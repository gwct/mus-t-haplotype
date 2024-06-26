############################################################
# For turtles web, 07.21
# This generates the file "trees.html"
############################################################

cat("Rendering trees.rmd/html\n")
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)
library(here)
setwd(here("docs", "scripts", "generators"))    
print(getwd())
output_dir = "../.."
render("../markdown/trees.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)

cat("Rendering trees_no_pahari.rmd/html\n")
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)
library(here)
setwd(here("docs", "scripts", "generators"))    
print(getwd())
output_dir = "../.."
render("../markdown/trees_no_pahari.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)

cat("Rendering trees_no_pahari_new_tree.rmd/html\n")
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
library(rmarkdown)
library(here)
setwd(here("docs", "scripts", "generators"))    
print(getwd())
output_dir = "../.."
render("../markdown/trees_no_pahari_new_tree.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)