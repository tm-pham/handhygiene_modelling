# ============================================================================ #
# Packages
# ============================================================================ #
# Packages necessary
# .libPaths("/hpc/local/CentOS7/julius_id/R_libs/4.0.1/")
check.and.install.pkgs <- function(pkgs){
  new.packages <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
  if(length(new.packages)) install.packages(new.packages, 
                                            dependencies = TRUE, 
                                            # lib = "/hpc/local/CentOS7/julius_id/R_libs/4.0.1/",
                                            repos = "http://cran.us.r-project.org")
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}
check.and.install.pkgs(c("ggplot2",
                         "reshape2", 
                         "dplyr",
                         "tidyr",
                         "gridExtra")) 
