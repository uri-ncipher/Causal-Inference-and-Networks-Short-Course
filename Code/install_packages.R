# install_packages.R
# Run this script once to install all required packages.


# --- Helper: install only if not already installed ---
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste0("'", pkg, "' is already installed. Skipping."))
  }
}

# --- Required packages ---
packages <- c("rstudioapi", "lme4", "plyr", "igraph", "numDeriv", "gtools", 
              "doParallel", "foreach", "knitr", "kableExtra", "ggplot2", "forcats", 
              "geepack", "tidyverse", "cobalt", "scales", "readxl", 
              "purrr", "ergm", "sna", "network", "coda", "dplyr")


# --- Install ---
invisible(lapply(packages, install_if_missing))

message("Done! All packages are installed.")





