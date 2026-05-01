# install_packages.R

# Set CRAN mirror
options(repos = c(CRAN = "http://cran.rstudio.com/"))

# Install BiocManager (needed even if some Bioc packages are pre-installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager")
}
library(BiocManager)

# List of required packages NOT expected in rocker/shiny-verse base
# Assumes shiny, tidyverse, dplyr, ggplot2, stringr, data.table, reshape2 are present
pkgs <- c(
 # 'lme4', 'pbkrtest', 'car', 'rstatix', 'ragg', # Dependencies likely covered by shiny-verse
 # 'tidyverse', # Included in shiny-verse
 "ggpubr", "plumber", # Keep ggpubr, plumber might not be included
 # 'rstudioapi', # Likely included
 # 'dplyr', 'stringr', 'ggplot2', # Included in shiny-verse
 # 'grid', # Base R
 "ggrepel", "gridGraphics",
 # 'shiny', # Included in shiny-verse
 "shinyFiles", "bslib",
 "spsComps", "DT", "shinyjs", "shinycssloaders",
 "data.table",
 "plotly", "ggiraph", "writexl", "fontawesome",
 "fst", "R.utils",
 "qtl2", # Bioconductor package
 "reshape2", # Added missing package,
 "qtl2fst", # Added missing package
 "promises", "future", "future.apply", # V2: async rendering support + prewarm
 "shinybusy", "waiter", # V2: Modal spinner support
 # V3: S3-backed data layer and Posit Connect deployment
 "paws.storage", # AWS S3 client (minimal paws slice)
 "filelock",     # Cross-process cache index locking
 "rsconnect", "renv" # Posit Connect publish + dependency pinning
)

# Install packages using BiocManager
print("--- Installing additional required packages ---")
BiocManager::install(pkgs, update = FALSE, ask = FALSE)

# Final check (only check packages we explicitly tried to install)
print("--- Final package check ---")
installed_pkgs <- installed.packages()[, "Package"]
required_pkgs_unique <- unique(pkgs)
missing_pkgs <- required_pkgs_unique[!(required_pkgs_unique %in% installed_pkgs)]

if (length(missing_pkgs) > 0) {
 stop(paste("Failed to install required packages:", paste(missing_pkgs, collapse = ", ")))
} else {
 print("--- Additional required packages installed successfully ---")
}
