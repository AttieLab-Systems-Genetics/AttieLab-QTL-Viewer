# scripts/convert_pheno_to_fst.R

# This script converts the large phenotype CSV file to a more efficient FST format.
# FST (Fast Serialization of Data Frames) allows for quick, on-demand column access
# without loading the entire file into memory. This is a one-time operation.

# Load required libraries
if (!requireNamespace("data.table", quietly = TRUE)) {
 install.packages("data.table")
}
if (!requireNamespace("fst", quietly = TRUE)) {
 install.packages("fst")
}
library(data.table)
library(fst)

# Configuration. Override with QTLAPP_DATA_ROOT on your host.
data_dir <- Sys.getenv("QTLAPP_DATA_ROOT", unset = "<your-data-root>")
csv_file_path <- file.path(data_dir, "pheno_with_covar_for_plot.csv")
fst_file_path <- file.path(data_dir, "pheno_with_covar_for_plot.fst")

# Conversion Process
if (!file.exists(csv_file_path)) {
 stop(paste("Source CSV file not found at:", csv_file_path))
}

message("Reading large CSV file...")
pheno_data <- data.table::fread(csv_file_path)
message(paste("CSV file loaded with", nrow(pheno_data), "rows and", ncol(pheno_data), "columns."))

message("Converting to FST format...")
fst::write_fst(pheno_data, fst_file_path, compress = 100)

message(paste("Successfully converted CSV to FST at:", fst_file_path))
message("The application will now use this optimized file for profile plots.")

# Show file size comparison
csv_size <- file.size(csv_file_path) / (1024^2) # MB
fst_size <- file.size(fst_file_path) / (1024^2) # MB
message(sprintf("File size comparison: CSV %.1f MB -> FST %.1f MB", csv_size, fst_size))
