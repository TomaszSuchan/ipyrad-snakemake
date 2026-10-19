#!/usr/bin/env Rscript
# Script to run whoa (Where's my Heterozygotes at?) QC analysis
# This evaluates genotyping accuracy by examining heterozygote miscall rates

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(vcfR)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--vcf"), type = "character", default = NULL,
              help = "Input VCF file path", metavar = "character"),
  make_option(c("--posteriors"), type = "character", default = NULL,
              help = "Output path for posteriors RDS file", metavar = "character"),
  make_option(c("--plot"), type = "character", default = NULL,
              help = "Output path for plot PDF", metavar = "character"),
  make_option(c("--report"), type = "character", default = NULL,
              help = "Output path for text report", metavar = "character"),
  make_option(c("--min-bin"), type = "integer", default = 1000,
              help = "Minimum bin size parameter for whoa [default=%default]", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$vcf) || is.null(opt$posteriors) || is.null(opt$plot) || is.null(opt$report)) {
  print_help(opt_parser)
  stop("All output file paths must be specified.\n", call. = FALSE)
}

# Install whoa if not already installed
if (!requireNamespace("whoa", quietly = TRUE)) {
  message("Installing whoa package from GitHub...")
  suppressPackageStartupMessages(library(devtools))
  devtools::install_github("eriqande/whoa", quiet = TRUE, upgrade = "never")
}

# Load whoa
suppressPackageStartupMessages(library(whoa))

# Read VCF file
message("Reading VCF file: ", opt$vcf)
vcf <- read.vcfR(opt$vcf, verbose = FALSE)

# Extract genotype matrix and convert to whoa format
message("Extracting genotypes and computing read depths...")
gt_matrix <- extract.gt(vcf, element = "GT")
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

# Check if we have depth information
if (all(is.na(dp_matrix))) {
    stop("No depth (DP) information found in VCF. Whoa requires read depth data.")
}

# Compute expected and observed genotype frequencies for visualization
message("Computing genotype frequencies...")
geno_freqs <- exp_and_obs_geno_freqs(vcf)

# Infer heterozygote miscall rates
# NOTE: infer_m takes the vcfR object directly, not geno_freqs
message("Inferring heterozygote miscall rates (this may take a while)...")
posteriors_result <- infer_m(vcf, minBin = opt$`min-bin`)
posteriors <- posteriors_result$m_posteriors

# Save posteriors as RDS
message("Saving posteriors to: ", opt$posteriors)
saveRDS(posteriors, file = opt$posteriors)

# Generate plot
message("Generating plot...")
pdf(opt$plot, width = 10, height = 6)
posteriors_plot(posteriors)
dev.off()

# Generate text report
message("Writing report to: ", opt$report)
sink(opt$report)
cat(strrep("=", 80), "\n")
cat("WHOA: Genotyping Accuracy QC Report\n")
cat(strrep("=", 80), "\n\n")

cat("Input VCF:", opt$vcf, "\n")
cat("minBin parameter:", opt$`min-bin`, "\n\n")

cat(strrep("=", 80), "\n")
cat("Summary Statistics\n")
cat(strrep("=", 80), "\n\n")

# Print summary table
print(posteriors)

cat("\nWHOA_POSTERIORS_TSV\n")
if (!is.null(posteriors) && nrow(posteriors) > 0) {
    write.table(posteriors, file=stdout(), sep="\t", row.names=FALSE, quote=FALSE)
}
cat("\n")

cat("\n")
cat(strrep("=", 80), "\n")
cat("Interpretation:\n")
cat(strrep("=", 80), "\n")
cat("- mean_m: Mean heterozygote miscall rate estimate\n")
cat("- lo95: Lower 95% credible interval\n")
cat("- hi95: Upper 95% credible interval\n")
cat("- n_bin: Number of genotypes in this read depth bin\n")
cat("- mean_dp: Mean read depth for genotypes in this bin\n\n")

cat("Higher miscall rates at lower read depths indicate that genotypes\n")
cat("with fewer reads are less reliable. This information can guide\n")
cat("filtering decisions for downstream analyses.\n")

sink()

message("Whoa analysis completed successfully!")
