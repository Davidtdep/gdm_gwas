# ==============================================================================
# TITLE: Gestational Diabetes GWAS Analysis
# DATE: July 2026
# ==============================================================================

# ==============================================================================
# 0. LIBRARIES
# ==============================================================================
# Data manipulation and cleaning
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(forcats)

# Data import/export
library(readxl)
library(writexl)
library(XML)

# Visualization
library(ggplot2)
library(circlize)
library(ggrepel)

# Bioinformatics and Literature APIs
library(rentrez)
library(rcrossref)
library(roadoi)

# Statistical / Utility tools
library(pcaMethods)
library(stringdist)

# ==============================================================================
# 1. CONFIGURATION & PATHS
# ==============================================================================
# Users should update these paths to match their local working directory
INPUT_DIR  <- "~/"
OUTPUT_DIR <- "~/"

# Example alternative if using an RStudio project structure:
# INPUT_DIR <- "data/" 
# OUTPUT_DIR <- "results/"

# ==============================================================================
# 2. DATA INPUT
# ==============================================================================
# Import raw GWAS data
data_file_path <- paste0(INPUT_DIR, "gwas.tsv")
data <- read.delim(data_file_path, header = TRUE, dec = ".")

# Import significant data (if already generated from a previous session)
# data_significant_path <- paste0(INPUT_DIR, "data_significant.xlsx")
# data_significant <- read_excel(data_significant_path)

# ==============================================================================
# 3. DATA CLEANING & PREPROCESSING
# ==============================================================================

# Filter all single-allele SNPs (e.g., rs12345-A)
data_filtered <- data %>%
  filter(str_detect(STRONGEST.SNP.RISK.ALLELE, "^rs\\d+-[A-Za-z]$"))

# Create a new "ALLELE" column by extracting the allele letter for each SNP
data_filtered$ALLELE <- sub(".*-", "", data_filtered$STRONGEST.SNP.RISK.ALLELE)

# ------------------------------------------------------------------------------
# Helper Function: extract_ci_bounds
# Purpose: Extracts lower and upper bounds of a confidence interval string.
# Input: A character string representing the CI (e.g., "[1.05-1.12]").
# Output: A numeric vector of length 2: c(lower_bound, upper_bound).
# ------------------------------------------------------------------------------
extract_ci_bounds <- function(ci_text) {
  # Handle empty, missing, or non-reported cases
  if (is.na(ci_text) || ci_text == "" || ci_text == "[NR]" || ci_text == "NR") {
    return(c(NA, NA))
  }
  
  # Updated regular expression to capture cases missing the closing bracket.
  # The \\]? makes the closing bracket optional.
  pattern <- "\\[([-0-9.]+)[-–]([-0-9.]+)\\]?"
  match <- regexec(pattern, ci_text)
  result <- regmatches(ci_text, match)
  
  if (length(result) > 0 && length(result[[1]]) >= 3) {
    lower <- as.numeric(result[[1]][2])
    upper <- as.numeric(result[[1]][3])
    return(c(lower, upper))
  }
  
  # Return NA if no match is found
  return(c(NA, NA))
}

# ------------------------------------------------------------------------------
# Wrapper Function: extract_confidence_intervals
# Purpose: Iterates through the dataset to append CI bounds as new columns.
# ------------------------------------------------------------------------------
extract_confidence_intervals <- function(data) {
  # Initialize columns for the lower and upper bounds
  data$CI_lower <- NA
  data$CI_upper <- NA
  
  # Apply the extraction function row by row
  for (i in 1:nrow(data)) {
    bounds <- extract_ci_bounds(data$X95..CI..TEXT.[i])
    data$CI_lower[i] <- bounds[1]
    data$CI_upper[i] <- bounds[2]
  }
  
  return(data)
}

# Apply the CI extraction function to the filtered dataset
data_filtered <- extract_confidence_intervals(data_filtered)

# Filter alleles with a confidence interval strictly greater than 1 or less than 1
data_significant <- subset(data_filtered, CI_lower > 1 | CI_upper < 1)

# Remove rows containing NAs in the confidence interval columns
data_significant <- data_significant %>%
  filter(!is.na(CI_lower) & !is.na(CI_upper))

# ==============================================================================
# 4. EXPORTS
# ==============================================================================

# Extract and export the list of unique significant SNPs
SNPs <- unique(data_significant$SNPS)

# Optional: Export the unique SNPs to a CSV file for downstream use
# write.csv(SNPs, file = paste0(OUTPUT_DIR, "unique_snps.csv"), row.names = FALSE)

# Optional: Export the full significant dataset
# write_xlsx(data_significant, paste0(OUTPUT_DIR, "Supplementary_Material_1.xlsx"))


# ==============================================================================
# 5. DATA MATCHING (ALLELE FREQUENCIES)
# ==============================================================================

# Load dataframe with allele frequencies
# Note: Update path to match your local INPUT_DIR if necessary
snps_fqs_path <- paste0(INPUT_DIR, "SNPs_FQs.csv")
SNPs_FQs <- read.csv(snps_fqs_path)

# Clean environment to free memory, keeping only essential dataframes
# (Note: Use caution with rm() in shared reproducible scripts)
rm(list = setdiff(ls(), c("data_significant", "SNPs_FQs")))

# Define target Colombian population groups
population_cols <- c("PLQ", "CHG", "ATQCES", "ATQPGC", "CLM")

# Initialize new columns for populations with NAs
data_significant[, population_cols] <- NA

# ------------------------------------------------------------------------------
# 5.1 Format Frequency Data
# ------------------------------------------------------------------------------
# Prepare the allele frequency table.
# Pivots to long format and selects exclusively the required columns to ensure 
# clean joining without transferring extra metadata.
SNPs_FQs_long <- SNPs_FQs %>%
  rename(SNPS = rs) %>%
  pivot_longer(
    cols = c(A, C, T, G),
    names_to = "ALLELE",
    values_to = "Frequency"
  ) %>%
  select(SNPS, Population.Group, ALLELE, Frequency)

# ------------------------------------------------------------------------------
# 5.2 Merge Frequencies into Significant Data
# ------------------------------------------------------------------------------
# Update 'data_significant' with matching population frequencies
data_significant <- data_significant %>%
  # Clean key columns to ensure a perfect join
  mutate(
    SNPS = trimws(SNPS),
    ALLELE = toupper(trimws(ALLELE))
  ) %>%
  # Pivot target population columns to long format
  pivot_longer(
    cols = all_of(population_cols),
    names_to = "Population.Group"
  ) %>%
  # Standardize the new population column to match the frequency table format
  mutate(
    Population.Group = toupper(trimws(Population.Group))
  ) %>%
  # Join with the filtered frequency table
  left_join(
    SNPs_FQs_long,
    by = c("SNPS", "ALLELE", "Population.Group")
  ) %>%
  # Reconstruct the dataframe back to wide format
  pivot_wider(
    names_from = "Population.Group",
    values_from = "Frequency"
  )

# Remove the redundant 'value' column generated during pivot
if("value" %in% names(data_significant)) {
  data_significant$value <- NULL
}

# Convert population columns to numeric
data_significant[, population_cols] <- lapply(
  data_significant[, population_cols], 
  function(x) as.numeric(as.character(x))
)

# ------------------------------------------------------------------------------
# 5.3 Filter and Calculate Group Means
# ------------------------------------------------------------------------------
# Filter out SNPs that lack established frequencies across ALL studied populations
# (e.g., reduces dataset by removing SNPs with NAs in all 5 population columns)
data_significant <- data_significant %>%
  filter(!if_all(all_of(population_cols), is.na))

# Create new columns calculating the mean for African and European origin groups
data_significant <- data_significant %>%
  mutate(
    African_mean = rowMeans(select(., all_of(c("PLQ", "CHG"))), na.rm = TRUE),
    European_mean = rowMeans(select(., all_of(c("ATQCES", "ATQPGC", "CLM"))), na.rm = TRUE)
  )

# ==============================================================================
# 6. STATISTICAL ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
# 6.1 Correlation Analysis
# ------------------------------------------------------------------------------

# Normality Test (Shapiro-Wilk) across populations
lapply(data_significant[population_cols], shapiro.test)

# Pearson Correlation Test
# Calculate the correlation matrix (using complete observations)
correlation_matrix <- data_significant %>%
  select(all_of(population_cols)) %>%
  cor(method = "pearson", use = "complete.obs")

# Prepare correlation data for ggplot2 plotting
rounded_corr_matrix <- round(correlation_matrix, 2)

# Isolate the lower triangle by setting the upper triangle to NA
rounded_corr_matrix[upper.tri(rounded_corr_matrix)] <- NA

# Melt the matrix into a long format for heatmap generation
long_corr <- melt(rounded_corr_matrix, na.rm = TRUE)

# Generate Correlation Heatmap
ggplot(data = long_corr, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black", size = 4) +
  scale_fill_gradient2(
    high = "#E41A1C", 
    low = "#377EB8", 
    mid = "white",
    midpoint = 0, 
    limit = c(-1, 1), 
    name = "Pearson\nCorrelation"
  ) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed()

# ------------------------------------------------------------------------------
# 6.2 Principal Component Analysis (PCA)
# ------------------------------------------------------------------------------

# Extract and transpose frequency data for PCA
freq_data <- t(data_significant[, population_cols])

# Perform PCA using probabilistic PCA (ppca)
pca_result <- pca(freq_data, method = "ppca", nPcs = 2)

# Convert the PCA result into a dataframe for visualization
scores <- as.data.frame(scores(pca_result))
scores$Population <- rownames(scores)

# Calculate percentage of variance explained by each principal component
var_explained <- pca_result@R2 * 100

# Generate PCA Plot
ggplot() +
  # Custom grid lines at 0.25 intervals
  geom_hline(yintercept = seq(-2, 2, by = 1), color = "gray90", linewidth = 0.2) +
  geom_vline(xintercept = seq(-2, 2, by = 1), color = "gray90", linewidth = 0.2) +
  
  # Dashed lines for origin axes (x=0, y=0)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Scatter points
  geom_point(data = scores, aes(x = PC1, y = PC2, color = Population), size = 4.5) +
  
  # Repelled text labels to prevent overlap
  geom_text_repel(
    data = scores, 
    aes(x = PC1, y = PC2, label = Population),
    size = 3, 
    box.padding = 0.45, 
    point.padding = 0.5, 
    segment.color = "black", 
    min.segment.length = 0
  ) +
  
  scale_color_brewer(palette = "Set1") +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black")
  )


# ==============================================================================
# 7. DATA VISUALIZATION
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Allele Frequency Jitter Plot (African vs. European Means)
# ------------------------------------------------------------------------------

# Step 1: Prepare data by grouping based on ancestry predominance (> 50%)
data_significant <- data_significant %>%
  mutate(
    group = case_when(
      African_mean > 0.5 & European_mean < 0.5 ~ "Group A",
      African_mean < 0.5 & European_mean > 0.5 ~ "Group B",
      TRUE                                     ~ "Other" # Remaining points
    )
  )

# Step 2: Generate publication-quality jitter plot
# Note: The mathematical thresholds used above are 0.5 (50%), but the legend 
# labels below state 65% / 35% based on original code. Adjust if necessary.
ggplot(data_significant, aes(x = African_mean, y = European_mean, color = group)) +
  
  # Add reference dashed lines at 0.5
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  
  # Add jittered points to prevent overplotting
  geom_jitter(width = 0.04, height = 0.04, alpha = 0.7, size = 2.5) +
  
  # Define a professional, colorblind-friendly palette
  scale_color_manual(
    values = c(
      "Group A" = "red",     # Safe Orange/Red
      "Group B" = "#0072B2", # Safe Blue
      "Other"   = "grey70"   # Neutral grey for the rest
    ),
    labels = c(
      "Group A" = "PLQ > 65% and European < 35%",
      "Group B" = "PLQ < 35% and European > 65%",
      "Other"   = "Other patterns"
    )
  ) +
  
  # Clear and concise axis labels
  labs(
    x = "African (San Basilio de Palenque)",
    y = "European"
  ) +
  
  # Apply a classic theme with fine-tuned adjustments for a professional look
  theme_classic(base_size = 12) + 
  theme(
    legend.position = "none",
    legend.title = element_blank(), 
    legend.text = element_text(size = 10),
    axis.text = element_text(color = "black", size = 10),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5) 
  )

# ------------------------------------------------------------------------------
# 7.2 Chromosome Distribution Bar Plot (p vs. q arms)
# ------------------------------------------------------------------------------

# Step 1: Prepare chromosome data
# Extract chromosome arms (p or q) from the REGION string
data_significant$Arm <- sub("^\\d+([pq]).*", "\\1", data_significant$REGION)

# Remove trailing decimals from chromosome IDs (e.g., "1.0" -> "1")
data_significant$CHR_ID <- gsub("\\.0*$", "", data_significant$CHR_ID)

# Create a numeric version of CHR_ID for proper sequential sorting
data_significant <- data_significant %>%
  mutate(CHR_ID_numeric = as.numeric(CHR_ID))

# Step 2: Generate horizontal bar plot
ggplot(data_significant, aes(x = reorder(CHR_ID, -CHR_ID_numeric), fill = Arm)) +
  
  geom_bar(width = 0.7) + 
  
  # Pastel color palette for chromosomal arms
  scale_fill_manual(
    name = "Arm",
    values = c("p" = "#FDBF6F", "q" = "#A6CEE3")
  ) +
  
  labs(
    x = "Chromosome",
    y = "N. of Associations"
  ) +
  
  # Professional minimalist theme with flipped coordinates
  theme_minimal(base_size = 12) +
  coord_flip() +
  theme(
    # Dashed lines for the associations axis (X-axis before flip)
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed", linewidth = 0.4),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 11),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# 7.3 Top 20 Associated Genes Bar Plot
# ------------------------------------------------------------------------------

# 1. Extract unique genes.
# NOTE: To concatenate both columns before finding unique values, use c():
# gene_vector <- unique(c(data_significant$MAPPED_GENE, data_significant$REPORTED.GENE.S.))
# The logic below is preserved exactly as originally written.
gene_vector <- unique(data_significant$MAPPED_GENE, data_significant$REPORTED.GENE.S.)

# 2. Split genes using a regular expression (splits by ", " OR by " - ").
# Returns a list where each element is a vector of genes.
split_genes_list <- strsplit(gene_vector, ",\\s*|\\s+-\\s+")

# 3. Flatten the list into a single vector of all genes.
all_genes_vector <- unlist(split_genes_list)

# 4. Count the frequency of each gene.
gene_counts <- table(all_genes_vector)

# 5. Convert the table to a dataframe for easier manipulation.
gene_counts_df <- as.data.frame(gene_counts)
names(gene_counts_df) <- c("Gene", "Count")

# 6. Sort the dataframe by count in descending order.
gene_counts_df_sorted <- gene_counts_df[order(gene_counts_df$Count, decreasing = TRUE), ]

# 7. Clean row names for neatness.
rownames(gene_counts_df_sorted) <- NULL

# 8. Subset the top 20 genes.
Top20_genes <- gene_counts_df_sorted[1:20, ]

# Generate horizontal bar plot for the top 20 genes
ggplot(Top20_genes, aes(x = Count, y = fct_reorder(Gene, Count))) +
  geom_col(fill = "#A6CEE3") +
  labs(
    y = "Gene", 
    x = "N. of Associations"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed", linewidth = 0.4)
  )


# ==============================================================================
# 8. COHORT DEMOGRAPHICS BY ANCESTRY
# ==============================================================================

# ------------------------------------------------------------------------------
# 8.1 Define Raw Study Demographics
# ------------------------------------------------------------------------------
# Create the dataset mapping study IDs to ancestry and case/control counts
results <- tribble(
  ~PUBMEDID, ~Ancestry,           ~Cases, ~Controls,
  
  # ID 22233651: Explicitly Korean cohort
  22233651,  "Korean",            468,    1242,
  
  # ID 35220425 (Pervjakova et al.): Complete data retrieved from the paper
  35220425,  "African American",  91,     985,
  35220425,  "East Asian",        867,    2984,
  35220425,  "European",          3780,   344902,
  35220425,  "Hispanic American", 174,    619,
  35220425,  "South Asian",       573,    3761,
  
  # ID 40325049: Using the "Superset" row (total Cases/Controls)
  # "Females" rows are ignored as they are a subset of the superset.
  40325049,  "Chinese",           12024,  67845
)

# Calculate total individuals per cohort
results$Total <- results$Cases + results$Controls

# ------------------------------------------------------------------------------
# 8.2 Group and Summarize by Unified Ancestry
# ------------------------------------------------------------------------------
ancestry_totals <- results %>%
  mutate(Ancestry_Unified = case_when(
    grepl("European", Ancestry, ignore.case = TRUE) ~ "European",
    # Aggregate Korean, Chinese, and East Asian into a single group
    grepl("Chinese|Korean|East Asian", Ancestry, ignore.case = TRUE) ~ "East Asian",
    grepl("South Asian", Ancestry, ignore.case = TRUE) ~ "South Asian",
    grepl("African", Ancestry, ignore.case = TRUE) ~ "African/African American",
    grepl("Hispanic", Ancestry, ignore.case = TRUE) ~ "Hispanic/Latino",
    TRUE ~ Ancestry
  )) %>%
  group_by(Ancestry_Unified) %>%
  summarise(
    Total_Cases = sum(Cases),
    Total_Controls = sum(Controls),
    Grand_Total = sum(Total)
  ) %>%
  arrange(desc(Grand_Total))

# Display the aggregated demographics
print(ancestry_totals)


# ==============================================================================
# 9. ANCESTRY DISTRIBUTION CHORD DIAGRAM
# ==============================================================================

# Sort the ancestry groups by the total number of individuals (descending)
ancestry_totals <- ancestry_totals %>% arrange(desc(Grand_Total))

# ------------------------------------------------------------------------------
# 9.1 Prepare Link Data for the Chord Diagram
# ------------------------------------------------------------------------------
links <- data.frame(
  # Each ancestry appears twice (once mapping to cases, once mapping to controls)
  from = rep(ancestry_totals$Ancestry_Unified, 2),  
  to = c(rep("Cases", nrow(ancestry_totals)), rep("Controls", nrow(ancestry_totals))),
  value = c(ancestry_totals$Total_Cases, ancestry_totals$Total_Controls)
)

# Remove connections with a value of zero to prevent rendering empty lines
links <- links[links$value > 0, ]

# ------------------------------------------------------------------------------
# 9.2 Define Ordering and Aesthetics
# ------------------------------------------------------------------------------

# Define manual factor levels for ancestry groups and destination sectors
order_manual <- c("European", "East Asian", "South Asian", "African/African American", "Hispanic/Latino")
order_cases_controls <- c("Cases", "Controls")

# Apply the desired sorting to the link dataframe
links$from <- factor(links$from, levels = order_manual)
links$to <- factor(links$to, levels = order_cases_controls)

# Define custom hex colors for each group
grid.col <- c(
  "European" = "#FF9966", 
  "East Asian" = "#006D6F", 
  "South Asian" = "#1F305E",
  "African/African American" = "#FF0800",
  "Hispanic/Latino" = "#B284BE",
  "Cases" = "black", 
  "Controls" = "darkgray"
)

# Reverse the order to display the ancestry sectors counter-clockwise
orden_inverso <- rev(order_manual)

# Update the 'from' column with the reversed order and sort the dataframe
links$from <- factor(links$from, levels = orden_inverso)
links <- links[order(links$from), ]

# Establish custom transparency for each ancestry group's chords
# (0 = completely opaque, 1 = completely transparent)
transparency_values <- c(
  "European" = 0.5, 
  "East Asian" = 0.5, 
  "South Asian" = 0.5,
  "African/African American" = 0.5,
  "Hispanic/Latino" = 0.5
)

# ------------------------------------------------------------------------------
# 9.3 Render the Chord Diagram
# ------------------------------------------------------------------------------

# Clear any pre-existing circos plots to prevent overlap/rendering bugs
circos.clear()

# Generate the chord diagram
chordDiagram(
  links,
  grid.col = grid.col, 
  transparency = transparency_values[as.character(links$from)], 
  annotationTrack = "grid",               # Adds grid lines for the sectors
  preAllocateTracks = list(track.height = 0.1) # Pre-allocates track height for labels
)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================