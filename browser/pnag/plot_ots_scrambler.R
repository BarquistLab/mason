library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(viridis)
library(writexl)

# Source shared utilities
source("./pnag/r_utils.R")

# Load data
path_output <- commandArgs(trailingOnly = TRUE)[1]
screen <- commandArgs(trailingOnly = TRUE)[2]

output_df <- read_csv(paste0(path_output, "/result_table.csv"))
df_plot <- read_csv(paste0(path_output, "/df_plot.csv"))

# Rename columns to match mason convention (Python uses hyphens/spaces)
df_plot <- df_plot %>% rename(
  off_target_type = `off-target type`,
  nr_mismatches = num_mismatch
)


# Generate kableExtra HTML table (no color formatting, just styled)
table_out <- kable(output_df, format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed")) %>%
  column_spec(1, bold = TRUE)

table_out %>% save_kable(paste0(path_output, "/result_table.html"))

# Re-save result_table as CSV and Excel from R-loaded data
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)
write_xlsx(output_df, paste0(path_output, "/result_table.xlsx"))


# Create bar plots using shared functions
df_plot$counts <- as.numeric(df_plot$counts)

# Make off_target_type a factor
ot_levels <- c("OT in transcriptome", "OT in TIR regions")
if ("OT in HMP microbiome" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in HMP microbiome")
}
if ("OT in human transcriptome" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in human transcriptome")
}
if ("OT in TIR of essential genes" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in TIR of essential genes")
}
df_plot$off_target_type <- factor(df_plot$off_target_type, levels = ot_levels)

# Set ASO factor order matching output_df
df_plot$ASO <- factor(df_plot$ASO, levels = unique(output_df$ASO))

wplot <- nrow(output_df) + 5

plot_ot_dodged(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
               "Number of off-targets in TIR regions", "inferno",
               output_prefix = paste0(path_output, "/plot_ots_start_regions"),
               wplot = wplot)

plot_ot_dodged(df_plot[df_plot$off_target_type == "OT in transcriptome", ],
               "Number of off-targets in whole transcriptome", "viridis",
               output_prefix = paste0(path_output, "/plot_ots_whole_transcriptome"),
               wplot = nrow(output_df) + 6)

# Screening off-target plots
if (!is.na(screen) && screen == "microbiome" && "OT in HMP microbiome" %in% df_plot$off_target_type) {
  df_hmp <- df_plot[df_plot$off_target_type == "OT in HMP microbiome", ]
  df_hmp$ASO <- factor(df_hmp$ASO, levels = unique(output_df$ASO))
  plot_ot_single(df_hmp, "Number of off-targets in HMP microbiome (0 mismatches)", "#440154",
                 output_prefix = paste0(path_output, "/plot_ots_hmp"),
                 wplot = wplot)
}

if (!is.na(screen) && screen == "human" && "OT in human transcriptome" %in% df_plot$off_target_type) {
  df_human <- df_plot[df_plot$off_target_type == "OT in human transcriptome", ]
  df_human$ASO <- factor(df_human$ASO, levels = unique(output_df$ASO))
  plot_ot_single(df_human, "Number of off-targets in human transcriptome (0 mms)", "#31688e",
                 output_prefix = paste0(path_output, "/plot_ots_human"),
                 wplot = wplot)
}

if (!is.na(screen) && screen == "essential_genes" && "OT in TIR of essential genes" %in% df_plot$off_target_type) {
  df_ess <- df_plot[df_plot$off_target_type == "OT in TIR of essential genes", ]
  df_ess$ASO <- factor(df_ess$ASO, levels = unique(output_df$ASO))
  plot_ot_dodged(df_ess, "Off-targets in TIR of essential genes", "inferno",
                 output_prefix = paste0(path_output, "/plot_ots_essential_genes"),
                 wplot = wplot)
}
