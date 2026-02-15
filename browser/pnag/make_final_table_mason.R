library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(viridis)
library(writexl)

# Source shared utilities
source("./pnag/r_utils.R")

print("before making tabls")
# Load data
path_output <- commandArgs(trailingOnly = TRUE)[1]
screen <- commandArgs(trailingOnly = TRUE)[2]
use_ml <- commandArgs(trailingOnly = TRUE)[3]

output_df <- read_csv(paste0(path_output, "/result_table.csv"))
df_plot <- read_csv(paste0(path_output, "/df_plot.csv"))

if (!is.na(use_ml) && use_ml == "yes") {
  ml_res <- read_csv(paste0(path_output, "/saved_table_ml.csv"))
  output_df$MIC_pred <- ml_res$MIC_pred
  # rank MIC_pred from low to high
  output_df$MIC_ranked <- rank(output_df$MIC_pred, ties.method = "min")
}

# reorder columns so that MIC columns are next to each other. [1] "ASO"              "gene"             "ASO_seq"          "SC_bases"
#  [5] "%_SC_bases"       "Tm"               "pur_perc"         "long_pur_stretch"
#  [9] "Mw"               "OT_TIR_0mm"       "OT_TIR_1mm"       "OT_TIR_2mm"
# [13] "OT_TIR_3mm"       "OT_tot_0mm"       "OT_tot_1mm"       "OT_tot_2mm"
# [17] "OT_tot_3mm"       "MIC_pred"         "MIC_ranked"

#output_df <- output_df[, c(1:9, 14, 15, 10:13, 16:18)]

# make table with kableextra
table_out <- kable(output_df, format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("bordered","hover", "condensed", "centered")) %>%
  column_spec(1, bold = TRUE) %>%
  # if column 4 is <5, make it red
  column_spec(4:5, color = "black",
              background = ifelse(output_df[["%_SC_bases"]] > 60, "red",
                                                      ifelse(output_df[["%_SC_bases"]] > 50, "salmon",
                                                       ifelse(output_df[["%_SC_bases"]] > 40, "yellow",
                                                        ifelse(output_df[["%_SC_bases"]] > 34, "yellow",
                                                               "lightgreen"))))) %>%
    column_spec(6, color = "black", background = ifelse(output_df[["Tm"]] <= 35, "red",
                                                            ifelse(output_df[["Tm"]] <= 40, "yellow",
                                                                       "lightgreen"))) %>%
    column_spec(7, color = "black", background = ifelse(output_df$pur_perc < 30, "lightgreen",
                                                        ifelse(output_df$pur_perc < 51, "white",
                                                               "yellow")))

# Conditionally add MIC_ranked coloring
if ("MIC_ranked" %in% names(output_df)) {
  table_out <- table_out %>%
    column_spec(which(names(output_df) == "MIC_ranked"), color = "black",
                background = colorRampPalette(c("green", "yellow", "red"))(100)[as.numeric(cut(output_df$MIC_ranked, breaks = 100))])
}



print("before saving html")
# save table as html


table_out %>% save_kable(paste0(path_output, "/result_table.html"))


print("before saving csv")
# Save dataframe as CSV
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)
print("after")
write_xlsx(output_df, paste0(path_output, "/result_table.xlsx"))
print("after saving xlsx")

# Create bar plots using shared functions
df_plot$counts <- as.numeric(df_plot$counts)

# make df_plot$off_target_type a factor and order it
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

aso_short <- function(x) gsub(".*_(ASO_\\d+)", "\\1", x)
df_plot$ASO <- factor(aso_short(df_plot$ASO), levels = unique(aso_short(output_df$ASO)))

wplot <- nrow(output_df) + 5

plot_ot_dodged(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
               "Number of off-targets in TIR regions", "inferno",
               output_prefix = paste0(path_output, "/plot_ots_start_regions"),
               wplot = wplot)

plot_ot_dodged(df_plot[df_plot$off_target_type == "OT in transcriptome", ],
               "Number of off-targets in whole transcriptome", "viridis",
               legend_title = "# off-targets",
               output_prefix = paste0(path_output, "/plot_ots_whole_transcriptome"),
               wplot = nrow(output_df) + 6)

# Screening off-target plots
if (!is.na(screen) && screen == "microbiome" && "OT in HMP microbiome" %in% df_plot$off_target_type) {
  df_hmp <- df_plot[df_plot$off_target_type == "OT in HMP microbiome", ]
  df_hmp$ASO <- factor(df_hmp$ASO, levels = unique(aso_short(output_df$ASO)))
  plot_ot_single(df_hmp, "Number of off-targets in HMP microbiome (0 mismatches)", "#440154",
                 output_prefix = paste0(path_output, "/plot_ots_hmp"),
                 wplot = wplot)
}

if (!is.na(screen) && screen == "human" && "OT in human transcriptome" %in% df_plot$off_target_type) {
  df_human <- df_plot[df_plot$off_target_type == "OT in human transcriptome", ]
  df_human$ASO <- factor(df_human$ASO, levels = unique(aso_short(output_df$ASO)))
  plot_ot_single(df_human, "Number of off-targets in human transcriptome (0 mms)", "#31688e",
                 output_prefix = paste0(path_output, "/plot_ots_human"),
                 wplot = wplot)
}

if (!is.na(screen) && screen == "essential_genes" && "OT in TIR of essential genes" %in% df_plot$off_target_type) {
  df_ess <- df_plot[df_plot$off_target_type == "OT in TIR of essential genes", ]
  df_ess$ASO <- factor(df_ess$ASO, levels = unique(aso_short(output_df$ASO)))
  plot_ot_dodged(df_ess, "Off-targets in TIR of essential genes", "inferno",
                 output_prefix = paste0(path_output, "/plot_ots_essential_genes"),
                 wplot = wplot)
}