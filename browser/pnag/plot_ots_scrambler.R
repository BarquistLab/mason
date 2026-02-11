library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(viridis)
library(writexl)


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


# Create bar plots using ggplot2 (matching mason style)
df_plot$counts <- as.numeric(df_plot$counts)

# Make off_target_type a factor
ot_levels <- c("OT in transcriptome", "OT in TIR regions")
if ("OT in HMP microbiome" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in HMP microbiome")
}
if ("OT in human genome" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in human genome")
}
df_plot$off_target_type <- factor(df_plot$off_target_type, levels = ot_levels)

# Set ASO factor order matching output_df
df_plot$ASO <- factor(df_plot$ASO, levels = unique(output_df$ASO))

# Plot for OTs in TIR regions
p <- ggplot(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3)))) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Number of off-targets in TIR regions") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
  guides(fill = guide_legend(title = "# mismatches")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill = "white", linewidth = .5),
        legend.direction = "horizontal",
        legend.title = element_text(size = 15),
        legend.position = "top",
        legend.margin = margin(6, 10, 6, 6)) +
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "inferno") +
  geom_text(aes(label = counts), position = position_dodge(width = 0.9), vjust = -0.25, size = 5)

wplot <- nrow(output_df) + 5

ggsave(paste0(path_output, "/plot_ots_start_regions.png"), p, width = wplot,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_start_regions.svg"), p, width = wplot,
       limitsize = FALSE)

# Plot for OTs in whole transcriptome
p_whole <- ggplot(df_plot[df_plot$off_target_type == "OT in transcriptome", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3)))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
  ggtitle("Number of off-targets in whole transcriptome") +
  guides(fill = guide_legend(title = "# mismatches")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill = alpha('white', 0.7), linewidth = .5),
        legend.title = element_text(size = 15),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(6, 10, 6, 6)) +
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "viridis") +
  geom_text(aes(label = counts), position = position_dodge(width = 0.9), vjust = -0.25, size = 5)

wplot <- nrow(output_df) + 6

ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.png"), p_whole, width = wplot,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.svg"), p_whole, width = wplot,
       limitsize = FALSE)

# Microbiome off-target plot
if (!is.na(screen) && screen == "microbiome" && "OT in HMP microbiome" %in% df_plot$off_target_type) {
  df_hmp <- df_plot[df_plot$off_target_type == "OT in HMP microbiome", ]
  df_hmp$counts <- as.numeric(df_hmp$counts)
  df_hmp$ASO <- factor(df_hmp$ASO, levels = unique(output_df$ASO))

  p_hmp <- ggplot(df_hmp, aes(x = ASO, y = counts)) +
    geom_bar(stat = "identity", fill = "#440154") +
    ggtitle("Number of off-targets in HMP microbiome (0 mismatches)") +
    labs(x = "ASO sequence", y = "Number of off-targets") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 25, hjust = 0.5, face = "bold")) +
    geom_text(aes(label = counts), vjust = -0.25, size = 7)

  wplot_hmp <- nrow(output_df) + 5
  ggsave(paste0(path_output, "/plot_ots_hmp.png"), p_hmp, width = wplot_hmp, limitsize = FALSE)
  ggsave(paste0(path_output, "/plot_ots_hmp.svg"), p_hmp, width = wplot_hmp, limitsize = FALSE)
}

# Human genome off-target plot
if (!is.na(screen) && screen == "human" && "OT in human genome" %in% df_plot$off_target_type) {
  df_human <- df_plot[df_plot$off_target_type == "OT in human genome", ]
  df_human$counts <- as.numeric(df_human$counts)
  df_human$ASO <- factor(df_human$ASO, levels = unique(output_df$ASO))

  p_human <- ggplot(df_human, aes(x = ASO, y = counts)) +
    geom_bar(stat = "identity", fill = "#31688e") +
    ggtitle("Number of off-targets in human genome (0 mismatches)") +
    labs(x = "ASO sequence", y = "Number of off-targets") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 25, hjust = 0.5, face = "bold")) +
    geom_text(aes(label = counts), vjust = -0.25, size = 7)

  wplot_human <- nrow(output_df) + 5
  ggsave(paste0(path_output, "/plot_ots_human.png"), p_human, width = wplot_human, limitsize = FALSE)
  ggsave(paste0(path_output, "/plot_ots_human.svg"), p_human, width = wplot_human, limitsize = FALSE)
}
