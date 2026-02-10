library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(viridis)
library(writexl)


# Load data
path_output <- commandArgs(trailingOnly = TRUE)[1]

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
df_plot$off_target_type <- factor(df_plot$off_target_type,
                                  levels = c("OT in transcriptome", "OT in TIR regions"))

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
