library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(viridis)
library(writexl)



print("before making tabls")
# Load data
path_output <- commandArgs(trailingOnly = TRUE)[1]
screen <- commandArgs(trailingOnly = TRUE)[2]
#path_output <- "./browser/pnag/static/data/2026_02_10_10_44_13/b0185/outputs"

output_df <- read_csv(paste0(path_output, "/result_table.csv"))
df_plot <- read_csv(paste0(path_output, "/df_plot.csv"))
ml_res <- read_csv(paste0(path_output, "/saved_table_ml.csv"))


output_df$MIC_pred <- ml_res$MIC_pred
# rank MIC_pred from low to high
output_df$MIC_ranked <- rank(output_df$MIC_pred, ties.method = "min")

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
                                                               "yellow"))) %>%
  # make ranked MIC column colored by rank. use coloring green to yellow to red, with 100 breaks, and reverse the order so that low MIC is green and high MIC is red
    column_spec(which(names(output_df) == "MIC_ranked"), color = "black", background = colorRampPalette(c("green", "yellow", "red"))(100)[as.numeric(cut(output_df$MIC_ranked, breaks = 100))])



print("before saving html")
# save table as html


table_out %>% save_kable(paste0(path_output, "/result_table.html"))


print("before saving csv")
# Save dataframe as CSV
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)
print("after")
write_xlsx(output_df, paste0(path_output, "/result_table.xlsx"))
print("after saving xlsx")

# Create bar plot using ggplot2
df_plot$counts <- as.numeric(df_plot$counts)

# make df_plot$off_target_type a factor and order it
ot_levels <- c("OT in transcriptome", "OT in TIR regions")
if ("OT in HMP microbiome" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in HMP microbiome")
}
if ("OT in human genome" %in% df_plot$off_target_type) {
  ot_levels <- c(ot_levels, "OT in human genome")
}
df_plot$off_target_type <- factor(df_plot$off_target_type, levels = ot_levels)

# plot for only ots in TIR regions
df_plot$ASO <- factor(gsub(".*_(ASO_\\d+)", "\\1", df_plot$ASO), levels = unique(output_df$ASO))

p <- ggplot(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3)))) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Number of off-targets in TIR regions") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
# define title of legend
    guides(fill = guide_legend(title = "# mismatches")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=13),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size= 25, hjust = 0.5, face = "bold"),
        # add increased legend fontsize and adjust position (topright of plot inside plot area
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="white", linewidth = .5),
        legend.direction = "horizontal",
        legend.title = element_text(size=15),
        # make legend above plot, not inside plot area
        legend.position = "top",
        legend.margin = margin(6, 10, 6, 6)) +
  # reverse viridis
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "inferno") +
  # add counts to bars
    geom_text(aes(label = counts), position = position_dodge(width = 0.9), vjust = -0.25, size = 5)


# make width as much as nr ow rows in output_df
wplot <- nrow(output_df) + 5

ggsave(paste0(path_output, "/plot_ots_start_regions.png"), p, width = wplot,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_start_regions.svg"), p, width = wplot,
       limitsize = FALSE)

# do plot for ots in whole transcriptome
p_whole <- ggplot(df_plot[df_plot$off_target_type == "OT in transcriptome", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3))) ) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
  ggtitle("Number of off-targets in whole transcriptome") +
  guides(fill = guide_legend(title = "# off-targets")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=13),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        # put title in the middle of the plot and make bigger and bold
        plot.title = element_text(size= 25, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill=alpha('white', 0.7), linewidth = .5),
        legend.title = element_text(size=15),
        # make legend above plot, not inside plot area
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
  df_hmp$ASO <- factor(gsub(".*_(ASO_\\d+)", "\\1", df_hmp$ASO), levels = unique(output_df$ASO))

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
  df_human$ASO <- factor(gsub(".*_(ASO_\\d+)", "\\1", df_human$ASO), levels = unique(output_df$ASO))

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