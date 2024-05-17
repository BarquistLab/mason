library(ggplot2)
library(dplyr)
library(ggpubr)

print("hello world")

path_output <- commandArgs(trailingOnly = TRUE)[1]
  #"./data/sandra_rn4220_2/outputs"
target <- commandArgs(trailingOnly = TRUE)[2]
         # ""
# Load data
ot_table <- paste0(path_output, "/offtargets_fulltranscripts_sorted.tab")
           #"./data/sandra_rn4220_2/outputs/offtargets_fulltranscripts_sorted.tab"


all_off_targets <- read.table(ot_table, sep = "\t", header = TRUE)


# change probe_id colname to ASO_id
colnames(all_off_targets)[6] <- "ASO_id"
colnames(all_off_targets)[4] <- "position_from_CDS_start"

# Filter mismatches with >+7 cons. matches. use dplyr::filter
all_off_targets <- all_off_targets %>% filter(longest_stretch > 6) %>%
  # add TIR column to dataframe
    mutate(TIR = ifelse(position_from_CDS_start %in% -20:5, "TIR", "not in TIR")) %>%
  # sort by ASO_id (asc), then TIR (desc),then longest_stretch (desc), then position_from_CDS_start (asc)
    arrange(ASO_id, TIR, desc(longest_stretch), position_from_CDS_start)


# Create output table
write.csv(all_off_targets, file = paste0(path_output, "/offtargets_fulltranscripts_sorted.csv"), row.names = FALSE)
write.csv(all_off_targets[all_off_targets$TIR == "TIR", ],
          file = paste0(path_output, "/offtargets_startregions_sorted.csv"), row.names = FALSE)



output_df <- read.table(paste0(path_output, "/result_table.tsv"), sep = "\t", header = TRUE)
# read.table("./data/sandra_rn4220_2/outputs/result_table.tsv", sep = "\t", header = TRUE)

# Add columns for df_plot
df_plot <- data.frame(ASO = character(0), "off_target_type" = character(0),
                      "transcripts" = character(0), counts = numeric(0),
                      "target_sequence" = character(0))

for (i in unique(all_off_targets$ASO)) {
  target_seq <- all_off_targets[all_off_targets$ASO == i, ]$probe_seq[1]
  aso_n <-  gsub(".*_(ASO.*)", "\\1", i)
  ot_aso <- all_off_targets[all_off_targets$ASO == i, ]
  num_tot_ot <- nrow(ot_aso) - ifelse(target == "", 0, 1)
  num_tir_ot <- nrow(ot_aso[ot_aso$TIR == "TIR", ]) - ifelse(target == "", 0, 1)
  df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                        "transcripts" = "whole transcriptome", counts = num_tot_ot,
                                        "target sequence" = target_seq))
  df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                        "transcripts" = "start regions", counts = num_tir_ot,
                                        "target sequence" = target_seq))
  # assign values to output_df
  print(aso_n)
  print(i)
  output_df[output_df$ASO == aso_n, "OT_tot"] <- num_tot_ot
  output_df[output_df$ASO == aso_n, "OT_TIR"] <- num_tir_ot
}

print(output_df)
print(df_plot)
# remove location column if no target gene was given

if(target == "") {output_df$location <- NULL}

table_out <- ggtexttable(output_df, rows = NULL, theme = ttheme("light"))

hplot <- ifelse(nrow(output_df) >49, 50, nrow(output_df))
ggsave(paste0(path_output, "/result_table.png"), table_out, width = 11, height= hplot,
       limitsize = FALSE)
ggsave(paste0(path_output, "/result_table.pdf"), table_out, width = 11, height= hplot,
       limitsize = FALSE)


# Save dataframe as CSV
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)


# Create bar plot using ggplot2
df_plot$counts <- as.numeric(df_plot$counts)

# make df_plot$off_target_type a factor and order it
df_plot$off_target_type <- factor(df_plot$off_target_type, levels = c("OT in transcriptome", "OT in TIR regions"))

p <- ggplot(df_plot, aes(x = ASO, y = counts, fill = off_target_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
# define title of legend
    guides(fill = guide_legend(title = "Off-target type")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=13),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        # add increased legend fontsize and adjust position (topright of plot inside plot area
        legend.text = element_text(size=12),
        legend.background = element_rect(fill="white", linewidth = .5),
        legend.title = element_text(size=12),
        legend.position = c(0.2, 0.9),
        legend.margin = margin(6, 10, 6, 6)) +
  scale_fill_manual(values = c("OT in TIR regions" = "skyblue","OT in transcriptome" = "steelblue" ))+
  # add counts to bars
    geom_text(aes(label = counts), position = position_dodge(width = 0.9), vjust = -0.25, size = 5)
p

wplot <- ifelse(nrow(df_plot) >49, 50, 4 + nrow(df_plot)/3)
ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.png"), p, width = wplot,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.pdf"), p, width = wplot,
       limitsize = FALSE)



# Remove temporary files
file.remove(paste0(path_output, "/result_table.tsv", sep = ""))
