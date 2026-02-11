
Sys.setenv(HOME = "/home/jakob")


library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(rmelting)
library(viridis)
library(writexl)

# Source shared utilities
source("./pnag/r_utils.R")

path_output <- commandArgs(trailingOnly = TRUE)[1]
screen <- commandArgs(trailingOnly = TRUE)[2]


# Load data
ot_table <- paste0(path_output, "/offtargets_fulltranscripts_sorted.tab")

print("ot_table")
print(ot_table)

all_off_targets <- read_table(ot_table, col_names = TRUE)



print(head(all_off_targets))

# change probe_id colname to ASO_id
colnames(all_off_targets)[6] <- "ASO"
colnames(all_off_targets)[4] <- "position_from_CDS_start"


# Filter mismatches with >+7 cons. matches. use dplyr::filter
all_off_targets <- all_off_targets %>% filter(longest_stretch > 6) %>%
  # add TIR column to dataframe
    mutate(TIR = ifelse(position_from_CDS_start %in% -20:5, "TIR", "not in TIR")) %>%
  # sort by ASO (asc), then TIR (desc),then longest_stretch (desc), then position_from_CDS_start (asc)
    arrange(ASO, TIR, desc(longest_stretch), position_from_CDS_start)

print(head(all_off_targets))


# Create output table
write_csv(all_off_targets, file = paste0(path_output, "/offtargets_fulltranscripts_sorted.csv"))
write_csv(all_off_targets[all_off_targets$TIR == "TIR", ],
          file = paste0(path_output, "/offtargets_startregions_sorted.csv"))

# write excel file
library(writexl)
write_xlsx(all_off_targets, paste0(path_output, "/offtargets_fulltranscripts_sorted.xlsx"))



output_df <- read.table(paste0(path_output, "/result_table.tsv"), sep = "\t", header = TRUE)


# Count off-targets using shared function (index_by_name=FALSE for checker)
ot_results <- count_offtargets_by_mismatch(all_off_targets, output_df, index_by_name = FALSE)
output_df <- ot_results$output_df
df_plot <- ot_results$df_plot

## Microbiome off-targets (0mm exact matches against HMP database)
if (!is.na(screen) && screen == "microbiome") {
  hmp_file <- paste0(path_output, "/offtargets_microbiome_sorted.tab")
  if (file.exists(hmp_file)) {
    hmp_ot <- read_tsv(hmp_file, col_names = TRUE)
    unique_asos <- unique(hmp_ot$probe_id)
    output_df[["OT_HMP_0mm"]] <- 0
    for (aso_n in unique_asos) {
      ot_aso <- hmp_ot[hmp_ot$probe_id == aso_n, ]
      n_0mm <- sum(ot_aso$num_mismatch == 0)
      output_df[output_df$ASO == aso_n, "OT_HMP_0mm"] <- n_0mm
      df_plot <- rbind(df_plot, data.frame(
        ASO = aso_n,
        off_target_type = "OT in HMP microbiome",
        transcripts = "HMP microbiome",
        counts = n_0mm,
        target.sequence = ot_aso$probe_seq[1],
        nr_mismatches = 0,
        stringsAsFactors = FALSE
      ))
    }
    hmp_ot_clean <- hmp_ot %>% mutate(ASO = probe_id) %>% select(-probe_id)
    write_csv(hmp_ot_clean, file = paste0(path_output, "/offtargets_hmp_sorted.csv"))
    write_xlsx(hmp_ot_clean, paste0(path_output, "/offtargets_hmp_sorted.xlsx"))
  }
}

## Human genome off-targets (0mm exact matches against GRCh38)
if (!is.na(screen) && screen == "human") {
  human_file <- paste0(path_output, "/offtargets_human_sorted.tab")
  if (file.exists(human_file)) {
    human_ot <- read_tsv(human_file, col_names = TRUE)
    unique_asos <- unique(human_ot$probe_id)
    output_df[["OT_GRCh38_0mm"]] <- 0
    for (aso_n in unique_asos) {
      ot_aso <- human_ot[human_ot$probe_id == aso_n, ]
      n_0mm <- sum(ot_aso$num_mismatch == 0)
      output_df[output_df$ASO == aso_n, "OT_GRCh38_0mm"] <- n_0mm
      df_plot <- rbind(df_plot, data.frame(
        ASO = aso_n,
        off_target_type = "OT in human genome",
        transcripts = "human genome",
        counts = n_0mm,
        target.sequence = ot_aso$probe_seq[1],
        nr_mismatches = 0,
        stringsAsFactors = FALSE
      ))
    }
    human_ot_clean <- human_ot %>% mutate(ASO = probe_id) %>% select(-probe_id)
    write_csv(human_ot_clean, file = paste0(path_output, "/offtargets_human_sorted.csv"))
    write_xlsx(human_ot_clean, paste0(path_output, "/offtargets_human_sorted.xlsx"))
  }
}

print(output_df)
print(df_plot)
# remove location column if no target gene was given


# calculate melting temperature for each ASO
output_df[["Tm (°C)"]] <- sapply(output_df$ASO_seq, function(x) {
  melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`
})

# calculate Mw for each ASO
output_df[["Mw"]] <- sapply(output_df$ASO_seq, calculate_pna_mw)


output_df[["%_SC_bases"]] <- round(output_df$SC_bases / nchar(output_df$ASO_seq) * 100, 2)

# re-arrange output df by ASO 	ASO_seq	SC_bases  % SC bases  Tm (°C)	pur_perc 	long_pur_stretch 	OT_TIR_0mm 	OT_TIR_1mm
# OT_TIR_2mm 	OT_TIR_3mm 	OT_tot_0mm 	OT_tot_1mm 	OT_tot_2mm 	OT_tot_3mm
select_cols <- c("ASO", "ASO_seq", "SC_bases", "%_SC_bases", "Tm (°C)", "pur_perc", "long_pur_stretch",
                 "Mw",
                 "OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm", "OT_tot_0mm", "OT_tot_1mm",
                 "OT_tot_2mm", "OT_tot_3mm")
if ("OT_HMP_0mm" %in% names(output_df)) select_cols <- c(select_cols, "OT_HMP_0mm")
if ("OT_GRCh38_0mm" %in% names(output_df)) select_cols <- c(select_cols, "OT_GRCh38_0mm")
output_df <- output_df %>% select(all_of(select_cols))



print("before making tabls")
# make table with kableextra
table_out <- kable(output_df, format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("bordered","hover", "condensed"),) %>%
  column_spec(1, bold = TRUE) %>%
  # if column 3 is <5, make it red
  column_spec(3:4, color = "black",
              background = ifelse(output_df[["%_SC_bases"]] > 60, "red",
                                                      ifelse(output_df[["%_SC_bases"]] > 50, "salmon",
                                                       ifelse(output_df[["%_SC_bases"]] > 40, "yellow",
                                                        ifelse(output_df[["%_SC_bases"]] > 34, "yellow",
                                                               "lightgreen"))))) %>%
    column_spec(5, color = "black", background = ifelse(output_df[["Tm (°C)"]] <= 35, "red",
                                                            ifelse(output_df[["Tm (°C)"]] <= 40, "yellow",
                                                                       "lightgreen"))) %>%
    column_spec(6, color = "black", background = ifelse(output_df$pur_perc < 30, "lightgreen",
                                                        ifelse(output_df$pur_perc < 51, "white",
                                                               "yellow")))


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
df_plot$ASO <- factor(df_plot$ASO, levels = unique(output_df$ASO))

p <- ggplot(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3)))) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Number of off-targets in TIR regions") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
# define title of legend
    guides(fill = guide_legend(title = "# mismatches")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#, size=13),
        #axis.text.y = element_text(size=15),
        #axis.title = element_text(size=20),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        # add increased legend fontsize and adjust position (topright of plot inside plot area
        #legend.text = element_text(size=12),
        legend.background = element_rect(fill="white", linewidth = .5),
        legend.direction = "horizontal",
        #legend.title = element_text(size=12),
        # make legend above plot, not inside plot area
        legend.position = "top",
        legend.margin = margin(6, 10, 6, 6)) +
  # reverse viridis
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "inferno") +
  # add counts to bars
    geom_text(aes(label = counts), position = position_dodge(width = 0.9), vjust = -0.25)


ggsave(paste0(path_output, "/plot_ots_start_regions.png"), p,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_start_regions.svg"), p,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_start_regions.pdf"), p,
       limitsize = FALSE)

# do plot for ots in whole transcriptome
p_whole <- ggplot(df_plot[df_plot$off_target_type == "OT in transcriptome", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3))) ) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
  ggtitle("Number of off-targets in whole transcriptome") +
  guides(fill = guide_legend(title = "# mismatches")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        # put title in the middle of the plot and make bigger and bold
        plot.title = element_text(hjust = 0.5, face = "bold"),
       # legend.text = element_text(size=12),
        legend.background = element_rect(fill=alpha('white', 0.7), linewidth = .5),
       # legend.title = element_text(size=12),
        # make legend above plot, not inside plot area
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(6, 10, 6, 6)) +
  scale_fill_viridis(discrete = TRUE, direction = -1, option = "viridis") +
    geom_text(aes(label = counts), position = position_dodge(width = 0.9), vjust = -0.25)



wplot <- nrow(output_df) + 6


ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.png"), p_whole,
       limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.svg"), p_whole,
         limitsize = FALSE)
ggsave(paste0(path_output, "/plot_ots_whole_transcriptome.pdf"), p_whole,
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_text(aes(label = counts), vjust = -0.25)

  ggsave(paste0(path_output, "/plot_ots_hmp.png"), p_hmp, limitsize = FALSE)
  ggsave(paste0(path_output, "/plot_ots_hmp.svg"), p_hmp, limitsize = FALSE)
  ggsave(paste0(path_output, "/plot_ots_hmp.pdf"), p_hmp, limitsize = FALSE)
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_text(aes(label = counts), vjust = -0.25)

  ggsave(paste0(path_output, "/plot_ots_human.png"), p_human, limitsize = FALSE)
  ggsave(paste0(path_output, "/plot_ots_human.svg"), p_human, limitsize = FALSE)
  ggsave(paste0(path_output, "/plot_ots_human.pdf"), p_human, limitsize = FALSE)
}

# Remove temporary files
#file.remove(paste0(path_output, "/result_table.tsv", sep = ""))
