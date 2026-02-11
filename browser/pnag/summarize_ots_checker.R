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
target_gene <- commandArgs(trailingOnly = TRUE)[3]
use_ml <- commandArgs(trailingOnly = TRUE)[4]


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
write_xlsx(all_off_targets, paste0(path_output, "/offtargets_fulltranscripts_sorted.xlsx"))



output_df <- read.table(paste0(path_output, "/result_table.tsv"), sep = "\t", header = TRUE)


# Count off-targets using shared function (index_by_name=FALSE for checker)
ot_results <- count_offtargets_by_mismatch(all_off_targets, output_df, index_by_name = FALSE)
output_df <- ot_results$output_df
df_plot <- ot_results$df_plot

## Microbiome/human off-targets (0mm exact matches)
if (!is.na(screen) && screen %in% c("microbiome", "human")) {
  screen_results <- process_screen_offtargets(screen, path_output, output_df, df_plot,
                                               index_by_name = FALSE)
  output_df <- screen_results$output_df
  df_plot <- screen_results$df_plot
}

print(output_df)
print(df_plot)


# calculate melting temperature for each ASO
output_df[["Tm"]] <- sapply(output_df$ASO_seq, function(x) {
  melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`
})

# calculate Mw for each ASO
output_df[["Mw"]] <- sapply(output_df$ASO_seq, calculate_pna_mw)


output_df[["%_SC_bases"]] <- round(output_df$SC_bases / nchar(output_df$ASO_seq) * 100, 2)


# Branch: target gene mode vs standard checker mode
if (!is.na(target_gene) && nchar(target_gene) > 0) {
  # === TARGET MODE: MASON-style output ===

  # Generate Tm barplot (like melting.R)
  tm_df <- output_df
  tm_df$ASO_label <- tm_df$ASO

  g <- tm_df %>% ggplot(aes(x = ASO_label, y = Tm)) +
    geom_bar(stat = "identity", fill = "#31688e") +
    ggtitle("Predicted melting temperature (Tm)") +
    labs(x = "ASO sequence", y = paste0("Tm (\u00B0C)")) +
    theme_classic() +
    coord_cartesian(ylim = c(min(tm_df$Tm) - 5, max(tm_df$Tm) + 5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 25, hjust = 0.5, face = "bold")) +
    geom_text(aes(label = round(Tm, 1)), vjust = -0.25, size = 7)

  wplot <- nrow(tm_df) + 5
  ggsave(paste0(path_output, "/tm.png"), g, width = wplot, limitsize = FALSE)
  ggsave(paste0(path_output, "/tm.svg"), g, width = wplot, limitsize = FALSE)

  # Add gene column
  output_df[["gene"]] <- target_gene

  # Select columns in MASON order
  base_cols <- c("ASO", "gene", "ASO_seq", "SC_bases", "%_SC_bases", "Tm", "pur_perc", "long_pur_stretch",
                 "Mw", "OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm",
                 "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm")
  if ("OT_HMP_0mm" %in% names(output_df)) base_cols <- c(base_cols, "OT_HMP_0mm")
  if ("OT_GRCh38_0mm" %in% names(output_df)) base_cols <- c(base_cols, "OT_GRCh38_0mm")
  output_df <- output_df %>% select(all_of(base_cols))

  # Prepare saved_table_ml.csv (CAI, MFE features for ML)
  cai <- read.csv(paste0(path_output, "/cai_value.txt"), header = FALSE) %>% as.numeric() %>% unlist()
  mfe <- read.csv(paste0(path_output, "/mfe_values.txt"), header = FALSE)[1,1] %>% as.numeric() %>% unlist()

  saved_table_ml <- output_df %>%
    mutate(CAI = cai) %>%
    mutate(MFE = mfe) %>%
    select(ASO, Tm, `%_SC_bases`, CAI, OT_TIR_1mm, MFE, pur_perc)

  colnames(saved_table_ml) <- c("ASO", "Tm", "sc_bases", "CAI", "upec_tir_off_targets_1mm",
                                "MFE_UPEC", "purine_percentage")
  write.csv(saved_table_ml, file = paste0(path_output, "/saved_table_ml.csv"), row.names = FALSE)

  # Write result_table.csv and df_plot.csv for make_final_table_mason.R
  write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)
  write.csv(df_plot, file = paste0(path_output, "/df_plot.csv"), row.names = FALSE)

} else {
  # === STANDARD CHECKER MODE ===

  # Use Tm (°C) column name for standard mode
  colnames(output_df)[colnames(output_df) == "Tm"] <- "Tm (°C)"

  # re-arrange output df
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

  # Create bar plots using shared functions
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

  df_plot$ASO <- factor(df_plot$ASO, levels = unique(output_df$ASO))

  # Checker uses default ggplot text sizes (no explicit size args) and saves PDFs
  checker_text <- list(axis_text = 11, axis_title = 14, plot_title = 16,
                       legend_text = 11, legend_title = 11, geom_text = 3.5)

  plot_ot_dodged(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
                 "Number of off-targets in TIR regions", "inferno",
                 output_prefix = paste0(path_output, "/plot_ots_start_regions"),
                 text_size = checker_text, save_pdf = TRUE)

  plot_ot_dodged(df_plot[df_plot$off_target_type == "OT in transcriptome", ],
                 "Number of off-targets in whole transcriptome", "viridis",
                 output_prefix = paste0(path_output, "/plot_ots_whole_transcriptome"),
                 text_size = checker_text, save_pdf = TRUE)

  # Screening off-target plots
  checker_single_text <- list(axis_text = 11, axis_title = 14, plot_title = 16, geom_text = 3.5)

  if (!is.na(screen) && screen == "microbiome" && "OT in HMP microbiome" %in% df_plot$off_target_type) {
    df_hmp <- df_plot[df_plot$off_target_type == "OT in HMP microbiome", ]
    df_hmp$ASO <- factor(df_hmp$ASO, levels = unique(output_df$ASO))
    plot_ot_single(df_hmp, "Number of off-targets in HMP microbiome (0 mismatches)", "#440154",
                   output_prefix = paste0(path_output, "/plot_ots_hmp"),
                   text_size = checker_single_text, save_pdf = TRUE)
  }

  if (!is.na(screen) && screen == "human" && "OT in human genome" %in% df_plot$off_target_type) {
    df_human <- df_plot[df_plot$off_target_type == "OT in human genome", ]
    df_human$ASO <- factor(df_human$ASO, levels = unique(output_df$ASO))
    plot_ot_single(df_human, "Number of off-targets in human genome (0 mismatches)", "#31688e",
                   output_prefix = paste0(path_output, "/plot_ots_human"),
                   text_size = checker_single_text, save_pdf = TRUE)
  }
}
