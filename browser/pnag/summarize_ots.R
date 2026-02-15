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
target_gene <- commandArgs(trailingOnly = TRUE)[2]
screen <- commandArgs(trailingOnly = TRUE)[3]



# Load data

ot_table <- paste0(path_output, "/offtargets_fulltranscripts_sorted.tab")
print(ot_table)
all_off_targets <- read_table(ot_table, col_names = TRUE)



print(head(all_off_targets))



# Filter mismatches with >+7 cons. matches. use dplyr::filter
all_off_targets <- all_off_targets %>% filter(longest_stretch > 6) %>%
  mutate(position_from_CDS_start = trans_coord - 31) %>%
  mutate(ASO = probe_id) %>%
  # unselect probe_id column and trans_coord column
    select(-probe_id, -trans_coord) %>%
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


# Count off-targets using shared function (index_by_name=TRUE for mason)
ot_results <- count_offtargets_by_mismatch(all_off_targets, output_df, index_by_name = TRUE)
output_df <- ot_results$output_df
df_plot <- ot_results$df_plot

## Microbiome/human off-targets (0mm exact matches)
if (!is.na(screen) && screen %in% c("microbiome", "human")) {
  screen_results <- process_screen_offtargets(screen, path_output, output_df, df_plot,
                                               index_by_name = TRUE)
  output_df <- screen_results$output_df
  df_plot <- screen_results$df_plot
}

## Essential gene off-targets
if (!is.na(screen) && screen == "essential_genes") {
  ess_results <- process_essential_gene_offtargets(path_output, all_off_targets, output_df, df_plot,
                                                    index_by_name = TRUE)
  output_df <- ess_results$output_df
  df_plot <- ess_results$df_plot
}

print(output_df)
print(df_plot)
# remove location column if no target gene was given


# calculate Mw for each ASO
output_df[["Mw"]] <- sapply(output_df$ASO_seq, calculate_pna_mw)


output_df[["%_SC_bases"]] <- round(output_df$SC_bases / nchar(output_df$ASO_seq) , 2)

# re-arrange output df by ASO 	ASO_seq	SC_bases  % SC bases  Tm (°C)	pur_perc 	long_pur_stretch 	OT_TIR_0mm 	OT_TIR_1mm
# OT_TIR_2mm 	OT_TIR_3mm 	OT_tot_0mm 	OT_tot_1mm 	OT_tot_2mm 	OT_tot_3mm
output_df["gene"] <- gsub("_ASO_.*", "", rownames(output_df))

# remove rownames
rownames(output_df) <- NULL

# make output_df a tibble
output_df <- as_tibble(output_df)


base_cols <- c("ASO", "gene", "ASO_seq", "SC_bases", "%_SC_bases", "Tm", "pur_perc", "long_pur_stretch",
               "Mw", "OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm",
               "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm")
if ("OT_HMP_0mm" %in% names(output_df)) {
  base_cols <- c(base_cols, "OT_HMP_0mm")
}
if ("OT_GRCh38_0mm" %in% names(output_df)) {
  base_cols <- c(base_cols, "OT_GRCh38_0mm")
}
if ("OT_ess_TIR_0mm" %in% names(output_df)) {
  base_cols <- c(base_cols, "OT_ess_TIR_0mm", "OT_ess_TIR_1mm", "OT_ess_TIR_2mm", "OT_ess_TIR_3mm")
}
output_df <- output_df %>% select(all_of(base_cols))

# get CAI from gene (a file with just this number)
cai <- read.csv(paste0(path_output, "/cai_value.txt"), header = FALSE) %>% as.numeric() %>% unlist()

mfe <- read.csv(paste0(path_output, "/mfe_values.txt"), header = FALSE) [1,1] %>% as.numeric() %>% unlist()

saved_table_ml <- output_df %>%
  mutate(CAI = cai) %>%
    mutate(MFE = mfe) %>%
  # select only ASO, Tm, % SC bases, pur_perc, OT_TIR_1mm, CAI, MFE
    select(ASO, `Tm`, `%_SC_bases`, CAI, OT_TIR_1mm, MFE, pur_perc)

# rename columns to 'Tm', 'sc_bases', 'CAI', 'upec_tir_off_targets_1mm',
#                      'MFE_UPEC', 'purine_percentage'
colnames(saved_table_ml) <- c("ASO", "Tm", "sc_bases", "CAI", "upec_tir_off_targets_1mm",
                              "MFE_UPEC", "purine_percentage")
# write as csv file
write.csv(saved_table_ml, file = paste0(path_output, "/saved_table_ml.csv"), row.names = FALSE)

# als write output_df as csv
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)

# aslo write df_plot as csv
write.csv(df_plot, file = paste0(path_output, "/df_plot.csv"), row.names = FALSE)


