
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
target_gene <- commandArgs(trailingOnly = TRUE)[2]
screen <- commandArgs(trailingOnly = TRUE)[3]


# path_output <- "./browser/pnag/static/data/2024_11_18_16_12_24/b0185/outputs"


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
library(writexl)
write_xlsx(all_off_targets, paste0(path_output, "/offtargets_fulltranscripts_sorted.xlsx"))



output_df <- read.table(paste0(path_output, "/result_table.tsv"), sep = "\t", header = TRUE)


# Count off-targets using shared function (index_by_name=TRUE for mason)
ot_results <- count_offtargets_by_mismatch(all_off_targets, output_df, index_by_name = TRUE)
output_df <- ot_results$output_df
df_plot <- ot_results$df_plot

## Microbiome off-targets (0mm exact matches against HMP database)
if (!is.na(screen) && screen == "microbiome") {
  hmp_file <- paste0(path_output, "/offtargets_microbiome_sorted.tab")
  if (file.exists(hmp_file)) {
    hmp_ot <- read_table(hmp_file, col_names = TRUE)

    # Count 0mm off-targets per ASO and add to output_df
    unique_asos <- unique(hmp_ot$probe_id)
    output_df[["OT_HMP_0mm"]] <- 0
    for (aso_n in unique_asos) {
      ot_aso <- hmp_ot[hmp_ot$probe_id == aso_n, ]
      n_0mm <- sum(ot_aso$num_mismatch == 0)
      output_df[aso_n, "OT_HMP_0mm"] <- n_0mm

      # Append to df_plot for microbiome
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

    # Export cleaned microbiome table
    hmp_ot_clean <- hmp_ot %>%
      mutate(ASO = probe_id) %>%
      select(-probe_id)
    write_csv(hmp_ot_clean, file = paste0(path_output, "/offtargets_hmp_sorted.csv"))
    write_xlsx(hmp_ot_clean, paste0(path_output, "/offtargets_hmp_sorted.xlsx"))
  }
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




# Remove temporary files
#file.remove(paste0(path_output, "/result_table.tsv", sep = ""))
