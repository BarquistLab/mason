
Sys.setenv(HOME = "/home/jakob")


library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(rmelting)
library(viridis)
library(writexl)




print("hello world")

path_output <- commandArgs(trailingOnly = TRUE)[1]
target_gene <- commandArgs(trailingOnly = TRUE)[2]


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



# columns of output_df are : ASO	ASO_seq	SC_bases	pur_perc	long_pur_stretch	OT_TIR_0mm	OT_TIR_1mm	OT_TIR_2mm	OT_TIR_3mm
# OT_tot_0mm	OT_tot_1mm	OT_tot_2mm	OT_tot_3mm

# Add columns for df_plot
df_plot <- data.frame(ASO = character(0), "off_target_type" = character(0),
                      "transcripts" = character(0), counts = numeric(0),
                      "target_sequence" = character(0), "nr_mismatches" = numeric(0))

for (i in unique(all_off_targets$ASO)) {
  target_seq <- all_off_targets[all_off_targets$ASO == i, ]$probe_seq[1]
  aso_n <-  i
  ot_aso <- all_off_targets[all_off_targets$ASO == i, ]

  num_tot_ot_0mm <- nrow(ot_aso[ot_aso$num_mismatch == 0, ])
  num_tot_ot_1mm <- nrow(ot_aso[ot_aso$num_mismatch == 1, ])
  num_tot_ot_2mm <- nrow(ot_aso[ot_aso$num_mismatch == 2, ])
  num_tot_ot_3mm <- nrow(ot_aso[ot_aso$num_mismatch == 3, ])

  num_tir_ot_0mm <- nrow(ot_aso[ot_aso$TIR == "TIR" & ot_aso$num_mismatch == 0, ])
  num_tir_ot_1mm <- nrow(ot_aso[ot_aso$TIR == "TIR" & ot_aso$num_mismatch == 1, ])
  num_tir_ot_2mm <- nrow(ot_aso[ot_aso$TIR == "TIR" & ot_aso$num_mismatch == 2, ])
  num_tir_ot_3mm <- nrow(ot_aso[ot_aso$TIR == "TIR" & ot_aso$num_mismatch == 3, ])

  df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                        "transcripts" = "whole transcriptome", counts = num_tot_ot_0mm,
                                        "target sequence" = target_seq, "nr_mismatches" =0))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                        "transcripts" = "whole transcriptome", counts = num_tot_ot_1mm,
                                        "target sequence" = target_seq,"nr_mismatches" = 1))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                        "transcripts" = "whole transcriptome", counts = num_tot_ot_2mm,
                                        "target sequence" = target_seq, "nr_mismatches" =2))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                        "transcripts" = "whole transcriptome", counts = num_tot_ot_3mm,
                                        "target sequence" = target_seq, "nr_mismatches" =3))



  df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                        "transcripts" = "start regions", counts = num_tir_ot_0mm,
                                        "target sequence" = target_seq,"nr_mismatches" = 0))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                        "transcripts" = "start regions", counts = num_tir_ot_1mm,
                                        "target sequence" = target_seq, "nr_mismatches" =1))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                        "transcripts" = "start regions", counts = num_tir_ot_2mm,
                                        "target sequence" = target_seq, "nr_mismatches" =2))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                        "transcripts" = "start regions", counts = num_tir_ot_3mm,
                                        "target sequence" = target_seq,"nr_mismatches" = 3))


  # assign values to output_df
  print(aso_n)
  print(i)
  output_df[aso_n, "OT_TIR_0mm"] <- num_tir_ot_0mm
    output_df[aso_n, "OT_TIR_1mm"] <- num_tir_ot_1mm
    output_df[aso_n, "OT_TIR_2mm"] <- num_tir_ot_2mm
    output_df[aso_n, "OT_TIR_3mm"] <- num_tir_ot_3mm

  output_df[aso_n, "OT_tot_0mm"] <- num_tot_ot_0mm
    output_df[aso_n, "OT_tot_1mm"] <- num_tot_ot_1mm
    output_df[aso_n, "OT_tot_2mm"] <- num_tot_ot_2mm
    output_df[aso_n, "OT_tot_3mm"] <- num_tot_ot_3mm
}

print(output_df)
print(df_plot)
# remove location column if no target gene was given


# calculate Mw for each ASO
calculate_pna_mw <- function(seq) {
  # Define the molecular weight of each nucleotide (nucleobase only)
  mw_nucleotide <- c(A=135.13, T=125.06, G=152.12, C=111.07)
  # Define the molecular weight of the peptide bond
  mw_peptidebond <- 9*1.01 + 2*14.01 + 2*16.00 + 4*12.01
  # Split the DNA sequence into individual nucleotides
  nucleotides <- strsplit(seq, "")[[1]]
  # Calculate the molecular weight of the sequence, accounting for the peptide bond
  mw <- sum(mw_nucleotide[nucleotides]) + mw_peptidebond * (length(nucleotides) - 1) + 1.01*2
  return(mw)
}

output_df[["Mw"]] <- sapply(output_df$ASO_seq, calculate_pna_mw)


output_df[["%_SC_bases"]] <- round(output_df$SC_bases / nchar(output_df$ASO_seq) * 100, 2)

# re-arrange output df by ASO 	ASO_seq	SC_bases  % SC bases  Tm (°C)	pur_perc 	long_pur_stretch 	OT_TIR_0mm 	OT_TIR_1mm
# OT_TIR_2mm 	OT_TIR_3mm 	OT_tot_0mm 	OT_tot_1mm 	OT_tot_2mm 	OT_tot_3mm
output_df["gene"] <- gsub("_ASO_.*", "", rownames(output_df))

# remove rownames
rownames(output_df) <- NULL

# make output_df a tibble
output_df <- as_tibble(output_df)


output_df <- output_df %>% select(ASO, gene, ASO_seq, SC_bases, `%_SC_bases`, `Tm`, pur_perc, long_pur_stretch,
                                  Mw,
                                  OT_TIR_0mm, OT_TIR_1mm, OT_TIR_2mm, OT_TIR_3mm, OT_tot_0mm, OT_tot_1mm,
                                  OT_tot_2mm, OT_tot_3mm)

# get CAI from gene (a file with just this number)
cai <- read.csv(paste0(path_output, "/cai_value.txt"), header = FALSE) %>% as.numeric() %>% unlist()

mfe <- read.csv(paste0(path_output, "/mfe_values.txt"), header = FALSE) [3,1] %>% as.numeric() %>% unlist()

saved_table_ml <- output_df %>%
  mutate(CAI = cai) %>%
    mutate(MFE = mfe) %>%
  # select only Tm, % SC bases, pur_perc, OT_TIR_1mm, CAI, MFE
    select(`Tm`, `%_SC_bases`, CAI, OT_TIR_1mm, MFE, pur_perc)

# rename columns to 'Tm', 'sc_bases', 'CAI', 'upec_tir_off_targets_1mm',
#                      'MFE_UPEC', 'purine_percentage'
colnames(saved_table_ml) <- c("Tm", "sc_bases", "CAI", "upec_tir_off_targets_1mm",
                              "MFE_UPEC", "purine_percentage")
# write as csv file
write.csv(saved_table_ml, file = paste0(path_output, "/saved_table_ml.csv"), row.names = FALSE)

# als write output_df as csv
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)

# aslo write df_plot as csv
write.csv(df_plot, file = paste0(path_output, "/df_plot.csv"), row.names = FALSE)




# Remove temporary files
#file.remove(paste0(path_output, "/result_table.tsv", sep = ""))
