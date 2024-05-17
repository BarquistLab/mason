library(ggplot2)
library(dplyr)
library(readr)
library(kableExtra)
library(rmelting)
library(viridis)




print("hello world")

path_output <- commandArgs(trailingOnly = TRUE)[1]
#path_output <- "./browser/pnag/static/data/2024_05_05_09_49_00/outputs"


# Load data
ot_table <- paste0(path_output, "/offtargets_fulltranscripts_sorted.tab")

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
  output_df[output_df$ASO == aso_n, "OT_TIR_0mm"] <- num_tir_ot_0mm
    output_df[output_df$ASO == aso_n, "OT_TIR_1mm"] <- num_tir_ot_1mm
    output_df[output_df$ASO == aso_n, "OT_TIR_2mm"] <- num_tir_ot_2mm
    output_df[output_df$ASO == aso_n, "OT_TIR_3mm"] <- num_tir_ot_3mm

  output_df[output_df$ASO == aso_n, "OT_tot_0mm"] <- num_tot_ot_0mm
    output_df[output_df$ASO == aso_n, "OT_tot_1mm"] <- num_tot_ot_1mm
    output_df[output_df$ASO == aso_n, "OT_tot_2mm"] <- num_tot_ot_2mm
    output_df[output_df$ASO == aso_n, "OT_tot_3mm"] <- num_tot_ot_3mm
}

print(output_df)
print(df_plot)
# remove location column if no target gene was given


# calculate melting temperature for each ASO
output_df[["Tm (°C)"]] <- sapply(output_df$ASO_seq, function(x) {
  melting(sequence = x, nucleic.acid.conc = 0.000008, hybridisation.type = "rnarna", Na.conc = 0.1)$Results$`Melting temperature (C)`
})

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
output_df <- output_df %>% select(ASO, ASO_seq, SC_bases, `%_SC_bases`, `Tm (°C)`, pur_perc, long_pur_stretch,
                                  Mw,
                                  OT_TIR_0mm, OT_TIR_1mm, OT_TIR_2mm, OT_TIR_3mm, OT_tot_0mm, OT_tot_1mm,
                                  OT_tot_2mm, OT_tot_3mm)


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


# save table as html
table_out %>% save_kable(paste0(path_output, "/result_table.html"))



# Save dataframe as CSV
write.csv(output_df, file = paste0(path_output, "/result_table.csv"), row.names = FALSE)

write_xlsx(output_df, paste0(path_output, "/result_table.xlsx"))


# Create bar plot using ggplot2
df_plot$counts <- as.numeric(df_plot$counts)

# make df_plot$off_target_type a factor and order it
df_plot$off_target_type <- factor(df_plot$off_target_type, levels = c("OT in transcriptome", "OT in TIR regions"))

# plot for only ots in TIR regions
df_plot$ASO <- factor(df_plot$ASO, levels = unique(output_df$ASO))

p <- ggplot(df_plot[df_plot$off_target_type == "OT in TIR regions", ],
            aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3)))) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Number of off-targets in TIR regions") +
  labs(x = "ASO sequence", y = "Number of off-targets") +
  theme_classic() +
# define title of legend
    guides(fill = guide_legend(title = "# off-targets")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=13),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        plot.title = element_text(size= 25, hjust = 0.5, face = "bold"),
        # add increased legend fontsize and adjust position (topright of plot inside plot area
        legend.text = element_text(size=12),
        legend.background = element_rect(fill="white", linewidth = .5),
        legend.direction = "horizontal",
        legend.title = element_text(size=12),
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
        legend.text = element_text(size=12),
        legend.background = element_rect(fill=alpha('white', 0.7), linewidth = .5),
        legend.title = element_text(size=12),
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



# Remove temporary files
#file.remove(paste0(path_output, "/result_table.tsv", sep = ""))
