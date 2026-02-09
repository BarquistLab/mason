# r_utils.R — shared utility functions for R scripts in MASON

calculate_pna_mw <- function(seq) {
  # Calculate average molecular weight of an aegPNA oligomer (H-PNA-NH2 form).
  # Monomer residue masses = full monomer - H2O (condensation), derived from:
  #   backbone: N-(2-aminoethyl)glycine + methylenecarbonyl linker
  #   nucleobase: attached at N9 (purines) or N1 (pyrimidines)
  mw_residue <- c(A=275.27, C=251.25, G=291.27, T=266.26)
  nucleotides <- strsplit(seq, "")[[1]]
  # Terminal groups: H (N-terminus) + NH2 (C-terminus) = 17.03
  mw <- sum(mw_residue[nucleotides]) + 17.03
  return(mw)
}

count_offtargets_by_mismatch <- function(all_off_targets, output_df, index_by_name = TRUE) {
  # Count off-targets by mismatch level (0-3mm) for both total and TIR regions.
  # Builds df_plot and updates output_df with OT counts.
  #
  # Args:
  #   all_off_targets: filtered off-target data frame with ASO, num_mismatch, TIR, probe_seq columns
  #   output_df: result table data frame to update with OT counts
  #   index_by_name: if TRUE, index output_df by row name (mason); if FALSE, filter by ASO column (checker)
  #
  # Returns:
  #   list(output_df, df_plot)

  df_plot <- data.frame(ASO = character(0), "off_target_type" = character(0),
                        "transcripts" = character(0), counts = numeric(0),
                        "target_sequence" = character(0), "nr_mismatches" = numeric(0))

  for (i in unique(all_off_targets$ASO)) {
    target_seq <- all_off_targets[all_off_targets$ASO == i, ]$probe_seq[1]
    aso_n <- i
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
                                          "target sequence" = target_seq, "nr_mismatches" = 0))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                          "transcripts" = "whole transcriptome", counts = num_tot_ot_1mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 1))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                          "transcripts" = "whole transcriptome", counts = num_tot_ot_2mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 2))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in transcriptome",
                                          "transcripts" = "whole transcriptome", counts = num_tot_ot_3mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 3))

    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                          "transcripts" = "start regions", counts = num_tir_ot_0mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 0))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                          "transcripts" = "start regions", counts = num_tir_ot_1mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 1))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                          "transcripts" = "start regions", counts = num_tir_ot_2mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 2))
    df_plot <- rbind(df_plot, data.frame(ASO = aso_n, "off_target_type" = "OT in TIR regions",
                                          "transcripts" = "start regions", counts = num_tir_ot_3mm,
                                          "target sequence" = target_seq, "nr_mismatches" = 3))

    # assign values to output_df
    if (index_by_name) {
      output_df[aso_n, "OT_TIR_0mm"] <- num_tir_ot_0mm
      output_df[aso_n, "OT_TIR_1mm"] <- num_tir_ot_1mm
      output_df[aso_n, "OT_TIR_2mm"] <- num_tir_ot_2mm
      output_df[aso_n, "OT_TIR_3mm"] <- num_tir_ot_3mm
      output_df[aso_n, "OT_tot_0mm"] <- num_tot_ot_0mm
      output_df[aso_n, "OT_tot_1mm"] <- num_tot_ot_1mm
      output_df[aso_n, "OT_tot_2mm"] <- num_tot_ot_2mm
      output_df[aso_n, "OT_tot_3mm"] <- num_tot_ot_3mm
    } else {
      output_df[output_df$ASO == aso_n, "OT_TIR_0mm"] <- num_tir_ot_0mm
      output_df[output_df$ASO == aso_n, "OT_TIR_1mm"] <- num_tir_ot_1mm
      output_df[output_df$ASO == aso_n, "OT_TIR_2mm"] <- num_tir_ot_2mm
      output_df[output_df$ASO == aso_n, "OT_TIR_3mm"] <- num_tir_ot_3mm
      output_df[output_df$ASO == aso_n, "OT_tot_0mm"] <- num_tot_ot_0mm
      output_df[output_df$ASO == aso_n, "OT_tot_1mm"] <- num_tot_ot_1mm
      output_df[output_df$ASO == aso_n, "OT_tot_2mm"] <- num_tot_ot_2mm
      output_df[output_df$ASO == aso_n, "OT_tot_3mm"] <- num_tot_ot_3mm
    }
  }

  return(list(output_df = output_df, df_plot = df_plot))
}
