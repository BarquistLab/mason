# r_utils.R — shared utility functions for R scripts in MASON

calculate_pna_mw <- function(seq) {
  # Calculate average molecular weight of an aegPNA oligomer (H-PNA-NH2 form).
  # Monomer residue masses = full monomer - H2O (condensation), derived from:
  #   backbone: N-(2-aminoethyl)glycine + methylenecarbonyl linker
  #   nucleobase: attached at N9 (purines) or N1 (pyrimidines)
  #
  # Args:
  #   seq: character string of nucleotides (A, C, G, T only)
  #
  # Returns:
  #   numeric molecular weight in Da, or NA if invalid nucleotides present

  if (!is.character(seq) || length(seq) != 1 || nchar(seq) == 0) {
    warning("Invalid input: seq must be a non-empty character string")
    return(NA_real_)
  }

  mw_residue <- c(A = 275.27, C = 251.25, G = 291.27, T = 266.26)
  nucleotides <- strsplit(toupper(seq), "")[[1]]


  # Validate nucleotides
  invalid <- nucleotides[!nucleotides %in% names(mw_residue)]
  if (length(invalid) > 0) {
    warning(sprintf("Invalid nucleotide(s) in sequence: %s", paste(unique(invalid), collapse = ", ")))
    return(NA_real_)
  }

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

  # Input validation
  required_cols <- c("ASO", "num_mismatch", "TIR", "probe_seq")
  missing_cols <- setdiff(required_cols, names(all_off_targets))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in all_off_targets: %s", paste(missing_cols, collapse = ", ")))
  }

  unique_asos <- unique(all_off_targets$ASO)
  n_asos <- length(unique_asos)
  mismatch_levels <- 0:3

  # Pre-allocate df_plot: 8 rows per ASO (4 mismatch levels × 2 regions)
  n_rows <- n_asos * 8
  df_plot <- data.frame(
    ASO = character(n_rows),
    off_target_type = character(n_rows),
    transcripts = character(n_rows),
    counts = numeric(n_rows),
    target.sequence = character(n_rows),
    nr_mismatches = numeric(n_rows),
    stringsAsFactors = FALSE
  )

  row_idx <- 1
  for (aso_n in unique_asos) {
    ot_aso <- all_off_targets[all_off_targets$ASO == aso_n, ]
    target_seq <- ot_aso$probe_seq[1]

    # Count off-targets for each mismatch level (vectorized)
    num_tot_ot <- sapply(mismatch_levels, function(mm) sum(ot_aso$num_mismatch == mm))
    num_tir_ot <- sapply(mismatch_levels, function(mm) sum(ot_aso$TIR == "TIR" & ot_aso$num_mismatch == mm))

    # Fill df_plot rows for transcriptome off-targets (4 rows)
    for (mm in mismatch_levels) {
      df_plot[row_idx, ] <- list(aso_n, "OT in transcriptome", "whole transcriptome",
                                  num_tot_ot[mm + 1], target_seq, mm)
      row_idx <- row_idx + 1
    }

    # Fill df_plot rows for TIR off-targets (4 rows)
    for (mm in mismatch_levels) {
      df_plot[row_idx, ] <- list(aso_n, "OT in TIR regions", "start regions",
                                  num_tir_ot[mm + 1], target_seq, mm)
      row_idx <- row_idx + 1
    }

    # Update output_df with OT counts
    ot_cols <- c("OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm",
                 "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm")
    ot_values <- c(num_tir_ot, num_tot_ot)

    if (index_by_name) {
      output_df[aso_n, ot_cols] <- ot_values
    } else {
      output_df[output_df$ASO == aso_n, ot_cols] <- ot_values
    }
  }

  return(list(output_df = output_df, df_plot = df_plot))
}
