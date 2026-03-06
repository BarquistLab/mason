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

calculate_self_binding_mfe <- function(sequences) {
  # Calculate self-binding MFE for ASO sequences using RNAfold.
  # Converts DNA (T) to RNA (U) before folding.
  #
  # Args:
  #   sequences: character vector of ASO sequences (DNA, with T's)
  #
  # Returns:
  #   numeric vector of MFE values (kcal/mol)

  rna_seqs <- gsub("T", "U", toupper(sequences))

  # Write sequences to temp FASTA
  tmp_fasta <- tempfile(fileext = ".fa")
  writeLines(paste0(">seq", seq_along(rna_seqs), "\n", rna_seqs), tmp_fasta)

  # Run RNAfold
  output <- system2("RNAfold", args = c("--noGU", "-T", "37", "--noPS"),
                     stdin = tmp_fasta, stdout = TRUE, stderr = FALSE)
  unlink(tmp_fasta)

  # Parse MFE from lines containing parenthesized energy values like "...((-2.30))"
  structure_lines <- grep("\\(\\s*-?[0-9.]+\\)", output, value = TRUE)
  mfe_values <- as.numeric(gsub(".*\\(\\s*(-?[0-9.]+)\\).*", "\\1", structure_lines))

  return(mfe_values)
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
    # Subtract 1 from 0mm counts to exclude the self-match (ASO binding its own target)
    num_tot_ot <- sapply(mismatch_levels, function(mm) sum(ot_aso$num_mismatch == mm))
    num_tot_ot[1] <- max(0, num_tot_ot[1] - 1)
    num_tir_ot <- sapply(mismatch_levels, function(mm) sum(ot_aso$TIR == "TIR" & ot_aso$num_mismatch == mm))
    num_tir_ot[1] <- max(0, num_tir_ot[1] - 1)

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

plot_ot_dodged <- function(df_plot_subset, title, viridis_option, legend_title = "# mismatches",
                           output_prefix, wplot = NULL, text_size = NULL, save_pdf = FALSE) {
  # Create a dodged bar plot for off-targets with mismatch-level fill.
  # Args:
  #   df_plot_subset: data frame filtered to one off_target_type
  #   title: plot title
  #   viridis_option: viridis palette option ("inferno", "viridis", etc.)
  #   legend_title: legend title string
  #   output_prefix: file path prefix (without extension)
  #   wplot: plot width (NULL for default ggsave width)
  #   text_size: list with axis_text, axis_title, plot_title, legend_text, legend_title, geom_text sizes
  #              (NULL uses defaults matching mason style)
  #   save_pdf: if TRUE, also save as PDF

  if (is.null(text_size)) {
    text_size <- list(axis_text = 13, axis_title = 20, plot_title = 25,
                      legend_text = 15, legend_title = 15, geom_text = 5)
  }

  p <- ggplot(df_plot_subset,
              aes(x = ASO, y = counts, fill = factor(nr_mismatches, levels = c(0, 1, 2, 3)))) +
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle(title) +
    labs(x = "ASO sequence", y = "Number of off-targets") +
    theme_classic() +
    guides(fill = guide_legend(title = legend_title)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = text_size$axis_text),
          axis.text.y = element_text(size = text_size$axis_text + 2),
          axis.title = element_text(size = text_size$axis_title),
          plot.title = element_text(size = text_size$plot_title, hjust = 0.5, face = "bold"),
          legend.text = element_text(size = text_size$legend_text),
          legend.background = element_rect(fill = "white", linewidth = .5),
          legend.direction = "horizontal",
          legend.title = element_text(size = text_size$legend_title),
          legend.position = "top",
          legend.margin = margin(6, 10, 6, 6)) +
    scale_fill_viridis(discrete = TRUE, direction = -1, option = viridis_option) +
    geom_text(aes(label = counts), position = position_dodge(width = 0.9),
              vjust = -0.25, size = text_size$geom_text)

  width_arg <- if (!is.null(wplot)) list(width = wplot) else list()
  do.call(ggsave, c(list(paste0(output_prefix, ".png"), p, limitsize = FALSE), width_arg))
  do.call(ggsave, c(list(paste0(output_prefix, ".svg"), p, limitsize = FALSE), width_arg))
  if (save_pdf) {
    do.call(ggsave, c(list(paste0(output_prefix, ".pdf"), p, limitsize = FALSE), width_arg))
  }

  invisible(p)
}

plot_ot_single <- function(df_subset, title, fill_color, output_prefix,
                           wplot = NULL, text_size = NULL, save_pdf = FALSE) {
  # Create a single-color bar plot for HMP/human off-targets.
  # Args:
  #   df_subset: data frame filtered to one screening type
  #   title: plot title
  #   fill_color: bar fill color
  #   output_prefix: file path prefix (without extension)
  #   wplot: plot width (NULL for default ggsave width)
  #   text_size: list with axis_text, axis_title, plot_title, geom_text sizes
  #   save_pdf: if TRUE, also save as PDF

  if (is.null(text_size)) {
    text_size <- list(axis_text = 13, axis_title = 20, plot_title = 25, geom_text = 7)
  }

  df_subset$counts <- as.numeric(df_subset$counts)

  p <- ggplot(df_subset, aes(x = ASO, y = counts)) +
    geom_bar(stat = "identity", fill = fill_color) +
    ggtitle(title) +
    labs(x = "ASO sequence", y = "Number of off-targets") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = text_size$axis_text),
          axis.text.y = element_text(size = text_size$axis_text + 2),
          axis.title = element_text(size = text_size$axis_title),
          plot.title = element_text(size = text_size$plot_title, hjust = 0.5, face = "bold")) +
    geom_text(aes(label = counts), vjust = -0.25, size = text_size$geom_text)

  width_arg <- if (!is.null(wplot)) list(width = wplot) else list()
  do.call(ggsave, c(list(paste0(output_prefix, ".png"), p, limitsize = FALSE), width_arg))
  do.call(ggsave, c(list(paste0(output_prefix, ".svg"), p, limitsize = FALSE), width_arg))
  if (save_pdf) {
    do.call(ggsave, c(list(paste0(output_prefix, ".pdf"), p, limitsize = FALSE), width_arg))
  }

  invisible(p)
}

process_essential_gene_offtargets <- function(path_output, all_off_targets, output_df, df_plot,
                                              index_by_name = TRUE) {
  # Filter TIR off-targets to user-specified essential genes, count by mismatch level,
  # update output_df and df_plot, and export filtered tables.
  #
  # Args:
  #   path_output: output directory path (must contain essential_genes.txt)
  #   all_off_targets: off-target data frame with ASO, num_mismatch, TIR, locus_tag columns
  #   output_df: result table data frame
  #   df_plot: plotting data frame
  #   index_by_name: if TRUE, index output_df by row name (mason); if FALSE, filter by ASO column
  #
  # Returns:
  #   list(output_df, df_plot)

  essential_file <- file.path(path_output, "essential_genes.txt")
  if (!file.exists(essential_file)) {
    return(list(output_df = output_df, df_plot = df_plot))
  }

  essential_genes <- readLines(essential_file)
  essential_genes <- trimws(essential_genes)
  essential_genes <- essential_genes[nchar(essential_genes) > 0]

  if (length(essential_genes) == 0) {
    return(list(output_df = output_df, df_plot = df_plot))
  }

  # Filter to TIR off-targets in essential genes
  ess_ot <- all_off_targets[all_off_targets$TIR == "TIR" & all_off_targets$locus_tag %in% essential_genes, ]

  unique_asos <- unique(all_off_targets$ASO)
  mismatch_levels <- 0:3
  ot_cols <- c("OT_ess_TIR_0mm", "OT_ess_TIR_1mm", "OT_ess_TIR_2mm", "OT_ess_TIR_3mm")

  # Initialize columns
  for (col in ot_cols) {
    output_df[[col]] <- 0
  }

  for (aso_n in unique_asos) {
    ot_aso <- ess_ot[ess_ot$ASO == aso_n, ]
    ot_aso_all <- all_off_targets[all_off_targets$ASO == aso_n, ]
    target_seq <- ot_aso_all$probe_seq[1]

    counts <- sapply(mismatch_levels, function(mm) sum(ot_aso$num_mismatch == mm))

    # Update output_df
    if (index_by_name) {
      output_df[aso_n, ot_cols] <- counts
    } else {
      output_df[output_df$ASO == aso_n, ot_cols] <- counts
    }

    # Append to df_plot
    for (mm in mismatch_levels) {
      df_plot <- rbind(df_plot, data.frame(
        ASO = aso_n,
        off_target_type = "OT in TIR of essential genes",
        transcripts = "essential genes",
        counts = counts[mm + 1],
        target.sequence = as.character(target_seq),
        nr_mismatches = mm,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Export filtered essential gene off-targets
  if (nrow(ess_ot) > 0) {
    write_csv(ess_ot, file = paste0(path_output, "/offtargets_essential_genes_sorted.csv"))
    write_xlsx(ess_ot, paste0(path_output, "/offtargets_essential_genes_sorted.xlsx"))
  }

  return(list(output_df = output_df, df_plot = df_plot))
}

process_screen_offtargets <- function(screen_type, path_output, output_df, df_plot,
                                      index_by_name = TRUE) {
  # Process HMP microbiome or human transcriptome screening off-targets.
  # Reads the screening file, counts 0mm off-targets per ASO, appends to df_plot,
  # exports CSV/Excel, and updates output_df.
  #
  # Args:
  #   screen_type: "microbiome" or "human"
  #   path_output: output directory path
  #   output_df: result table data frame
  #   df_plot: plotting data frame
  #   index_by_name: if TRUE, index output_df by row name (mason); if FALSE, filter by ASO column
  #
  # Returns:
  #   list(output_df, df_plot)

  if (screen_type == "microbiome") {
    screen_file <- paste0(path_output, "/offtargets_microbiome_sorted.tab")
    ot_col <- "OT_HMP_0mm"
    ot_label <- "OT in HMP microbiome"
    ot_transcripts <- "HMP microbiome"
    csv_name <- "/offtargets_hmp_sorted.csv"
    xlsx_name <- "/offtargets_hmp_sorted.xlsx"
  } else if (screen_type == "human") {
    screen_file <- paste0(path_output, "/offtargets_human_sorted.tab")
    ot_col <- "OT_GRCh38_0mm"
    ot_label <- "OT in human transcriptome"
    ot_transcripts <- "human transcriptome"
    csv_name <- "/offtargets_human_sorted.csv"
    xlsx_name <- "/offtargets_human_sorted.xlsx"
  } else {
    return(list(output_df = output_df, df_plot = df_plot))
  }

  if (!file.exists(screen_file)) {
    return(list(output_df = output_df, df_plot = df_plot))
  }

  screen_ot <- read_tsv(screen_file, col_names = TRUE)
  unique_asos <- unique(screen_ot$probe_id)
  output_df[[ot_col]] <- 0

  for (aso_n in unique_asos) {
    ot_aso <- screen_ot[screen_ot$probe_id == aso_n, ]
    n_0mm <- sum(ot_aso$num_mismatch == 0)

    if (index_by_name) {
      output_df[aso_n, ot_col] <- n_0mm
    } else {
      output_df[output_df$ASO == aso_n, ot_col] <- n_0mm
    }

    df_plot <- rbind(df_plot, data.frame(
      ASO = aso_n,
      off_target_type = ot_label,
      transcripts = ot_transcripts,
      counts = n_0mm,
      target.sequence = ot_aso$probe_seq[1],
      nr_mismatches = 0,
      stringsAsFactors = FALSE
    ))
  }

  # Export cleaned screening table
  screen_ot_clean <- screen_ot %>%
    mutate(ASO = probe_id) %>%
    select(-probe_id)
  write_csv(screen_ot_clean, file = paste0(path_output, csv_name))
  write_xlsx(screen_ot_clean, paste0(path_output, xlsx_name))

  return(list(output_df = output_df, df_plot = df_plot))
}
