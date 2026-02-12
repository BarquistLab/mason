import sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pna_utils import calculate_purine_stats, calculate_self_complementarity, best_sd_match

res_path = sys.argv[2]
input_pna = sys.argv[1]

# Detect target mode: args 3 and 4 are startreg and mfe FASTA paths
target_mode = len(sys.argv) > 3
if target_mode:
    startreg_fasta = sys.argv[3]
    mfe_fasta = sys.argv[4]

    # Read the 46nt TIR and 60nt MFE sequences
    startreg_record = next(SeqIO.parse(startreg_fasta, "fasta"))
    mfe_record = next(SeqIO.parse(mfe_fasta, "fasta"))
    tir_seq_46 = startreg_record.seq  # 46nt: -30 to +16
    tir_seq_60 = mfe_record.seq       # 60nt: -30 to +30

    # Detect Shine-Dalgarno region (same window as make_pnas.py)
    sd_result = best_sd_match(str(tir_seq_46[18:28]))
    if sd_result:
        sd_motif, sd_rel_start = sd_result
        sd_pos = sd_rel_start + 18  # 0-indexed in 46nt window
        sd_len = len(sd_motif)
    else:
        sd_pos = None
        sd_len = 0

seqs = []
aso_seqs = []
warnings = []
varna_positions = []

if target_mode:
    output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "target_seq", "location",
                                       "SC_bases", "pur_perc", "long_pur_stretch",
                                       "OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm",
                                       "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm"])
    heatmap_list = []
    heatmap_annot = []
    heatmap_labels = []
    tir_mRNA = tir_seq_46.transcribe()
    seqlist = list(str(tir_mRNA))
    seqlist_pnas = list(str(tir_seq_46.complement()))
    count = 0

    for record in SeqIO.parse(input_pna, "fasta"):
        aso_raw = record.seq
        aso_name = record.id
        aso_len = len(aso_raw)

        maxcomp_raw, aso_target = calculate_self_complementarity(aso_raw)
        pur_perc, longest_purine_stretch_raw = calculate_purine_stats(aso_raw)

        # The ASO targets the reverse complement (i.e. the mRNA sequence)
        mRNA_target = str(aso_target)  # This is the mRNA sequence the ASO would bind

        # Search for exact match in the 60nt TIR window
        tir_str = str(tir_seq_60)
        match_pos = tir_str.find(mRNA_target)

        if match_pos == -1:
            warnings.append(f"{aso_name}: sequence does not match the translation initiation region of the target gene")
            continue

        # Position relative to start codon (position 30 in 0-indexed = start codon)
        loc_start = match_pos - 30
        loc_end = loc_start + aso_len
        location = f"{loc_start};{loc_end}"

        target_seq_str = str(tir_seq_60[match_pos:match_pos + aso_len].transcribe())

        added_row = pd.Series([aso_name, str(aso_raw), target_seq_str, location,
                                maxcomp_raw, pur_perc, longest_purine_stretch_raw,
                                None, None, None, None, None, None, None, None],
                              index=output_df.columns)
        output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)

        seqs.append(SeqIO.SeqRecord(aso_target, aso_name, description=""))
        aso_seqs.append(SeqIO.SeqRecord(aso_raw, aso_name, description=""))

        # VARNA positions: 1-indexed for VARNA highlight (within 60nt window)
        varna_positions.append((aso_name, match_pos + 1, match_pos + aso_len))

        # Heatmap data (within 46nt window)
        if match_pos + aso_len <= 46:
            row = [0] * 46
            for j in range(match_pos, match_pos + aso_len):
                row[j] = 1
            # Highlight start codon (positions 30-32 in 0-indexed)
            for j in range(30, 33):
                if j < 46:
                    row[j] = row[j] + 0.4
            # Highlight SD region (only if detected)
            if sd_pos is not None:
                for j in range(sd_pos, sd_pos + sd_len):
                    if j < 46:
                        row[j] = row[j] + 0.2
            heatmap_list.append(row)
            annot_row = [""] * 46
            for j in range(match_pos, match_pos + aso_len):
                annot_row[j] = seqlist_pnas[j]
            heatmap_annot.append(annot_row)
            heatmap_labels.append(aso_name)
        count += 1

    # Write warnings
    with open(res_path + "/outputs/warnings.txt", "w") as wf:
        for w in warnings:
            wf.write(w + "\n")

    # Write VARNA positions
    with open(res_path + "/outputs/varna_positions.tsv", "w") as vf:
        vf.write("aso_name\tstart_1idx\tend_1idx\n")
        for aso_name, start, end in varna_positions:
            vf.write(f"{aso_name}\t{start}\t{end}\n")

    # Write SD position for VARNA highlighting (1-indexed)
    with open(res_path + "/outputs/sd_position.txt", "w") as sf:
        if sd_pos is not None:
            sf.write(f"{sd_pos + 1}\t{sd_pos + sd_len}\n")

    if count > 0:
        SeqIO.write(seqs, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
        SeqIO.write(aso_seqs, res_path + "/reference_sequences/aso_sequences.fasta", "fasta")
        output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)

        # Generate heatmap
        if heatmap_list:
            heatmap_array = np.array(heatmap_list)
            fig, ax = plt.subplots(figsize=(15, max(count / 3, 2)))
            ax = sns.heatmap(heatmap_array, linewidths=.2, linecolor="black",
                             annot=np.array(heatmap_annot), fmt="", cmap="Blues",
                             cbar=False)
            ax.set_xticks(np.arange(0.5, len(seqlist), 1))
            ax.set_xticklabels(seqlist, rotation=0, fontsize=12)
            ax.set_xlabel("mRNA (5' to 3')", fontsize=12)
            ylabels = heatmap_labels
            ax.set_yticklabels(ylabels, rotation=0)
            ax.set_ylabel("ASO sequence", fontsize=12)
            plt.tight_layout()
            plt.savefig(res_path + '/outputs/heatmap.png', dpi=500)
            plt.savefig(res_path + '/outputs/heatmap.svg')
            plt.close()
    else:
        # No ASOs matched — write empty aso_targets.fasta so shell can detect
        open(res_path + "/reference_sequences/aso_targets.fasta", "w").close()

else:
    # Standard checker mode (no target gene)
    output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "SC_bases", "pur_perc", "long_pur_stretch",
                                       "OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm",
                                       "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm"])

    for record in SeqIO.parse(input_pna, "fasta"):
        aso_raw = record.seq

        maxcomp_raw, aso_target = calculate_self_complementarity(aso_raw)
        pur_perc, longest_purine_stretch_raw = calculate_purine_stats(aso_raw)

        aso_name = record.id

        added_row = pd.Series([aso_name, str(aso_raw), maxcomp_raw, pur_perc, longest_purine_stretch_raw,
                                None, None, None, None, None, None, None, None],
                              index=output_df.columns)
        output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)

        seqs.append(SeqIO.SeqRecord(aso_target, aso_name, description=""))
        aso_seqs.append(SeqIO.SeqRecord(aso_raw, aso_name, description=""))

    SeqIO.write(seqs, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
    SeqIO.write(aso_seqs, res_path + "/reference_sequences/shuffled_sequences.fasta", "fasta")
    output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)

print("finished modify_PNAs script!")
