import sys
from Bio import SeqIO
from difflib import SequenceMatcher
from Bio.SeqIO import SeqRecord
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import pandas as pd

# import used gene length:
length = int(sys.argv[1])
res_path = sys.argv[2]

# import fasta file:
target_regions = list(SeqIO.parse(res_path + "/reference_sequences" + "/targetgene_startreg.fasta", "fasta"))

list_aso_sequences = []
list_aso_targets = []
count = 1

# go through target regions and
for r in range(len(target_regions)):
    gene = re.sub("::.*", "", target_regions[r].id)
    seq = target_regions[r].seq
    print(seq)
    seqlist = [i for i in seq.transcribe().__str__()]  # add sequence list for heatmap
    seqlist_pnas = [i for i in seq.complement().__str__()]
    heatmap_list = []
    heatmap_annot = []
    output_df = pd.DataFrame(columns=["ASO_name", "ASO_sequence", "mRNA_target_seq", "position_from_start",
                                      "self_comp_bases", "melting_temperature",
                                      "OT_transcriptome", "OT_TIR_regions"])
    # identify location of SD regions:
    if seq[5:28].find("GAGG") != -1:
        sd_pos = seq[5:28].find("GAGG") + 5

    for i in range(5, len(seq)-length):
        s = seq[i:i+length]
        s_rev = s[::-1]
        aso = s.reverse_complement()

        maxcomp = max([SequenceMatcher(None, aso, s).find_longest_match(0, len(aso), 0, len(s)).size,
                       SequenceMatcher(None, aso, s_rev).find_longest_match(0, len(aso), 0, len(s_rev)).size])
        # design PNA for regions overlapping SC or SD region:
        if ((i in range(33-length, 31)) or (i in range(5, 33-length) and "GAGG" in s)) and maxcomp < length*0.5:
            aso_name = gene+"_ASO_"+str(count).zfill(3)
            # add to dataframe:
            added_row = pd.Series([aso_name, aso.__str__(), s.transcribe().__str__(), str(i-30) + ";" + str(i-30+length),
                                   maxcomp, None, None, None], index=output_df.columns)
            output_df = output_df.append(added_row, ignore_index=True)
            # save ASO sequences:
            list_aso_sequences += [SeqRecord(aso, id=aso_name, description="")]
            # for fasta with reverse as well to get mismatches with seqmap:
            list_aso_targets += [SeqRecord(s, id=aso_name, description="")]
            list_aso_targets += [SeqRecord(s_rev, id=aso_name + "_rev", description="")]
            heatmap_list += [i*[0] + length*[1] + (len(seq)-i-length)*[0]]
            heatmap_list[count-1][30:33] = [0.4 + i for i in heatmap_list[count-1][30:33]]
            if seq[5:28].find("GAGG") != -1:
                heatmap_list[count - 1][sd_pos:sd_pos+4] = [0.2 + i for i in heatmap_list[count - 1][sd_pos:sd_pos+4] ]
            heatmap_annot += [i*[""] + seqlist_pnas[i:i+length] + (len(seq)-i-length)*[""]]
            count += 1
    heatmap_array = np.array(heatmap_list)

    SeqIO.write(list_aso_sequences, res_path + "/reference_sequences/aso_sequences.fasta", "fasta")
    SeqIO.write(list_aso_targets, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
    # safe output df:
    print(output_df)
    output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)
    # run rscript:
    subprocess.run(["Rscript", "./pnag/melting.R", res_path + "/outputs/result_table.tsv"])

    # seaborn:
    fig, ax = plt.subplots(figsize=(15, count/3))
    ax = sns.heatmap(heatmap_array, linewidths=.2, linecolor="black",
                     annot=np.array(heatmap_annot), fmt="", cmap="Blues",
                     cbar=False)
    ax.set_xticks(np.arange(0.5, len(seqlist), 1))
    ax.set_xticklabels(seqlist, rotation=0, fontsize=12)
    ax.set_xlabel("mRNA (5' to 3')", fontsize=12)
    ylabels = ["ASO_" + str(i + 1).zfill(3) for i in range(count-1)]
    ax.set_yticklabels(ylabels, rotation=0)
    ax.set_ylabel("ASO sequence", fontsize=12)
    plt.tight_layout()
    plt.savefig(res_path + '/outputs/heatmap.png', dpi=500)

    count = 0




