import sys
from Bio import SeqIO
from cdifflib import CSequenceMatcher
from Bio.SeqIO import SeqRecord
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import pandas as pd
from pna_utils import calculate_purine_stats, calculate_self_complementarity

# import used gene length:
length = int(sys.argv[1])
res_path = sys.argv[2]
bases_before = int(sys.argv[3])

print("BASES BEFORE:")
print(bases_before)

# import fasta file:
target_regions = list(SeqIO.parse(res_path + "/reference_sequences" + "/targetgene_startreg.fasta", "fasta"))

list_aso_sequences = []
list_aso_targets = []
count = 1
sd = r"AGGAGG|GGAGG|GGAG|AAGGA|AGGA|GAGG|AGG|GGA"  # Shine Dalgarno sequence

# Longest-to-shortest order is still good style, but not sufficient by itself
SD_PATTERN = r"(AGGAGG|GGAGG|AAGGA|AGGA|GAGG|GGAG|AGG|GGA)"

def best_sd_match(seq_slice: str):
    # Overlapping matches via lookahead; capture the motif in group(1)
    pat = r"(?=(" + SD_PATTERN + r"))"
    matches = [(m.group(1), m.start()) for m in re.finditer(pat, seq_slice)]
    if not matches:
        return None
    # Pick the longest; if tied, pick the leftmost
    motif, start = max(matches, key=lambda t: (len(t[0]), -t[1]))
    return motif, start

# go through target regions and
for r in range(len(target_regions)):
    gene = re.sub("::.*", "", target_regions[r].id)
    seq = target_regions[r].seq
    seqlist = [i for i in seq.transcribe().__str__()]  # add sequence list for heatmap
    seqlist_pnas = [i for i in seq.complement().__str__()]
    heatmap_list = []
    heatmap_annot = []
    output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "target_seq", "location",
                                      "SC_bases", "pur_perc", "long_pur_stretch", "Tm",
                                      "OT_tot", "OT_TIR", "OT_GRCh38", "OT_HMP", "OT_cust"])
    # identify location of SD regions:

    if best_sd_match(str(seq[18:28])):
        sd_pos = best_sd_match(str(seq[18:28]))[1] + 18
        print("SD found in region:", str(seq[18:28]))
        print("SD position at:", sd_pos)
        print("SD sequence:", best_sd_match(str(seq[18:28]))[0])
        sd = best_sd_match(str(seq[18:28]))[0]

        if bases_before != 0:
            print("designing PNA in region before SD, as specified by user. bases:" , bases_before)
        else:
            for i in range(sd_pos - (length - len(sd)), sd_pos + 1):
                print("designing PNA overlapping SD region at pos:", i)
                s = seq[i:i + length]
                print(s, "=========================")
                aso = s.reverse_complement()

                maxcomp, _ = calculate_self_complementarity(aso)

                if maxcomp < length+1:
                    aso_name = gene + "_ASO_" + str(count).zfill(3)
                    pur_perc, longest_purine_stretch = calculate_purine_stats(aso)

                    # add to dataframe:
                    added_row = pd.Series([aso_name, aso.__str__(), s.transcribe().__str__(),
                                           str(i - 30) + ";" + str(i - 30 + length),
                                           maxcomp, pur_perc, longest_purine_stretch, None, None, None, None, None,
                                           None],
                                          index=output_df.columns)
                    # output_df = output_df.append(added_row, ignore_index=True)
                    # change above to concat and not append:
                    output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)
                    # save ASO sequences:
                    list_aso_sequences += [SeqRecord(aso, id=aso_name, description="")]
                    # for fasta with reverse as well to get mismatches with seqmap:
                    list_aso_targets += [SeqRecord(s, id=aso_name, description="")]

                    heatmap_list += [i * [0] + length * [1] + (len(seq) - i - length) * [0]]
                    heatmap_list[count - 1][30:33] = [0.4 + i for i in heatmap_list[count - 1][30:33]]
                    heatmap_list[count - 1][sd_pos:sd_pos + len(sd)] = [0.2 + i for i in
                                                                        heatmap_list[count - 1][sd_pos:sd_pos + len(sd)]]
                    heatmap_annot += [i * [""] + seqlist_pnas[i:i + length] + (len(seq) - i - length) * [""]]
                    count += 1
    else:
        print("No SD found in region:", str(seq[18:28]))





    for i in range(30-(length-3), 31):   # 30 is start of cds
        s = seq[i:i+length]
        aso = s.reverse_complement()

        maxcomp, _ = calculate_self_complementarity(aso)

        # design PNA for regions overlapping SC or SD region:
        if maxcomp < length+1:
            aso_name = gene+"_ASO_"+str(count).zfill(3)
            pur_perc, longest_purine_stretch = calculate_purine_stats(aso)

            # add to dataframe:
            added_row = pd.Series([aso_name, aso.__str__(), s.transcribe().__str__(), str(i-30) + ";" + str(i-30+length),
                                   maxcomp, pur_perc, longest_purine_stretch, None, None, None, None, None, None],
                                  index=output_df.columns)
            #output_df = output_df.append(added_row, ignore_index=True)
            # change above to concat and not append:
            output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)
            # save ASO sequences:
            list_aso_sequences += [SeqRecord(aso, id=aso_name, description="")]
            # for fasta with reverse as well to get mismatches with seqmap:
            list_aso_targets += [SeqRecord(s, id=aso_name, description="")]

            heatmap_list += [i*[0] + length*[1] + (len(seq)-i-length)*[0]]
            heatmap_list[count-1][30:33] = [0.4 + i for i in heatmap_list[count-1][30:33]]
            if seq[15:28].find(sd) != -1:
                heatmap_list[count - 1][sd_pos:sd_pos+len(sd)] = [0.2 + i for i in heatmap_list[count - 1][sd_pos:sd_pos+len(sd)] ]
            heatmap_annot += [i*[""] + seqlist_pnas[i:i+length] + (len(seq)-i-length)*[""]]
            count += 1
    heatmap_array = np.array(heatmap_list)

    SeqIO.write(list_aso_sequences, res_path + "/reference_sequences/aso_sequences.fasta", "fasta")
    SeqIO.write(list_aso_targets, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
    # safe output df:
    output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)
    # run rscript:
    print(res_path + "/outputs/result_table.tsv")
    subprocess.run(["Rscript", "./pnag/melting.R", res_path + "/outputs/result_table.tsv"])

    # seaborn:
    print(np.array(heatmap_annot))
    print(heatmap_array)

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
    plt.savefig(res_path + '/outputs/heatmap.svg')

    count = 0
