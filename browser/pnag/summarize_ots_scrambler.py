import os
import matplotlib.pyplot as plt
import re
import pandas as pd
import seaborn as sns
import numpy as np
import sys

ot_table = sys.argv[1] + "/offtargets_fulltranscripts_sorted.tab"
screen = sys.argv[2] if len(sys.argv) > 2 else "none"
all_off_targets = pd.read_table(ot_table, sep='\t', index_col=False)


all_off_targets["ASO"] = all_off_targets["probe_id"].replace("_rev", "", regex=True)
all_off_targets["target_seq"] = all_off_targets["target_seq"].replace("T", "U", regex=True)
all_off_targets["probe_seq"] = all_off_targets["probe_seq"].replace("T", "U", regex=True)
all_off_targets["binding_sequence"] = all_off_targets["binding_sequence"].replace("T", "U", regex=True)

all_off_targets = all_off_targets.rename(columns={"target_seq": "off_target_seq_mRNA", "probe_seq": "mRNA_target_seq",
                        "binding_sequence": "matching_sequence"})

# only use the mismatches with >+7 cons. mm:
all_off_targets = all_off_targets[all_off_targets["longest_stretch"] >= 7]


# add output df used for other things, e.g. Tm:
output_df = pd.read_csv(sys.argv[1] + "/result_table.tsv", sep="\t", index_col=None)

# remove first row of output_df

# add variable for whether it is in the TIR (-20 to +5 start)
all_off_targets.loc[all_off_targets["trans_coord"].isin(range(-20, 5)), "TIR"] = "TIR"
all_off_targets.loc[~all_off_targets["trans_coord"].isin(range(-20, 5)), "TIR"] = "not in TIR"

all_off_targets.to_csv(sys.argv[1] + "/offtargets_fulltranscripts_sorted.csv")
all_off_targets.to_excel(sys.argv[1] + "/offtargets_fulltranscripts_sorted.xlsx")

all_off_targets[all_off_targets["TIR"] == "TIR"].to_csv(sys.argv[1] + "/offtargets_startregions_sorted.csv")

df_plot = pd.DataFrame(columns=["ASO", "off-target type", "transcripts", "counts", "target sequence", "num_mismatch"])


# create dataframe:
for i in all_off_targets["ASO"].unique():
    target_seq = all_off_targets[all_off_targets["ASO"] == i].iloc[0, ]["mRNA_target_seq"]
    aso_n = re.sub(".*_(ASO.*)", "\\1", i)
    ot_aso = all_off_targets[all_off_targets["ASO"] == i]
    # get number of OTs in whole transcriptome:
    # 3mm:
    num_tot_ot_3mm = ot_aso[ot_aso["num_mismatch"] == 3].shape[0]
    num_tir_ot_3mm = ot_aso[(ot_aso["TIR"] == "TIR") & (ot_aso["num_mismatch"] == 3)].shape[0]
    # 2mm:
    num_tot_ot_2mm = ot_aso[ot_aso["num_mismatch"] == 2].shape[0]
    num_tir_ot_2mm = ot_aso[(ot_aso["TIR"] == "TIR") & (ot_aso["num_mismatch"] == 2)].shape[0]
    # 1mm:
    num_tot_ot_1mm = ot_aso[ot_aso["num_mismatch"] == 1].shape[0]
    num_tir_ot_1mm = ot_aso[(ot_aso["TIR"] == "TIR") & (ot_aso["num_mismatch"] == 1)].shape[0]
    # 0mm:
    num_tot_ot_0mm = ot_aso[ot_aso["num_mismatch"] == 0].shape[0]
    num_tir_ot_0mm = ot_aso[(ot_aso["TIR"] == "TIR") & (ot_aso["num_mismatch"] == 0)].shape[0]

    # use Use pandas.concat instead of append here:
    # 3mm:
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot_3mm,
                                                 target_seq, 3]], columns=df_plot.columns)])
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in TIR regions", "start regions", num_tir_ot_3mm,
                                                 target_seq, 3]], columns=df_plot.columns)])
    # 2mm:
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot_2mm,
                                                    target_seq, 2]], columns=df_plot.columns)])
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in TIR regions", "start regions", num_tir_ot_2mm,
                                                    target_seq, 2]], columns=df_plot.columns)])
    # 1mm:
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot_1mm,
                                                    target_seq, 1]], columns=df_plot.columns)])
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in TIR regions", "start regions", num_tir_ot_1mm,
                                                    target_seq, 1]], columns=df_plot.columns)])
    # 0mm:
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot_0mm,
                                                    target_seq, 0]], columns=df_plot.columns)])
    df_plot = pd.concat([df_plot, pd.DataFrame([[aso_n, "OT in TIR regions", "start regions", num_tir_ot_0mm,
                                                    target_seq, 0]], columns=df_plot.columns)])

    # change output of output_df in rows. i is "ASO" column
    # 3mm:
    output_df.loc[output_df["ASO"] == i, "OT_tot_3mm"] = int(num_tot_ot_3mm)
    output_df.loc[output_df["ASO"] == i, "OT_TIR_3mm"] = int(num_tir_ot_3mm)
    # 2mm:
    output_df.loc[output_df["ASO"] == i, "OT_tot_2mm"] = int(num_tot_ot_2mm)
    output_df.loc[output_df["ASO"] == i, "OT_TIR_2mm"] = int(num_tir_ot_2mm)
    # 1mm:
    output_df.loc[output_df["ASO"] == i, "OT_tot_1mm"] = int(num_tot_ot_1mm)
    output_df.loc[output_df["ASO"] == i, "OT_TIR_1mm"] = int(num_tir_ot_1mm)
    # 0mm:
    output_df.loc[output_df["ASO"] == i, "OT_tot_0mm"] = int(num_tot_ot_0mm)
    output_df.loc[output_df["ASO"] == i, "OT_TIR_0mm"] = int(num_tir_ot_0mm)


output_df["OT_tot_3mm"] = output_df["OT_tot_3mm"].astype(int)
output_df["OT_TIR_3mm"] = output_df["OT_TIR_3mm"].astype(int)
output_df["OT_tot_2mm"] = output_df["OT_tot_2mm"].astype(int)
output_df["OT_TIR_2mm"] = output_df["OT_TIR_2mm"].astype(int)
output_df["OT_tot_1mm"] = output_df["OT_tot_1mm"].astype(int)
output_df["OT_TIR_1mm"] = output_df["OT_TIR_1mm"].astype(int)
output_df["OT_tot_0mm"] = output_df["OT_tot_0mm"].astype(int)
output_df["OT_TIR_0mm"] = output_df["OT_TIR_0mm"].astype(int)



# sort the ASO table by OT_TIR_0mm, then OT_TIR_1mm, then OT_TIR_2mm, then OT_TIR_3mm, then OT_tot_0mm, then OT_tot_1mm,
# then OT_tot_2mm, then OT_tot_3mm:
output_df = output_df.sort_values(by=["OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm", "OT_tot_0mm",
                                        "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm"], ascending=True)


# move input PNA to the top. use concat instead of append:
output_df = pd.concat([output_df[output_df["ASO"] == "input_PNA"], output_df[output_df["ASO"] != "input_PNA"]])
# reset row indexes:
output_df = output_df.reset_index(drop=True)

# make outputdf to the first 11 rows, if it has more than 11 rows:
if output_df.shape[0] > 11:
    #print(output_df.shape[0])
    # restrict output df to rows 1:11. dont use index,
    output_df = output_df.iloc[0:11, ]

# restrict df_plot to to only "ASO" rows that contain the ASOs in output_df:
df_plot = df_plot[df_plot["ASO"].isin(output_df["ASO"])]

# HMP microbiome off-targets (0mm exact matches)
if screen == "microbiome":
    hmp_file = sys.argv[1] + "/offtargets_microbiome_sorted.tab"
    if os.path.exists(hmp_file):
        hmp_ot = pd.read_table(hmp_file, sep='\t', index_col=False)
        kept_asos = set(output_df["ASO"])
        output_df["OT_HMP_0mm"] = 0
        for aso_id in hmp_ot["probe_id"].unique():
            if aso_id not in kept_asos:
                continue
            ot_aso = hmp_ot[hmp_ot["probe_id"] == aso_id]
            n_0mm = int((ot_aso["num_mismatch"] == 0).sum())
            output_df.loc[output_df["ASO"] == aso_id, "OT_HMP_0mm"] = n_0mm
            df_plot = pd.concat([df_plot, pd.DataFrame([[aso_id, "OT in HMP microbiome", "HMP microbiome",
                                                          n_0mm, ot_aso.iloc[0]["probe_seq"], 0]], columns=df_plot.columns)])
        output_df["OT_HMP_0mm"] = output_df["OT_HMP_0mm"].astype(int)
        hmp_ot.to_csv(sys.argv[1] + "/offtargets_hmp_sorted.csv", index=False)
        hmp_ot.to_excel(sys.argv[1] + "/offtargets_hmp_sorted.xlsx", index=False)

# Human transcriptome off-targets (0mm exact matches)
if screen == "human":
    human_file = sys.argv[1] + "/offtargets_human_sorted.tab"
    if os.path.exists(human_file):
        human_ot = pd.read_table(human_file, sep='\t', index_col=False)
        kept_asos = set(output_df["ASO"])
        output_df["OT_GRCh38_0mm"] = 0
        for aso_id in human_ot["probe_id"].unique():
            if aso_id not in kept_asos:
                continue
            ot_aso = human_ot[human_ot["probe_id"] == aso_id]
            n_0mm = int((ot_aso["num_mismatch"] == 0).sum())
            output_df.loc[output_df["ASO"] == aso_id, "OT_GRCh38_0mm"] = n_0mm
            df_plot = pd.concat([df_plot, pd.DataFrame([[aso_id, "OT in human transcriptome", "human transcriptome",
                                                          n_0mm, ot_aso.iloc[0]["probe_seq"], 0]], columns=df_plot.columns)])
        output_df["OT_GRCh38_0mm"] = output_df["OT_GRCh38_0mm"].astype(int)
        human_ot.to_csv(sys.argv[1] + "/offtargets_human_sorted.csv", index=False)
        human_ot.to_excel(sys.argv[1] + "/offtargets_human_sorted.xlsx", index=False)

# export df_plot for R plotting script
df_plot.to_csv(sys.argv[1] + "/df_plot.csv", index=False)

# create a heatmap visualizing the sequence (atgcs in diff. colors) for all the sequences in output_df["ASO_seq"].
# use seaborn for this.
# we have to start by creating a dataframe with 1 column per base and 1 row per ASO:

# create a np. array with 1 row per ASO and 1 column per base. For the list, give numbers +1 and 0 to bases matching the
# first row. Then, use seaborn to create a heatmap with the values in the array.
heatmap_list = []
heatmap_annot = []
first_aso = output_df["ASO_seq"][0]
for aso in output_df["ASO_seq"]:
    # give numbers +1 and 0 to bases matching the first_aso
    heatmap_list += [[1 if aso[i] == first_aso[i] else 0 for i in range(len(aso))]]
    heatmap_annot += [[i for i in aso]]

heatmap_array = np.array(heatmap_list)
heatmap_annotations = np.array(heatmap_annot)

# create a heatmap. adjust to the number of ASOs (height) and the number of bases (width):
figwidth = len(heatmap_array[0])/2+2
figheight = len(heatmap_array)/2+0.5

sns.set_theme(style="ticks")

fig, ax = plt.subplots(figsize=(figwidth, figheight))
ax = sns.heatmap(heatmap_array, linewidths=2, linecolor="white",
                    annot=heatmap_annotations, fmt="",
                    # cmap=["#ff681f", "#4682b4"],
                    #cmap=["darkred", "steelblue"],
                    cmap=["#14c8ff", "#002864"],
                    cbar=True)
# increase fontsize annotations & make bold:
for t in ax.texts:
    t.set_fontsize(14)
    t.set_weight("bold")
# addx ticks (make them black)
ax.set_xticks(np.arange(0.5, len(aso), 1))
# add line to separate first row:
ax.axhline(y=1, color='white', linewidth=8)
# annotate x labels with number of base:
ax.set_xticklabels(np.arange(1, len(aso)+1, 1), rotation=0, fontsize=14)
ax.set_xlabel("ASO position (5' to 3')", fontsize=18, fontweight="bold")
# annotate y labels with ASO name:
ax.set_yticks(np.arange(0.5, len(heatmap_array), 1))
ax.set_yticklabels(output_df["ASO"], rotation=0, fontsize=14)
# add cbar ticks and labels:
cbar = ax.collections[0].colorbar
cbar.set_ticks([0.25, 0.75])
cbar.set_ticklabels(["mismatch", "match"])
cbar.ax.tick_params(labelsize=14)
# bold
ax.set_ylabel("ASO name", fontsize=18, fontweight="bold")
plt.tight_layout()

# save figure:
plt.savefig(sys.argv[1] + '/heatmap.png', dpi=500)
plt.savefig(sys.argv[1] + '/heatmap.svg')

# save df for download/visualization:
output_df.to_csv(sys.argv[1] + "/result_table.csv", sep=",", index=False)
output_df.to_excel(sys.argv[1] + "/result_table.xlsx")
os.remove(sys.argv[1] + "/result_table.tsv")

