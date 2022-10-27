import copy
import os
import matplotlib.pyplot as plt
import re
import pandas as pd
import seaborn as sns
import numpy as np
import sys
import six


other_ots = sys.argv[2]
print(other_ots)

colname = "none"
if other_ots == "human" or other_ots == "microbiome":
    if other_ots == "human":
        ot_screen = sys.argv[1] + "/offtargets_human_sorted.tab"
        screen_ot = pd.read_table(ot_screen, sep='\t', index_col=False)
        colname = "GRCh38"
        screen_ot["ASO"] = screen_ot["probe_id"].replace("_rev", "", regex=True)
        screen_ot.rename(columns={'locus_tag': 'transcript_name_NCBI', 'gene_name': 'description'}, inplace=True)
        screen_ot = screen_ot[screen_ot["longest_stretch"] >= 7]
        screen_ot.to_csv(sys.argv[1] + "/offtargets_human_sorted.csv")
        print("before excel")
        screen_ot.to_excel(sys.argv[1] + "/offtargets_human_sorted.xlsx")
        print("after excel")

    elif other_ots == "microbiome":
        ot_screen = sys.argv[1] + "/offtargets_microbiome_sorted.tab"
        screen_ot = pd.read_table(ot_screen, sep='\t', index_col=False)
        colname = "HMP"
        screen_ot["ASO"] = screen_ot["probe_id"].replace("_rev", "", regex=True)
        screen_ot.rename(columns={'locus_tag': 'organism_GenBank', 'gene_name': 'locus_tag'}, inplace=True)
        screen_ot = screen_ot[screen_ot["longest_stretch"] >= 7]
        screen_ot.to_csv(sys.argv[1] + "/offtargets_hmp_sorted.csv")
        screen_ot.to_excel(sys.argv[1] + "/offtargets_hmp_sorted.xlsx")


ot_table = sys.argv[1] + "/offtargets_fulltranscripts_sorted.tab"
all_off_targets = pd.read_table(ot_table, sep='\t', index_col=False)
all_off_targets["trans_coord"] = all_off_targets["trans_coord"] - 31


all_off_targets["ASO"] = all_off_targets["probe_id"].replace("_rev", "", regex=True)

# only use the mismatches with >+7 cons. mm:
all_off_targets = all_off_targets[all_off_targets["longest_stretch"] >= 7]



# add output df used for other things, e.g. Tm:
output_df = pd.read_csv(sys.argv[1] + "/result_table.tsv", sep="\t", index_col=None)

# add variable for whether it is in the TIR (-20 to +5 start)
all_off_targets.loc[all_off_targets["trans_coord"].isin(range(-20, 5)), "TIR"] = "TIR"
all_off_targets.loc[~all_off_targets["trans_coord"].isin(range(-20, 5)), "TIR"] = "not in TIR"

all_off_targets.to_csv(sys.argv[1] + "/offtargets_fulltranscripts_sorted.csv")
all_off_targets.to_excel(sys.argv[1] + "/offtargets_fulltranscripts_sorted.xlsx")

all_off_targets[all_off_targets["TIR"] == "TIR"].to_csv(sys.argv[1] + "/offtargets_startregions_sorted.csv")

df_plot = pd.DataFrame(columns=["ASO", "off-target type", "transcripts", "counts", "target sequence"])

# create dataframe:
for i in all_off_targets["ASO"].unique():
    target_seq = all_off_targets[all_off_targets["ASO"] == i].iloc[0, ]["probe_seq"]
    aso_n = re.sub(".*_(ASO.*)", "\\1", i)
    ot_aso = all_off_targets[all_off_targets["ASO"] == i]
    num_tot_ot = ot_aso.shape[0]-1
    num_tir_ot = ot_aso[ot_aso["TIR"] == "TIR"].shape[0]-1
    df_plot = df_plot.append(pd.Series([aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot,
                                        target_seq], index=df_plot.columns), ignore_index=True)
    df_plot = df_plot.append(pd.Series([aso_n, "OT in TIR regions", "start regions", num_tir_ot,
                                        target_seq], index=df_plot.columns), ignore_index=True)

    if other_ots == "human" or other_ots == "microbiome":
        ots_screen = screen_ot[screen_ot["ASO"] == i]
        num_scr_ot = ots_screen.shape[0]
        df_plot = df_plot.append(pd.Series([aso_n, "OT in " + colname, colname, num_scr_ot,
                                            target_seq], index=df_plot.columns), ignore_index=True)
        if other_ots == "human":
            output_df.loc[i, "OT_GRCh38"] = int(num_scr_ot)
        elif other_ots == "microbiome":
            output_df.loc[i, "OT_HMP"] = int(num_scr_ot)


    # change output:
    output_df.loc[ i, "OT_tot"] = int(num_tot_ot)
    output_df.loc[i, "OT_TIR"] = int(num_tir_ot)


# visualize:
def render_mpl_table(data, col_width=3.0, row_height=1, font_size=10,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=(25, data.shape[1]))
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0] % len(row_colors)])
    return ax


def show_on_single_plot(ax, maxlen):
    for p in ax.patches:
        _x = p.get_x() + p.get_width() / 2
        _y = p.get_y() + p.get_height() + maxlen/120
        value = '{:.0f}'.format(p.get_height())
        ax.text(_x, _y, value, ha="center", fontsize=15)


def create_ot_barplot(dataframe, title, filepath):
    # get the number of asos:
    sns.set_style("whitegrid")
    plt.figure(figsize=(16, 11))
    plt.xlabel('xlabel', fontsize=20, fontweight='bold')
    plt.ylabel('ylabel', fontsize=20, fontweight='bold')
    plt.title('titlelabel', fontsize=30, fontweight='bold')

    bp = sns.barplot(x="ASO", y="counts", hue="off-target type", data=dataframe,
                     palette='Blues_r')
    bp.set(xlabel='ASO sequence', ylabel='Number of off-targets')
    show_on_single_plot(bp, maxlen=dataframe["counts"].max())
    bp.set_title(title, fontsize=30, fontweight='bold')
    bp.set_xticklabels(bp.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=15)
    bp.legend(title="off-target type", fontsize=15, title_fontsize=20)
    bp.figure.savefig(filepath + ".png")
    bp.figure.savefig(filepath + ".svg")
    plt.clf()


create_ot_barplot(df_plot, "Critical off-targets of ASOs", sys.argv[1] + "/plot_ots_whole_transcriptome")


output_df["OT_tot"] = output_df["OT_tot"].astype(int)
output_df["OT_TIR"] = output_df["OT_TIR"].astype(int)


if colname == "none":
    output_df = output_df.drop(columns=["OT_HMP", "OT_cust", "OT_GRCh38"])
elif colname == "GRCh38":
    output_df = output_df.drop(columns=["OT_HMP", "OT_cust"])
    output_df["OT_GRCh38"] = output_df["OT_GRCh38"].astype(int)
elif colname == "HMP":
    output_df = output_df.drop(columns=["OT_cust", "OT_GRCh38"])
    output_df["OT_HMP"] = output_df["OT_HMP"].astype(int)
else:
    output_df = output_df.drop(columns=["OT_HMP", "OT_GRCh38"])

# add 1 if less than 20
# output_df["OT_TIR"] = np.where(output_df["location"].str.split(";", expand=True)[0].astype(int) < -20,
#                               output_df["OT_TIR"]+1, output_df["OT_TIR"])
# output_df["OT_tot"] = np.where(output_df["location"].str.split(";", expand=True)[0].astype(int) < -20,
#                               output_df["OT_tot"]+1, output_df["OT_tot"]+1)


output_df_fig = copy.deepcopy(output_df)
#output_df_fig["ASO_seq"] = output_df_fig["ASO_seq"].replace("^([ATGC]{10})([ATGC]+)", "\\1\n\\2", regex=True)
#output_df_fig["target_seq"] = output_df_fig["target_seq"].replace("^([AUGC]{10})([AUGC]+)", "\\1\n\\2", regex=True)

ax = render_mpl_table(output_df_fig, header_columns=0, col_width=4.0, row_height=4)
ax.figure.savefig(sys.argv[1] + "/result_table.png", bbox_inches='tight')
ax.figure.savefig(sys.argv[1] + "/result_table.svg", bbox_inches='tight')

# save df for download/visualization:
output_df.to_csv(sys.argv[1] + "/result_table.csv", sep=",")
output_df.to_excel(sys.argv[1] + "/result_table.xlsx")
os.remove(sys.argv[1] + "/result_table.tsv")

