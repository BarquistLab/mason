import matplotlib.pyplot as plt
import re
import pandas as pd
import seaborn as sns
import sys


all_off_targets = pd.read_table(sys.argv[1], sep='\t', index_col=False)

all_off_targets["ASO"] = all_off_targets["probe_id"].replace("_rev", "", regex=True)

# only use the mismatches with >+7 cons. mm:
all_off_targets = all_off_targets[all_off_targets["longest_stretch"] >= 7]


# add variable for whether it is in the TIR
all_off_targets.loc[all_off_targets["trans_coord"].isin(range(-20, 5)), "TIR"] = "TIR"
all_off_targets.loc[~all_off_targets["trans_coord"].isin(range(-20, 5)), "TIR"] = "not in TIR"

all_off_targets.to_csv(sys.argv[2] + "/offtargets_fulltranscripts_sorted.csv")
all_off_targets[all_off_targets["TIR"] == "TIR"].to_csv(sys.argv[2] + "/offtargets_startregions_sorted.csv")

df_plot = pd.DataFrame(columns=["ASO", "off-target type", "transcripts", "counts", "target sequence"])

# create dataframe:
for i in all_off_targets["ASO"].unique():
    target_seq = all_off_targets[all_off_targets["ASO"] == i].iloc[0, ]["probe_seq"]
    aso_n = re.sub(".*_(ASO.*)", "\\1", i)
    print(aso_n)
    ot_aso = all_off_targets[all_off_targets["ASO"] == i]
    num_tot_ot = ot_aso.shape[0]
    num_tir_ot = ot_aso[ot_aso["TIR"] == "TIR"].shape[0]
    df_plot = df_plot.append(pd.Series([aso_n, "OT in transcriptome", "whole transcriptome", num_tot_ot,
                                        target_seq], index=df_plot.columns), ignore_index=True)
    df_plot = df_plot.append(pd.Series([aso_n, "OT in TIR regions", "start regions", num_tir_ot,
                                        target_seq], index=df_plot.columns), ignore_index=True)


def show_on_single_plot(ax, maxlen):
    for p in ax.patches:
        _x = p.get_x() + p.get_width() / 2
        _y = p.get_y() + p.get_height() + maxlen/ 120
        value = '{:.0f}'.format(p.get_height())
        ax.text(_x, _y, value, ha="center")


def create_ot_barplot(dataframe, title, filepath):
    # get the number of asos:
    sns.set_style("whitegrid")
    plt.figure(figsize=(16, 11))
    plt.xlabel('xlabel', fontsize=16, fontweight='bold')
    plt.ylabel('ylabel', fontsize=16, fontweight='bold')
    plt.title('titlelabel', fontsize=16, fontweight='bold')

    bp = sns.barplot(x="ASO", y="counts", hue="off-target type", data=dataframe,
                     palette='Blues_r')
    bp.set(xlabel='ASO sequence', ylabel='Number of off-targets')
    show_on_single_plot(bp, maxlen=dataframe["counts"].max())
    bp.set_title(title, fontsize=20, fontweight='bold')
    bp.set_xticklabels(bp.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=14)
    bp.legend(title="off-target type", fontsize='large', title_fontsize=14)
    bp.figure.savefig(filepath)
    plt.clf()


create_ot_barplot(df_plot, "Critical off-targets of ASOs", sys.argv[2] + "/plot_ots_whole_transcriptome.png")

