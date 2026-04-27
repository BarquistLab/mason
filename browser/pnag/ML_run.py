# run ml model and prediction. load from pickle file
# also generate per-ASO SHAP force plots so the result table can link to them.
# The model predicts log2(MIC), so we relabel the SHAP plot axis (which shows
# the raw model output) with MIC values at log2-spaced tick positions.

import os
import pickle
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap


MIC_TICKS = (1.25, 2.5, 5, 10, 20, 40)

# SHAP hardcodes its red/blue palette in the matplotlib force plot. Replace only
# the two primary hex codes with darkorange / darkblue; leave any lighter shades
# (translucent fills, etc.) alone.
SHAP_COLOR_REMAP = {
    '#ff0d57': '#FF8C00',  # primary positive (red → darkorange)
    '#1e88e5': '#00008B',  # primary negative (blue → darkblue)
}


def recolor_shap_plot(fig):
    """Swap SHAP's two primary colors; keep everything else untouched."""
    def to_hex(color):
        try:
            return mcolors.to_hex(color).lower()
        except (ValueError, TypeError):
            return None

    for ax in fig.axes:
        for patch in ax.patches:
            new_fc = SHAP_COLOR_REMAP.get(to_hex(patch.get_facecolor()))
            if new_fc:
                patch.set_facecolor(new_fc)
            new_ec = SHAP_COLOR_REMAP.get(to_hex(patch.get_edgecolor()))
            if new_ec:
                patch.set_edgecolor(new_ec)
        for text in ax.texts:
            new_color = SHAP_COLOR_REMAP.get(to_hex(text.get_color()))
            if new_color:
                text.set_color(new_color)

# Shorter feature names for SHAP plots (other names are kept as-is)
FEATURE_RENAME = {
    'upec_tir_off_targets_1mm': 'OT_TIR_1mm',
    'purine_percentage': 'pur_%',
    'sc_bases': 'SC',
    'MFE_UPEC': 'MFE',
}


def relabel_axis_to_mic(fig):
    """Show MIC values on the log2(MIC) x-axis and convert prediction/base-value annotations."""
    log2_pos = np.log2(MIC_TICKS)
    num_re = re.compile(r'-?\d+\.\d+')
    pure_num_re = re.compile(r'^\s*-?\d+\.\d+\s*$')
    for ax in fig.axes:
        xlim = ax.get_xlim()
        visible = [(p, v) for p, v in zip(log2_pos, MIC_TICKS) if xlim[0] <= p <= xlim[1]]
        if visible:
            ax.set_xticks([p for p, _ in visible])
            ax.set_xticklabels([f'{v:g}' for _, v in visible])
        ax.set_xlabel('MIC (µM)', labelpad=6, fontsize=9, loc='left')
        ax.xaxis.set_label_position('top')
        # Convert the model-output annotations from log2 to MIC. Three cases:
        # - "f(x) = X.XX" / "base value = X.XX"     → substitute the number
        # - "X.XX" alone (the bold base/prediction) → also substitute
        # - "feature_name = X.XX"                   → leave alone
        for text in ax.texts:
            t = text.get_text()
            if pure_num_re.fullmatch(t):
                text.set_text(f'{2 ** float(t):.2f}')
            elif 'f(x)' in t or 'base value' in t.lower():
                new_t = num_re.sub(lambda m: f'{2 ** float(m.group(0)):.2f}', t)
                if new_t != t:
                    text.set_text(new_t)


filename = './pnag/static/rf_optimized_model_mason.sav'
rf = pickle.load(open(filename, 'rb'))

out_path = sys.argv[1]

data = pd.read_csv(out_path + "/saved_table_ml.csv")
asos = data["ASO"]
features = data.drop(columns=["ASO"])

y_pred = rf.predict(features)

# SHAP force plot per ASO — saved as PNG, linked from the result table.
# The model needs the original column names, but we display shorter ones in the plot.
explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(features)
features_display = features.rename(columns=FEATURE_RENAME)
shap_dir = os.path.join(out_path, 'shap')
os.makedirs(shap_dir, exist_ok=True)

mic_pred = np.exp2(y_pred)
for i, aso_name in enumerate(asos):
    shap.force_plot(
        explainer.expected_value,
        shap_values[i],
        features_display.iloc[i],
        matplotlib=True,
        show=False,
        figsize=(14, 3),
    )
    fig = plt.gcf()
    recolor_shap_plot(fig)
    relabel_axis_to_mic(fig)
    plt.savefig(os.path.join(shap_dir, f'shap_{aso_name}.svg'),
                bbox_inches='tight')
    plt.close('all')

features['MIC_pred'] = mic_pred
features.insert(0, "ASO", asos)
features.to_csv(out_path + "/saved_table_ml.csv", index=False)
