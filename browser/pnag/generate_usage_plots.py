"""Generate usage plots from usage.log.

Produces two SVG plots in static/data/:
  - usage_plot_overall.svg  (total jobs per week, bar chart)
  - usage_plot_by_tool.svg  (jobs per week per tool, line chart)

Run daily via cron or on-demand.
"""

import os
from datetime import datetime

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
USAGE_LOG = os.path.join(SCRIPT_DIR, 'static/data/usage.log')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'static/data')

TOOL_COLORS = {
    'mason': 'steelblue',
    'scrambler': 'deeppink',
    'checker': 'limegreen',
}


def parse_usage_log(path):
    records = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split(' | ')
            if len(parts) >= 2:
                try:
                    dt = datetime.strptime(parts[0], '%Y-%m-%d %H:%M:%S')
                    tool = parts[1].strip()
                    records.append({'date': dt, 'tool': tool})
                except ValueError:
                    continue
    return pd.DataFrame(records)


def make_overall_plot(df):
    weekly = df.resample('W', on='date').size().reset_index(name='jobs')
    # Fill missing weeks with 0 using the same week boundaries as resample
    all_weeks = pd.date_range(start=weekly['date'].min(), end=weekly['date'].max(), freq='W')
    weekly = weekly.set_index('date').reindex(all_weeks, fill_value=0).reset_index()
    weekly.columns = ['date', 'jobs']
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.bar(weekly['date'], weekly['jobs'], width=5, color='steelblue', edgecolor='none', alpha=0.6)
    ax.set_title('Total MASON usage per week', fontsize=12)
    ax.set_xlabel('Week')
    ax.set_ylabel('Number of jobs')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45, ha='right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'usage_plot_overall.svg'), format='svg')
    plt.close(fig)


def make_tool_plot(df):
    # Create a complete weekly date range so weeks with 0 usage show as 0
    # Use resample boundaries (Sunday end-of-week) for consistent alignment
    all_weekly = df.resample('W', on='date').size().reset_index(name='_dummy')
    all_weeks = all_weekly['date']

    fig, ax = plt.subplots(figsize=(9, 3.5))
    for tool, color in TOOL_COLORS.items():
        tool_df = df[df['tool'] == tool]
        weekly = tool_df.resample('W', on='date').size().reset_index(name='jobs')
        # Reindex to full week range, fill missing weeks with 0
        weekly = weekly.set_index('date').reindex(all_weeks, fill_value=0).reset_index()
        weekly.columns = ['date', 'jobs']
        label = tool.upper() if tool == 'mason' else tool.capitalize()
        ax.plot(weekly['date'], weekly['jobs'], color=color, linewidth=2,
                marker='o', markersize=4, label=label, alpha=0.6)
    ax.set_title('Usage per tool per week', fontsize=12)
    ax.set_xlabel('Week')
    ax.set_ylabel('Number of jobs')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45, ha='right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper left', frameon=False)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'usage_plot_by_tool.svg'), format='svg')
    plt.close(fig)


def main():
    if not os.path.exists(USAGE_LOG):
        print(f"Usage log not found: {USAGE_LOG}")
        return

    df = parse_usage_log(USAGE_LOG)
    if df.empty:
        print("No usage data found.")
        return

    make_overall_plot(df)
    make_tool_plot(df)
    print("Usage plots generated.")


if __name__ == '__main__':
    main()
