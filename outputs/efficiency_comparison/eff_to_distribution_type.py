#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def process_file(directory):
    file_path = os.path.join(directory, 'results.csv')
    
    df = pd.read_csv(file_path, delimiter=';')

    df['suffix'] = df['name'].apply(lambda x: x.split('_')[-1])
    df['prefix'] = df['name'].apply(lambda x: '_'.join(x.split('_')[:-2])).apply(lambda x: '_'.join(x.split('_')[1:]).replace("_", " "))

    df_large = df[df['suffix'] == 'large']

    mean_time = df_large.groupby('prefix')['time_ms'].mean()
    std_time = df_large.groupby('prefix')['time_ms'].std()

    return mean_time, std_time

def plot_subplot(ax, dir, title, color, marker):
    mean_time, std_time = process_file(dir)

    ax.errorbar(mean_time.index, mean_time.values, yerr=std_time.values, capsize=5, marker=marker, linestyle='-', color=color, label=title)

def main(parent_dir):
    subdirs = ["quasi-mcp-cpu"]
    titles = ["CPU quasi-MCP", "CUDA quasi-MCP", "CPU MCP", "CPU QMCP"]
    colors = ["#003f5c", "#7a5195", "#ef5675", "#ffa600"]
    markers = ["o", "v", "s", "p"]
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    for i, subdir in enumerate(subdirs):
        dir = os.path.join(parent_dir, subdir)
        plot_subplot(ax, dir, titles[i], colors[i], markers[i])

    ax.set_title('Execution Time by Coverage Distribution Type (large only)')
    ax.set_xlabel('Coverage Distribution Type')
    ax.set_ylabel('Mean Execution Time (ms)')
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process CSV files from specified parent directory.')
    parser.add_argument('parent_dir', type=str, help='Parent directory containing subdirectories for each dataset')
    
    args = parser.parse_args()
    main(args.parent_dir)

