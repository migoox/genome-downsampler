#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def process_file(directory):
    file_path = os.path.join(directory, 'results.csv')
    
    df = pd.read_csv(file_path, delimiter=';')

    mean_time = df.groupby('m')['time_ms'].mean()
    std_time = df.groupby('m')['time_ms'].std()

    return mean_time, std_time

def plot_influence_m(ax, directory, title, color, marker):
    mean_time, std_time = process_file(directory)
    
    ax.errorbar(mean_time.index, mean_time.values, yerr=std_time.values, marker=marker, capsize=5, linestyle='-', color=color, label=title)
    
    ax.set_title(title)
    ax.set_xlabel('m')
    ax.set_ylabel('Mean Execution Time (ms)')
    ax.legend()
    ax.grid(True)

def main(alpha_dir):
    subdirs = ["quasi-mcp-cpu", "quasi-mcp-cuda", "mcp-cpu", "qmcp-cpu"]
    titles = ["CPU quasi-MCP", "CUDA quasi-MCP", "CPU MCP", "CPU QMCP"]
    colors = ["#003f5c", "#7a5195", "#ef5675", "#ffa600"]
    markers = ["o", "v", "s", "p"]

    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    axs = axs.flatten()
    
    for i, subdir in enumerate(subdirs):
        dir_path = os.path.join(alpha_dir, subdir)
        plot_influence_m(axs[i], dir_path, titles[i], colors[i], markers[i])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot influence of m value on execution time from alpha directory.')
    parser.add_argument('alpha_dir', type=str, help='Directory containing subdirectories for each dataset')
    
    args = parser.parse_args()
    main(args.alpha_dir)

