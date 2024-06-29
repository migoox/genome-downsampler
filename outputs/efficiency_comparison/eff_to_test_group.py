#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def process_file(directory):
    # Construct the full file path
    file_path = os.path.join(directory, 'results.csv')
    
    # Load the CSV file
    df = pd.read_csv(file_path, delimiter=';')

    # Extract the suffix from the 'name' column
    df['suffix'] = df['name'].apply(lambda x: x.split('_')[-1])

    # Group by the suffix and calculate the mean of 'time_ms'
    mean_time = df.groupby('suffix')['time_ms'].mean()
    std_time = df.groupby('suffix')['time_ms'].std()

    # Ensure correct order of suffixes
    mean_time = mean_time.loc[['small', 'medium', 'large', 'real']]

    return mean_time, std_time

def plot_subplot(ax, dir, title, color, marker):
    # Process each directory
    mean_time, std_time = process_file(dir)

    # Plot the results on the provided axes
    ax.errorbar(mean_time.index, mean_time.values, yerr=std_time.values, capsize=5, marker=marker, linestyle='-', color=color, label=title)


def main(parent_dir):
    subdirs = ["quasi-mcp-cpu", "quasi-mcp-cuda", "mcp-cpu", "qmcp-cpu"]
    titles = ["CPU quasi-MCP", "CUDA quasi-MCP", "CPU MCP", "CPU QMCP"]
    colors = ["#003f5c", "#7a5195", "#ef5675", "#ffa600"]
    markers = ["o", "v", "s", "p"]
    
    # Prepare subplots
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))

    for i, subdir in enumerate(subdirs):
        dir = os.path.join(parent_dir, subdir)
        plot_subplot(ax, dir, titles[i], colors[i], markers[i])

    ax.set_title('')
    ax.set_xlabel('Test Group')
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

