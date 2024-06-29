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
    mean_time = mean_time.loc[['small', 'medium', 'large']]

    return mean_time, std_time

def plot_subplot(ax, alpha_dir, beta_dir, gamma_dir, title):
    # Process each directory
    mean_time_alpha, std_time_alpha = process_file(alpha_dir)
    mean_time_beta, std_time_beta = process_file(beta_dir)
    mean_time_gamma, std_time_gamma = process_file(gamma_dir)

    # Plot the results on the provided axes
    ax.errorbar(mean_time_alpha.index, mean_time_alpha.values, yerr=std_time_alpha.values, capsize=5, marker='o', linestyle='-', color='blue', label=r'$\alpha$')
    ax.errorbar(mean_time_beta.index, mean_time_beta.values, yerr=std_time_beta.values, capsize=5, marker='v', linestyle='-', color='orange', label=r'$\beta$')
    ax.errorbar(mean_time_gamma.index, mean_time_gamma.values, yerr=std_time_gamma.values, capsize=5, marker='s', linestyle='-', color='green', label=r'$\gamma$')

    ax.set_title(title)
    ax.set_xlabel('Size')
    ax.set_ylabel('Mean Execution Time (ms)')
    ax.grid(True)
    ax.legend()

def main(parent_dir):
    subdirs = ["quasi-mcp-cpu", "quasi-mcp-cuda", "mcp-cpu", "qmcp-cpu"]
    titles = ["CPU quasi-MCP", "CUDA quasi-MCP", "CPU MCP", "CPU QMCP"]
    
    # Prepare subplots
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    axs = axs.flatten()
    
    for i, subdir in enumerate(subdirs):
        alpha_dir = os.path.join(parent_dir, "alpha", subdir)
        beta_dir = os.path.join(parent_dir, 'beta', subdir)
        gamma_dir = os.path.join(parent_dir, 'gamma', subdir)

        plot_subplot(axs[i], alpha_dir, beta_dir, gamma_dir, titles[i])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process CSV files from specified parent directory.')
    parser.add_argument('parent_dir', type=str, help='Parent directory containing subdirectories for each dataset')
    
    args = parser.parse_args()
    main(args.parent_dir)

