#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import os

if len(sys.argv) != 3:
    print("Usage: ./cmp.py parent_dir filename")
    sys.exit(1)

parent_dir = sys.argv[1]
filename = sys.argv[2]
subdirs = ["quasi-mcp-cpu", "mcp-cpu"]
titles = ["CPU quasi-MCP", "CPU MCP"]
colors = ["#003f5c", "#ef5675"]

data = []

for subdir in subdirs:
    filepath = os.path.join(parent_dir, subdir, filename)
    x = []
    y1 = []
    y2 = []
    with open(filepath, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            x.append(float(cols[0]))
            y1.append(float(cols[1]))
            y2.append(float(cols[2]))
    data.append((x, y1, y2))

fig, axs = plt.subplots(2, 1, figsize=(14, 10))

for i, ax in enumerate(axs.flat):
    x, y1, y2 = data[i]
    ax.plot(x, y1, label='Input')
    ax.plot(x, y2, label=f'After Downsampling with {titles[i]}', color=colors[i])
    ax.set_xlabel('Reference Genome Index')
    ax.set_ylabel('Coverage')
    ax.set_title(titles[i])
    ax.legend()
    ax.grid(True)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()

