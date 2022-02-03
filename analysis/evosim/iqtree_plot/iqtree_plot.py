"""Plot results from fitting evolutionary parameters."""

import os
import re
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse

order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
labels = [('disorder_0%', '0red_D'), ('disorder_50%', '50red_D'), ('disorder_100%', '100red_D'),
          ('order_0%', '0red_O'), ('order_50%', '50red_O'), ('order_100%', '100red_O')]
Record = namedtuple('Record', ['plot_label', 'file_label', 'matrix', 'freqs'])

# Read WAG model
with open('../config/WAG.txt') as file:
    # Parse exchangeability matrix
    WAG_matrix = np.zeros((20, 20))
    for i in range(19):
        line = file.readline()
        for j, value in enumerate(line.split()):
            WAG_matrix[i + 1, j] = float(value)
            WAG_matrix[j, i + 1] = float(value)

    # Parse equilibrium frequencies
    for _ in range(2):
        line = file.readline()
    WAG_freqs = np.array([float(value) for value in line.split()])
rate = (WAG_freqs * (WAG_freqs * WAG_matrix).sum(axis=1)).sum()
WAG_matrix = WAG_matrix / rate  # Normalize average rate to 1

# Read IQ-TREE matrices
records = []
for plot_label, file_label in labels:
    with open(f'../iqtree_fit/out/{file_label}.iqtree') as file:
        # Move to exchangeability matrix and parse
        line = file.readline()
        while not line.startswith('Substitution parameters'):
            line = file.readline()
        for _ in range(2):
            line = file.readline()

        rows = []
        while line != '\n':
            rows.append([float(value) for value in line.split()])
            line = file.readline()

        # Move to equilibrium frequencies and parse
        for _ in range(3):
            line = file.readline()

        syms, freqs = [], []
        while line != '\n':
            match = re.search('pi\(([A-Z])\) = (0.[0-9]+)', line)
            syms.append(match.group(1))
            freqs.append(float(match.group(2)))
            line = file.readline()
        freqs = np.array(freqs)

        if syms != order:
            raise RuntimeError('Symbols in matrix are not in expected order.')

    # Make matrix and scale
    matrix = np.zeros((len(syms), len(syms)))
    for i, row in enumerate(rows[:-1]):
        for j, value in enumerate(row):
            matrix[i+1, j] = value
            matrix[j, i+1] = value

    rate = (freqs * (freqs * matrix).sum(axis=1)).sum()
    matrix = matrix / rate
    records.append(Record(plot_label, file_label, matrix, freqs))

# Make plots
if not os.path.exists('out/'):
    os.mkdir('out/')

# 1 HEATMAP
fig, axs = plt.subplots(2, 3, figsize=(8, 6))
vmax = max([record.matrix.max() for record in records])
for ax, record in zip(axs.ravel(), records):
    ax.imshow(record.matrix, vmax=vmax)
    ax.set_title(record.plot_label)
    ax.set_xticks(range(20))
    ax.set_xticklabels(order, fontdict={'fontsize': 7})
    ax.set_yticks(range(20))
    ax.set_yticklabels(order, fontdict={'fontsize': 7})
fig.colorbar(ScalarMappable(Normalize(0, vmax)), ax=axs, fraction=0.025)
plt.savefig('out/heatmap_all.png', bbox_inches='tight')
plt.close()

# 2 DOT PLOT
scale = 0.0325
fig, axs = plt.subplots(2, 3, figsize=(8, 6))
for ax, record in zip(axs.ravel(), records):
    ax.set_title(record.plot_label)
    ax.set_xlim(0, 20)
    ax.set_xticks(range(1, 21))
    ax.set_xticklabels(order, fontdict={'fontsize': 7})
    ax.set_ylim(0, 20)
    ax.set_yticks(range(1, 21))
    ax.set_yticklabels(order[::-1], fontdict={'fontsize': 7})
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_aspect(1)

    for i, row in enumerate(record.matrix):
        for j, value in enumerate(row):
            if j >= i:
                continue
            c = Circle((j+1, 20-i), scale*value, color='black')
            ax.add_patch(c)
plt.savefig('out/dotplot_all.png', bbox_inches='tight')
plt.close()

# 3 DOT COMPARISONS
scale = 0.03
pairs = [(records[1].plot_label, records[1].file_label, records[1].matrix, 'black', 'WAG', 'WAG', WAG_matrix, 'white'),
         (records[4].plot_label, records[4].file_label, records[4].matrix, 'grey', 'WAG', 'WAG', WAG_matrix, 'white'),
         (records[1].plot_label, records[1].file_label, records[1].matrix, 'black', records[4].plot_label, records[4].file_label, records[4].matrix, 'grey')]
for pl1, fl1, m1, c1, pl2, fl2, m2, c2 in pairs:
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.set_xlim(0, 20)
    ax.set_xticks(range(1, 21))
    ax.set_xticklabels(order, fontdict={'fontsize': 8})
    ax.set_ylim(0, 20)
    ax.set_yticks(range(1, 21))
    ax.set_yticklabels(order[::-1], fontdict={'fontsize': 8})
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_aspect(0.5)  # Scale vertical axis half of horizontal axis

    for i, row in enumerate(m1):
        for j, value in enumerate(row):
            if j >= i:
                continue
            c = Ellipse((j + 1, 20 - i), height=2*scale*value, width=scale*value, facecolor=c1, edgecolor='black')
            ax.add_patch(c)
    for i, row in enumerate(m2):
        for j, value in enumerate(row):
            if j >= i:
                continue
            c = Ellipse((j + 1.5, 20 - i), height=2*scale*value, width=scale*value, facecolor=c2, edgecolor='black')
            ax.add_patch(c)
    handles = [Line2D([], [], label=pl1, marker='o', markerfacecolor=c1, markeredgecolor='black', markersize=8, linestyle='None'),
               Line2D([], [], label=pl2, marker='o', markerfacecolor=c2, markeredgecolor='black', markersize=8, linestyle='None')]
    ax.legend(handles=handles, bbox_to_anchor=(1, 0.5), loc='center left')

    plt.savefig(f'out/dotplot_adj_{fl1}-{fl2}.png', bbox_inches='tight')
    plt.close()

# 4 RATIO DOT COMPARISONS
scale = 0.15
pairs = [(records[1].plot_label, records[1].file_label, records[1].matrix, 'WAG', 'WAG', WAG_matrix),
         (records[4].plot_label, records[4].file_label, records[4].matrix, 'WAG', 'WAG', WAG_matrix),
         (records[1].plot_label, records[1].file_label, records[1].matrix, records[4].plot_label, records[4].file_label, records[4].matrix)]
for pl1, fl1, m1, pl2, fl2, m2 in pairs:
    fig, ax = plt.subplots()
    ax.set_title(f'log10 ratio of {pl1} to {pl2}')
    ax.set_xlim(0, 20)
    ax.set_xticks(range(1, 21))
    ax.set_xticklabels(order, fontdict={'fontsize': 8})
    ax.set_ylim(0, 20)
    ax.set_yticks(range(1, 21))
    ax.set_yticklabels(order[::-1], fontdict={'fontsize': 8})
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_aspect(1)

    for i, row in enumerate(np.log10(m1/m2)):
        for j, value in enumerate(row):
            if j >= i:
                continue
            c = Circle((j+1, 20-i), scale*abs(value), facecolor='black' if value > 1 else 'white', edgecolor='black')
            ax.add_patch(c)

    # Make legend
    y1, y2 = 10, 15
    rs = [1, 2, 4]
    dy = (y2 - y1 - sum([2*scale*r for r in rs])) / len(rs)
    ys, y = [], y1
    for r in rs:
        ys.append(y + scale*r)
        y += 2*scale*r + dy
    for r, y in zip(rs, ys):
        ax.add_patch(Circle((22, y), scale*r, facecolor='grey', edgecolor='black', clip_on=False))
        ax.text(23, y, f'|log(R)| = {r}', size=8, va='center', clip_on=False)

    ax.add_patch(Circle((22, 8), 2.5*scale, facecolor='black', edgecolor='black', clip_on=False))
    ax.text(23, 8, 'R > 1', size=8, va='center', clip_on=False)
    ax.add_patch(Circle((22, 6.5), 2.5*scale, facecolor='white', edgecolor='black', clip_on=False))
    ax.text(23, 6.5, 'R < 1', size=8, va='center', clip_on=False)

    plt.savefig(f'out/dotplot_ratio_{fl1}-{fl2}.png')
    plt.close()

# 5 DIFFERENCE DOT COMPARISONS
scale = 0.06
pairs = [(records[1].plot_label, records[1].file_label, records[1].matrix, 'WAG', 'WAG', WAG_matrix),
         (records[4].plot_label, records[4].file_label, records[4].matrix, 'WAG', 'WAG', WAG_matrix),
         (records[1].plot_label, records[1].file_label, records[1].matrix, records[4].plot_label, records[4].file_label, records[4].matrix)]
for pl1, fl1, m1, pl2, fl2, m2 in pairs:
    fig, ax = plt.subplots()
    ax.set_title(f'Difference of {pl1} and {pl2}')
    ax.set_xlim(0, 20)
    ax.set_xticks(range(1, 21))
    ax.set_xticklabels(order, fontdict={'fontsize': 8})
    ax.set_ylim(0, 20)
    ax.set_yticks(range(1, 21))
    ax.set_yticklabels(order[::-1], fontdict={'fontsize': 8})
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_aspect(1)

    for i, row in enumerate(m1-m2):
        for j, value in enumerate(row):
            if j >= i:
                continue
            c = Circle((j+1, 20-i), scale*abs(value), facecolor='black' if value > 0 else 'white', edgecolor='black')
            ax.add_patch(c)

    # Make legend
    y1, y2 = 10, 15
    rs = [1, 4, 7]
    dy = (y2 - y1 - sum([2*scale*r for r in rs])) / len(rs)
    ys, y = [], y1
    for r in rs:
        ys.append(y + scale*r)
        y += 2*scale*r + dy
    for r, y in zip(rs, ys):
        ax.add_patch(Circle((22, y), scale*r, facecolor='grey', edgecolor='black', clip_on=False))
        ax.text(23, y, f'|D| = {r}', size=8, va='center', clip_on=False)

    ax.add_patch(Circle((22, 8), 4*scale, facecolor='black', edgecolor='black', clip_on=False))
    ax.text(23, 8, 'D > 0', size=8, va='center', clip_on=False)
    ax.add_patch(Circle((22, 6.5), 4*scale, facecolor='white', edgecolor='black', clip_on=False))
    ax.text(23, 6.5, 'D < 0', size=8, va='center', clip_on=False)

    plt.savefig(f'out/dotplot_diff_{fl1}-{fl2}.png')
    plt.close()

# 6 FREQUENCIES
width = 0.2
bars = [(records[1].plot_label, records[1].freqs, 'black'),
        (records[4].plot_label, records[4].freqs, 'grey'),
        ('WAG', WAG_freqs, 'white')]
plt.figure(figsize=(8, 4))
for i, (plot_label, freqs, color) in enumerate(bars):
    dx = -len(bars)//2 + i + (1.5 if len(bars)%2 == 0 else 1)
    plt.bar([x+width*dx for x in range(len(order))], freqs, label=plot_label, facecolor=color, edgecolor='black', width=width)
plt.xticks(range(len(order)), order)
plt.xlabel('Amino acid')
plt.ylabel('Frequency')
plt.legend()
plt.savefig('out/bar_freqs.png')
plt.close()

"""
DEPENDENCEIES
../iqtree_fit/iqtree_fit.py
    ../iqtree_fit/out/*.iqtree
../config/WAG.txt
"""