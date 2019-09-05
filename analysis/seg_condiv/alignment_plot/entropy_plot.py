"""Plot the alignments and their conservation scores by column."""

import math
import PyQt5.QtWidgets as Widgets
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
from Bio import AlignIO
from sys import exit

# CONSERVATION #
path = ''
gap_thresh = 0.5  # Maximum fraction of gaps in column before the column entropy is fixed at the value of gap_ent
gap_ent = 10  # Fixed "entropy" of gappy regions

MSA = AlignIO.read(path, 'fasta')
row_len = len(MSA[0, :])
column_len = len(MSA)

# Calculate distribution of symbols in each column
column_dists = []
for i in range(row_len):
    dist = {}
    for base in MSA[:, i]:  # Get counts
        dist[base] = dist.get(base, 0) + 1
    for base in dist:  # Convert counts to dist
        dist[base] = dist[base] / column_len
    column_dists.append(dist)

# Calculate entropy of each column
column_ents = []
for dist in column_dists:
    if dist.get('-', 0) > gap_thresh:
        ent = gap_ent
    else:
        ent = 0
        for base in dist:
            prob = dist[base]
            ent -= prob * math.log(prob, 2)
    column_ents.append(ent)

# PLOT #
width = 12
height = 12
scale = 30  # Factor by which to scale scores for visualization
colors_dark = {'A': '808080', 'C': 'e6e600', 'D': 'e60a0a', 'E': 'e60a0a', 'F': '3232aa',
               'G': 'c0c0c0', 'H': '8282d2', 'I': '0f820f', 'K': '145aff', 'L': '0f820f',
               'M': 'e6e600', 'N': '00dcdc', 'P': 'dc9682', 'Q': '00dcdc', 'R': '145aff',
               'S': 'fa9600', 'T': 'fa9600', 'V': '0f820f', 'W': 'b45ab4', 'Y': '3232aa',
               '-': 'ffffff', 'X': '000000'}
colors_lite = {'A': 'bfbfbf', 'C': 'f2f27f', 'D': 'f28484', 'E': 'f28484', 'F': '9898d4',
               'G': 'dfdfdf', 'H': 'c0c0e8', 'I': '86c086', 'K': '89acff', 'L': '86c086',
               'M': 'f2f27f', 'N': '7feded', 'P': 'edcac0', 'Q': '7feded', 'R': '89acff',
               'S': 'fcca7f', 'T': 'fcca7f', 'V': '86c086', 'W': 'd9acd9', 'Y': '9898d4',
               '-': 'ffffff', 'X': '7f7f7f'}
colors = colors_lite

# Window widgets
app = Widgets.QApplication([])
scene = Widgets.QGraphicsScene()
view = Widgets.QGraphicsView(scene)

# Style settings
pen = QtGui.QPen()
pen.setStyle(QtCore.Qt.NoPen)
brush = QtGui.QBrush()
brush.setStyle(QtCore.Qt.SolidPattern)
font = QtGui.QFont('Consolas', 7)

# Draw alignment
MSA = AlignIO.read(path, 'fasta')
for y, seq in enumerate(MSA):
    for x, char in enumerate(seq):
        rgb = tuple(int(colors[char][i:i+2], 16) for i in (0, 2, 4))
        brush.setColor(QtGui.QColor(*rgb))
        rect = scene.addRect(x * width, y * height, width, height, pen=pen, brush=brush)
        text = scene.addText(char, font=font)
        bound = text.boundingRect()
        text.setPos(x * width + (width - bound.width()) / 2, y * height + (height - bound.height()) / 2)

# Draw conservation scores
brush.setColor(QtGui.QColor(63, 63, 63))
y = (y + 1) * height
for i, ent in enumerate(column_ents):
    rect = scene.addRect(i * width, y, width, scale * ent, pen=pen, brush=brush)

view.show()
exit(app.exec_())
