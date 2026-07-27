import acts
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon, Circle
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection, LineCollection

def polyArea(face):
    x = [item[0] for item in face]
    y = [item[1] for item in face]
    # shoelace formula for computing area of polygons
    return 0.5 * abs(np.dot(x, np.roll(y,-1)) - np.dot(y, np.roll(x,-1)))


def plot(self, projection, ax=None):
    N = len(self.surfaces)
    poly_patches = [Polygon(face, closed=True) for face in self.surfaces]     # face = [[x1,y1], [x2,y2],...]   
    poly_collection = PatchCollection(poly_patches, alpha=0.5)
    poly_collection.set_facecolor(self.faceColors/255)
    
    # Plotting part 
    if ax is None:
        ax = plt.gca()
    ax.add_collection(poly_collection)
    
    if projection == "rz":
        ax.set_xlabel('z')
        ax.set_ylabel('r')

    if projection == "xy":
        ax.set_xlabel('x')
        ax.set_ylabel('y')

    
    for n_face, face in enumerate(self.surfaces):
        if polyArea(face) < 1e-8:
            ax.plot([item[0] for item in face], [item[1] for item in face], color=(self.faceColors[n_face][0]/255, self.faceColors[n_face][1]/255, self.faceColors[n_face][2]/255, 0.5), lw=1)


def plotTrack(self, projection, lineColor=None, linewidth=None, linestyle = None, ax=None):

    if lineColor == None:
        color = self.lineColor
    else:
        color = lineColor

    if linewidth == None:
        width = 1
    else: 
        width = linewidth

    if linestyle == None:
        style = 'solid'
    else:
        style = linestyle

    line_segments = [[segment[0], segment[1]] for segment in self.segments]

    line_collection = LineCollection(line_segments,linewidths=width, linestyles=style)
    line_collection.set_color(color/255)

    # Plotting part 
    if ax is None:
        ax = plt.gca()
    ax.add_collection(line_collection)
    
    if projection == "rz":
        ax.set_xlabel('z')
        ax.set_ylabel('r')

    if projection == "xy":
        ax.set_xlabel('x')
        ax.set_ylabel('y')


acts.PyVisualization.plot = plot
acts.PyVisualization.plotTrack = plotTrack