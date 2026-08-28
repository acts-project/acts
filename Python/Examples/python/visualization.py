import acts
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection


def polyArea(face):
    x = [item[0] for item in face]
    y = [item[1] for item in face]
    # shoelace formula for computing area of polygons
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def computeProjection(char, vector3):
    if char == "x":
        return vector3[0]
    elif char == "y":
        return vector3[1]
    elif char == "z":
        return vector3[2]
    elif char == "r":
        return np.sqrt(vector3[0] ** 2 + vector3[1] ** 2)
    else:
        print("Not a valid projection")


def interpolateCurve(points):
    """
    points has to be of the type np.array()[[x1,y1], [x2,y2], ...]
    """
    from scipy.interpolate import make_interp_spline

    x = points[:, 0]
    y = points[:, 1]

    # Parameterize by vertex index
    t = np.arange(len(points))

    # Cubic spline through all vertices
    spline_x = make_interp_spline(t, x, k=3)
    spline_y = make_interp_spline(t, y, k=3)

    # Evaluate densely
    tt = np.linspace(t[0], t[-1], 500)
    curve_x = spline_x(tt)
    curve_y = spline_y(tt)

    return curve_x, curve_y


# algorithm to visualize trackk
class TrackVisualizerAlg(acts.examples.IAlgorithm):
    def __init__(self, name, level, vis):
        acts.examples.IAlgorithm.__init__(self, name, level)

        self._vis = vis
        self.tracks = acts.examples.ReadDataHandle(
            self, acts.examples.ConstTrackContainer, "Tracks"
        )
        self.tracks.initialize("tracks")

    def execute(self, context):
        tracks = self.tracks(context.eventStore)
        print(f"Event {context.eventNumber}: {len(tracks)} tracks")
        for track in tracks:
            acts.EventDataView3D.drawTrack(
                self._vis, track, context.geoContext
            )  # draw track not a free function

        return acts.examples.ProcessCode.SUCCESS


class PyVisualization2D(acts.VisualizationBuffer):

    def plot(
        self, projection, filename, interpolate=False, linewidth=None, linestyle=None
    ):

        # Reduce font size and complexity to avoid raster overflow
        plt.rcParams["font.size"] = 8
        plt.rcParams["figure.figsize"] = (12, 10)

        fig, ax = plt.subplots()

        if linewidth == None:
            width = 1
        else:
            width = linewidth

        if linestyle == None:
            style = "solid"
        else:
            style = linestyle

        proj2D = [char for char in projection]

        if len(proj2D) > 2:
            print("Only 2D projection supported")

        surfaces2D = [
            [
                [computeProjection(proj2D[0], v), computeProjection(proj2D[1], v)]
                for v in surface
            ]
            for surface in self.surfaces
        ]

        poly_patches = [
            Polygon(face, closed=True) for face in surfaces2D
        ]  # face = [[x1,y1], [x2,y2],...]
        poly_collection = PatchCollection(poly_patches, alpha=0.5)
        poly_collection.set_facecolor(self.faceColors / 255)

        ax.set_xlabel(proj2D[0])
        ax.set_ylabel(proj2D[1])
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_xticks(np.linspace(-1000, 1000, 5))
        ax.set_yticks(np.linspace(-1000, 1000, 5))

        ax.add_collection(poly_collection)

        ax.autoscale()

        for n_face, face in enumerate(surfaces2D):
            if polyArea(face) < 1e-8:
                ax.plot(
                    [item[0] for item in face],
                    [item[1] for item in face],
                    color=(
                        self.faceColors[n_face][0] / 255,
                        self.faceColors[n_face][1] / 255,
                        self.faceColors[n_face][2] / 255,
                        0.5,
                    ),
                    lw=1,
                )

        # Check if there is a track to be drawn
        print(len(self.segments))
        if len(self.segments) != 0:
            line_segments = [
                [
                    [
                        computeProjection(proj2D[0], segment[0]),
                        computeProjection(proj2D[1], segment[0]),
                    ],
                    [
                        computeProjection(proj2D[0], segment[1]),
                        computeProjection(proj2D[1], segment[1]),
                    ],
                ]
                for segment in self.segments
            ]

            if interpolate == True:
                line_segments = np.asarray(line_segments)
                points = line_segments[:, 0:]
                shape = points.shape
                points = points.reshape(shape[0] * shape[1], -1)
                print(points)
                interp_x, interp_y = interpolateCurve(points)
                ax.plot(
                    interp_x, interp_y, color=self.lineColor[0] / 255, linewidth=width
                )

            else:
                line_collection = LineCollection(
                    line_segments, linewidths=width, linestyles=style
                )
                line_collection.set_color(self.lineColor / 255)
                ax.add_collection(line_collection)

        ax.relim()
        ax.autoscale_view()
        fig.savefig(filename)
