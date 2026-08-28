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

    elif char == "phi":
        x = vector3[0]
        y = vector3[1]
        r = np.sqrt(x**2 + y**2)
        alpha = np.arccos(abs(x) / r)
        if x >= 0 and y >= 0:
            return alpha
        elif x >= 0 and y < 0:
            return alpha + np.pi / 2
        elif x < 0 and y < 0:
            return alpha + np.pi
        else:
            return alpha + 3 / 2 * np.pi
    else:
        print("Not a valid projection")


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
                self._vis, track, context.simGeoContext
            )  # draw track not a free function

        return acts.examples.ProcessCode.SUCCESS


class PyVisualization2D(acts.VisualizationBuffer):

    def plot(
        self,
        projection,
        linewidth=None,
        linestyle=None,
        drawHitSensitives=True,
        **kwargs,
    ):
        import matplotlib.pyplot as plt

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

        # Applying optional ranges for plots
        surfaces = self.surfaces
        for k, val in kwargs.items():
            condition = lambda x: val[0] <= x <= val[1]

            # x range
            if k == "x":
                surfaces = [
                    surface
                    for surface in surfaces
                    if any(condition(v[0]) for v in surface)
                ]

            if k == "y":
                surfaces = [
                    surface
                    for surface in surfaces
                    if any(condition(v[1]) for v in surface)
                ]

            if k == "z":
                surfaces = [
                    surface
                    for surface in surfaces
                    if any(condition(v[2]) for v in surface)
                ]

            if k == "r":
                surfaces = [
                    surface
                    for surface in surfaces
                    if any(condition(np.sqrt(v[0] ** 2 + v[1] ** 2)) for v in surface)
                ]

            if k == "phi":
                surfaces = [
                    surface
                    for surface in surfaces
                    if any(condition(computeProjection("phi", v)) for v in surface)
                ]

        surfaces2D = [
            [
                [computeProjection(proj2D[0], v), computeProjection(proj2D[1], v)]
                for v in surface
            ]
        ]

        poly_patches = [
            Polygon(face, closed=True) for face in surfaces2D
        ]  # face = [[x1,y1], [x2,y2],...]
        poly_collection = PatchCollection(poly_patches, alpha=0.5)
        poly_collection.set_facecolor(self.faceColors / 255)
        poly_collection.set_edgecolor("black")

        # Plotting part
        ax.set_xlabel(proj2D[0])
        ax.set_ylabel(proj2D[1])
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set_xticks(np.linspace(-1000, 1000, 5))

        ax.add_collection(poly_collection)

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
            line_collection = LineCollection(
                line_segments, linewidths=width, linestyles=style
            )
            line_collection.set_color(self.lineColor / 255)
            ax.add_collection(line_collection)
            if drawHitSensitives == True:
                from matplotlib.path import Path

                hitSensitives = [
                    surface
                    for surface in surfaces2D
                    if any(
                        Path(surface).contains_point(segment[1])
                        for segment in line_segments
                    )
                ]
                hit_patches = [Polygon(face, closed=True) for face in hitSensitives]
                hit_collection = PatchCollection(hit_patches, alpha=0.5)
                hit_collection.set_facecolor("red")

                for n_face, face in enumerate(hitSensitives):
                    if polyArea(face) < 1e-8:
                        ax.plot(
                            [item[0] for item in face],
                            [item[1] for item in face],
                            color="red",
                            lw=1,
                        )

        ax.relim()
        ax.autoscale_view()
        fig.savefig(filename)
