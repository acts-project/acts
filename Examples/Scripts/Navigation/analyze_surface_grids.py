import numpy as np
import matplotlib.pyplot as plt
import argparse
import json
import math


# Define a super class of an axis, a regular and a vraiable type
class Axis:
    def __init__(self, name, bins, range):
        self.name = name
        self.bins = bins
        self.range = range


class RegularAxis(Axis):
    def __init__(self, name, bins, range):
        super().__init__(name, bins, range)

    def get_edges(self):
        return np.linspace(self.range[0], self.range[1], self.bins + 1)


class VariableAxis(Axis):
    def __init__(self, name, edges: np.ndarray):
        super().__init__(name, len(edges) - 1, (edges[0], edges[-1]))
        self.edges = edges
        self.bins = len(edges) - 1

    def get_edges(self):
        return self.edges


def sort_vertices(policy_type: str, vertices: list):
    """Sorts the vertices in counter-clockwise order around their centroid."""
    if policy_descr == "Disc":
        return vertices

    # Calculate the centroid of the vertices
    centroid = [
        sum(v[0] for v in vertices) / len(vertices),
        sum(v[1] for v in vertices) / len(vertices),
    ]

    # Sort the vertices based on the angle from the centroid
    sorted_vertices = sorted(
        vertices, key=lambda v: math.atan2(v[1] - centroid[1], v[0] - centroid[0])
    )

    return sorted_vertices


def correct_vertices(
    policy_type: str,
    vertices: list,
    wrap_threshold=1.5 * math.pi,
    overlap_threshold=0.1,
):
    """Corrects the vertices for wrapping around the -pi to pi boundary in the phi coordinate."""
    if policy_type == "Plane":
        return vertices, []
    if policy_type == "Disc" or policy_type == "Ring":
        # bring vertices from (r,phi) to (x,y) for drawing
        corrected_vertices = []
        for v in vertices:
            x = v[0] * math.cos(v[1])
            y = v[0] * math.sin(v[1])
            corrected_vertices.append([x, y])
        return corrected_vertices, []
    # make a deep copy of the vertices
    corrected_vertices = vertices.copy()
    # sort them in the second coordinate
    corrected_vertices.sort(key=lambda v: v[1])
    # check if first and last are of bigger than the threshold
    diff_phi = abs(corrected_vertices[0][1] - corrected_vertices[-1][1])
    half_phi = 0.5 * abs(corrected_vertices[0][1] + corrected_vertices[-1][1])
    # sort then in the first coordinate
    corrected_vertices.sort(key=lambda v: v[0])
    diff_z = abs(corrected_vertices[0][0] - corrected_vertices[-1][0])

    if diff_phi > wrap_threshold:
        corrected_vertices = []
        mirror_vertices = []
        # we have a wrap situation
        first_wrap_point = False
        for v in vertices:
            if v[1] < -overlap_threshold * half_phi:
                # add 2pi to negative phi values
                corr_phi = v[1] + 2 * math.pi
                corrected_vertices.append([v[0], corr_phi])
                mirror_vertices.append(v)
                if first_wrap_point == False:
                    first_wrap_point = True
            else:
                corrected_vertices.append(v)
                mirror_vertices.append([v[0], v[1] - 2 * math.pi])
        return corrected_vertices, mirror_vertices
    else:
        return vertices, []


def plot_rectangular_grid(
    x: Axis, y: Axis, grid_data: np.ndarray, add_text=True, add_lines=True
):
    """Helper method to plot the rectangular grid given the two axes and the grid data."""
    x_edges = x.get_edges()
    y_edges = y.get_edges()

    # 3. Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # use a colormap with white for zero entries
    cmap = plt.cm.magma_r
    cmap.set_bad("white")

    # Add grid lines if requested
    if add_lines:
        for edge in x_edges:
            ax.axvline(edge, color="black", linestyle="--", linewidth=0.5)
        for edge in y_edges:
            ax.axhline(edge, color="black", linestyle="--", linewidth=0.5)

    # Plot the grid
    plt.hist2d(
        x=np.repeat((x_edges[:-1] + x_edges[1:]) / 2, y.bins),
        y=np.tile((y_edges[:-1] + y_edges[1:]) / 2, x.bins),
        bins=[x_edges, y_edges],
        weights=grid_data.flatten() if grid_data is not None else None,
        cmap=cmap,
    )

    # 3. Add text annotations to the bins
    if add_text and grid_data is not None:
        for i in range(len(x_edges) - 1):
            for j in range(len(y_edges) - 1):
                # Get the count value for the current bin
                count = grid_data[i, j]
                # Calculate the center position of the bin
                center_x = (x_edges[i] + x_edges[i + 1]) / 2
                center_y = (y_edges[j] + y_edges[j + 1]) / 2

                # Only add text if the count is greater than 0
                if count > 0:
                    # the color is black or white depending on the background color
                    text_color = "white" if count > 2 else "black"
                    # Place the text using plt.text()
                    ax.text(
                        center_x,
                        center_y,
                        int(count),
                        color=text_color,
                        ha="center",
                        va="center",
                        fontsize=8,
                        fontweight="bold",
                    )
    # Add labels and title
    ax.set_xlabel(f"{x.name}")
    ax.set_ylabel(f"{y.name}")

    return fig, ax


def plot_polar_grid(
    r: Axis, phi: Axis, grid_data: np.ndarray, add_text=True, add_lines=True
):
    """Helper method to plot the polar grid given the two axes and the grid data."""
    r_edges = r.get_edges()
    phi_edges = phi.get_edges()

    # use a colormap with white for zero entries
    cmap = plt.cm.magma_r
    cmap.set_bad("white")

    # Draw the polar grid sectors in cartesian coordinates
    for i in range(len(r_edges) - 1):
        for j in range(len(phi_edges) - 1):
            # Define the sector vertices
            sector_vertices = []
            sector_vertices.append(
                [
                    r_edges[i] * math.cos(phi_edges[j]),
                    r_edges[i] * math.sin(phi_edges[j]),
                ]
            )
            sector_vertices.append(
                [
                    r_edges[i + 1] * math.cos(phi_edges[j]),
                    r_edges[i + 1] * math.sin(phi_edges[j]),
                ]
            )
            sector_vertices.append(
                [
                    r_edges[i + 1] * math.cos(phi_edges[j + 1]),
                    r_edges[i + 1] * math.sin(phi_edges[j + 1]),
                ]
            )
            sector_vertices.append(
                [
                    r_edges[i] * math.cos(phi_edges[j + 1]),
                    r_edges[i] * math.sin(phi_edges[j + 1]),
                ]
            )
            # Create a polygon for the sector
            polygon = plt.Polygon(
                sector_vertices,
                closed=True,
                facecolor=(
                    cmap(grid_data[i, j] / np.nanmax(grid_data))
                    if grid_data is not None
                    else "white"
                ),
                edgecolor=None,
                alpha=1.0,
            )
            # Add the polygon to the plot
            plt.gca().add_patch(polygon)

    # Add grid lines if requested
    if add_lines:
        for edge in r_edges:
            circle = plt.Circle(
                (0, 0), edge, color="black", fill=False, linestyle="--", linewidth=0.5
            )
            ax.add_artist(circle)
        for edge in phi_edges:
            x = [r.range[0] * math.cos(edge), r.range[1] * math.cos(edge)]
            y = [r.range[0] * math.sin(edge), r.range[1] * math.sin(edge)]
            ax.plot(x, y, color="black", linestyle="--", linewidth=0.5)

    # Plot the grid
    # plt.hist2d(
    #    x=np.repeat(
    #        np.array([r_edge * math.cos(phi_edge + 0.5 * (phi_edges[1] - phi_edges[0])) for r_edge in (r_edges[:-1] + r_edges[1:]) / 2 for phi_edge in phi_edges[:-1]]),
    #        1,
    #    ),
    #    y=np.repeat(
    #        np.array([r_edge * math.sin(phi_edge + 0.5 * (phi_edges[1] - phi_edges[0])) for r_edge in (r_edges[:-1] + r_edges[1:]) / 2 for phi_edge in phi_edges[:-1]]),
    #        1,
    #    ),
    #    bins=[r_edges, phi_edges],
    #    weights=grid_data.flatten() if grid_data is not None else None,
    #    cmap=cmap,
    # )
    if add_text and grid_data is not None:
        for i in range(len(r_edges) - 1):
            for j in range(len(phi_edges) - 1):
                # Get the count value for the current bin
                count = grid_data[i, j]
                # Calculate the center position of the bin
                center_x = (
                    (r_edges[i] + r_edges[i + 1])
                    / 2
                    * math.cos((phi_edges[j] + phi_edges[j + 1]) / 2)
                )
                center_y = (
                    (r_edges[i] + r_edges[i + 1])
                    / 2
                    * math.sin((phi_edges[j] + phi_edges[j + 1]) / 2)
                )
                # Only add text if the count is greater than 0
                if count > 0:
                    # the color is black or white depending on the background color
                    text_color = "white" if count > 2 else "black"
                    # Place the text using plt.text()
                    ax.text(
                        center_x,
                        center_y,
                        int(count),
                        color=text_color,
                        ha="center",
                        va="center",
                        fontsize=8,
                        fontweight="bold",
                    )

    # Add labels and title
    ax.set_xlabel(f"x [mm]")
    ax.set_ylabel(f"y [mm]")
    ax.set_aspect("equal", adjustable="box")

    return fig, ax


# plot a surface as a Polygon
def plot_surface(vertices: np.ndarray, ax, fill_color):
    from matplotlib.patches import Polygon

    polygon = Polygon(
        vertices, closed=True, facecolor=fill_color, edgecolor=fill_color, alpha=0.25
    )
    ax.add_patch(polygon)


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument(
        "-p", "--policy", type=str, default="", help="Input JSON policy file"
    )
    p.add_argument("--no-text", action="store_true", help="Switch off bin text display")
    p.add_argument("--no-grid", action="store_true", help="Switch off grid display")
    p.add_argument(
        "--no-lines", action="store_true", help="Switch off grid lines display"
    )
    p.add_argument(
        "--no-surfaces", action="store_true", help="Switch off surface display"
    )
    p.add_argument(
        "--random-surface-color", action="store_true", help="Use random surface colors"
    )

    args = p.parse_args()

    if args.policy != "":
        with open(args.policy, "r") as f:
            policy_descr = json.load(f)
            type_descr = policy_descr["type"]
            grid_descr = policy_descr["grid"]
            axes_descr = grid_descr["axes"]
            # Lets define the type, possible types are plane, ring, disc, cylinder
            policy_type = None
            if "Plane" in type_descr:
                policy_type = "Plane"
            if "Disc" in type_descr:
                policy_type = "Disc"
            if "Ring" in type_descr:
                policy_type = "Ring"
            if "Cylinder" in type_descr:
                policy_type = "Cylinder"
            if policy_type is None:
                raise ValueError(f"Unknown grid type: {type_descr}")

            # Define default axes
            axes = []
            # Loop over the axes descriptions and replace the default axes
            for i, axis_descr in enumerate(axes_descr):
                axis_type = axis_descr["type"]
                if axis_type == "Equidistant":
                    axis_bins = axis_descr["bins"]
                    axis_range = axis_descr["range"]
                    # Get the axis values, regular for the moment
                    axes.append(RegularAxis("Axis_name", axis_bins, axis_range))

            if len(axes) == 1:
                if policy_type == "Ring":
                    reference_range = policy_descr["projectedReferenceRange"]
                    axes.insert(0, RegularAxis("r [mm]", 1, reference_range))
                grid_data = np.full((axes[0].bins, axes[1].bins), np.nan)
                if not args.no_grid:
                    grid_data_descr = grid_descr["data"]
                    for entry_descr in grid_data_descr:
                        bin_descr = entry_descr[0]
                        value = len(entry_descr[1])
                        grid_data[0, bin_descr[0] - 1] = value
            elif len(axes) == 2:
                # Create an empty grid data for demonstration
                grid_data = np.full((axes[0].bins, axes[1].bins), np.nan)
                if not args.no_grid:
                    grid_data_descr = grid_descr["data"]
                    for entry_descr in grid_data_descr:
                        bin_descr = entry_descr[0]
                        value = len(entry_descr[1])
                        grid_data[bin_descr[0] - 1, bin_descr[1] - 1] = value
            else:
                raise ValueError("Only 1D and 2D grids are supported in this example.")

            # Now plot the grid according to its type
            if policy_type == "Plane":
                axes[0].name = "x [mm]"
                axes[1].name = "y [mm]"
                fig, ax = plot_rectangular_grid(
                    axes[0],
                    axes[1],
                    grid_data,
                    add_text=not args.no_text,
                    add_lines=not args.no_lines,
                )
            elif policy_type == "Cylinder":
                axes[0].name = "phi [rad]"
                axes[1].name = "z [mm]"
                fig, ax = plot_rectangular_grid(
                    axes[0],
                    axes[1],
                    grid_data,
                    add_text=not args.no_text,
                    add_lines=not args.no_lines,
                )

            elif policy_type == "Disc" or policy_type == "Ring":

                # Make a cartesian view first where the polar grid will be sitting
                cart_axes = [
                    RegularAxis(
                        "x [mm]", axes[0].bins, (-axes[0].range[1], axes[0].range[1])
                    ),
                    RegularAxis(
                        "y [mm]", axes[1].bins, (-axes[0].range[1], axes[0].range[1])
                    ),
                ]
                empty_data = np.full((axes[0].bins, axes[1].bins), np.nan)
                fig, ax = plot_rectangular_grid(
                    cart_axes[0], cart_axes[1], empty_data, False, False
                )
                plot_polar_grid(
                    axes[0],
                    axes[1],
                    grid_data,
                    add_text=not args.no_text,
                    add_lines=not args.no_lines,
                )

            # Draw the projected surface points
            if not args.no_surfaces and "projectedSurfaces" in policy_descr:
                # plt.subplot(projection=None)
                for surface_vertices in policy_descr["projectedSurfaces"]:
                    # check special cylinder treatment
                    surface_vertices, wrap_vertices = correct_vertices(
                        policy_type,
                        surface_vertices,
                    )
                    color = "blue"
                    # Take a random color for each surface
                    if args.random_surface_color:
                        color = np.random.rand(
                            3,
                        )
                    plot_surface(
                        np.array(sort_vertices(policy_type, surface_vertices)),
                        ax,
                        fill_color=color,
                    )

                    # if len(projected_wrap) > 0:
                    #    plot_surface(np.array(sort_vertices(projected_wrap)), ax)

            # Add a colorbar to the plot
            if not args.no_grid:
                # Set range for color scale to min 1 to max value + 2
                plt.clim(0, np.nanmax(grid_data))
                plt.colorbar(label="Counts", ax=ax)
            plt.show()
