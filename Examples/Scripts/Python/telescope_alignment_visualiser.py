### Visualises eigenmodes of the alignment
# problem by drawing the original and distorted
# detector geometries in 2D projections.
# Intended for singular / weak mode diagnostics.
# Disclaimer: AI generated.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial.transform import Rotation

# ============================================================
# Configuration
# ============================================================

N_LAYERS = 9

# arbitrary scales for the layer spacing
# and distortions - used for visualisation
# purposes, adjust to taste.

LAYER_SPACING = 5.0

SENSOR_HALF_X = 2.0
SENSOR_HALF_Z = 2.0

TRANSLATION_SCALE = 2.0
ROTATION_SCALE = 30.0

DOF_NAMES = ["tx", "ty", "tz", "rx", "ry", "rz"]

# Projection axes: (horizontal coordinate, vertical coordinate)
PROJECTIONS = [
    ("x-y", (0, 1), "x [a.u.]", "y [a.u.]"),
    ("z-y", (2, 1), "z [a.u.]", "y [a.u.]"),
    ("x-z", (0, 2), "x [a.u.]", "z [a.u.]"),
]


# ============================================================
# 1. Convert sparse eigenvector to a 9 x 6 array
# ============================================================


def unpack_eigenvector(mode):
    """
    mode.components: eigenvector coefficients (floats)
    mode.labels:     1-based global parameter labels (1..54)

    Parameter order per layer: tx, ty, tz, rx, ry, rz.
    """
    labels = np.asarray(mode.labels, dtype=int) - 1
    values = np.asarray(mode.components, dtype=float)

    if len(labels) != len(values):
        raise ValueError("labels and components must have equal lengths")
    if np.any((labels < 0) | (labels >= 54)):
        raise ValueError("Parameter labels must be between 1 and 54")
    if len(np.unique(labels)) != len(labels):
        raise ValueError("Duplicate parameter labels")

    alignment = np.zeros((N_LAYERS, 6))
    alignment.ravel()[labels] = values
    return alignment


# ============================================================
# 2. Calculate nominal and shifted sensor planes
# ============================================================


def transform_plane(alignment, layer):
    """Return nominal and shifted 3D corners for one layer."""

    hx, hz = SENSOR_HALF_X, SENSOR_HALF_Z

    # Local sensor corners: sensor lies in the x-z plane.
    corners = np.array(
        [
            [-hx, 0, -hz],
            [hx, 0, -hz],
            [hx, 0, hz],
            [-hx, 0, hz],
        ],
        dtype=float,
    )

    y0 = layer * LAYER_SPACING
    centre = np.array([0.0, y0, 0.0])

    nominal = corners + centre

    # Translation and rotation are independently exaggerated
    # for visualisation; they are not physical corrections.
    translation = alignment[layer, :3] * TRANSLATION_SCALE
    angles = alignment[layer, 3:] * ROTATION_SCALE

    rotation = Rotation.from_euler("xyz", angles, degrees=True)

    shifted = rotation.apply(corners) + centre + translation

    return nominal, shifted


# ============================================================
# 3. Plot one eigenmode
# ============================================================


def plot_mode(mode, mode_index):
    alignment = unpack_eigenvector(mode)

    # Compute each plane only once; reuse for plotting and bounds.
    planes = [transform_plane(alignment, layer) for layer in range(N_LAYERS)]

    nominal_planes = [p[0] for p in planes]
    shifted_planes = [p[1] for p in planes]

    # Include every nominal and shifted corner when determining bounds.
    all_points = np.concatenate(nominal_planes + shifted_planes)
    margin = 1.0

    lower = all_points.min(axis=0) - margin
    upper = all_points.max(axis=0) + margin

    fig, axes = plt.subplots(1, 3, figsize=(17, 8))
    colors = plt.get_cmap("tab10")(np.arange(N_LAYERS))

    for ax, (name, (ix, iy), xlabel, ylabel) in zip(axes, PROJECTIONS):
        for layer, (nominal, shifted) in enumerate(planes):
            # Nominal geometry: dashed grey outlines.
            ax.add_patch(
                Polygon(
                    nominal[:, [ix, iy]],
                    closed=True,
                    fill=False,
                    edgecolor="0.55",
                    linestyle="--",
                    alpha=0.7,
                )
            )

            # Shifted geometry: solid layer-coloured outlines.
            ax.add_patch(
                Polygon(
                    shifted[:, [ix, iy]],
                    closed=True,
                    fill=False,
                    edgecolor=colors[layer],
                    linewidth=2,
                )
            )

            # Mark centres in the two longitudinal projections.
            if name != "x-z":
                ax.plot(
                    nominal[:, [ix, iy]].mean(axis=0)[0],
                    nominal[:, [ix, iy]].mean(axis=0)[1],
                    "o",
                    color="0.55",
                    markersize=3,
                )
                ax.plot(
                    shifted[:, [ix, iy]].mean(axis=0)[0],
                    shifted[:, [ix, iy]].mean(axis=0)[1],
                    "o",
                    color=colors[layer],
                    markersize=4,
                )

        ax.set(
            xlabel=xlabel,
            ylabel=ylabel,
            title=f"Projection: {name}",
        )
        ax.grid(True, alpha=0.2)
        ax.set_aspect("auto")

        # Projection-specific limits, derived from all corners.
        ax.set_xlim(lower[ix], upper[ix])
        ax.set_ylim(lower[iy], upper[iy])
        ax.set_autoscale_on(False)

    # Indicate beam direction in the longitudinal views.
    for ax in axes[:2]:
        ax.annotate(
            "Beam (+y)",
            xy=(0.95, 0.97),
            xycoords="axes fraction",
            ha="right",
            va="top",
        )

    # Legend.
    from matplotlib.lines import Line2D

    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="0.55",
                linestyle="--",
                label="Nominal geometry",
            ),
            Line2D(
                [0],
                [0],
                color="black",
                linestyle="-",
                linewidth=2,
                label="Shifted geometry",
            ),
        ],
        loc="lower center",
        ncol=2,
        frameon=False,
    )

    # Largest eigenvector coefficients.
    flat = alignment.ravel()
    dominant_indices = np.argsort(np.abs(flat))[::-1][:4]

    dominant = ", ".join(
        f"L{i // 6 + 1}:{DOF_NAMES[i % 6]}={flat[i]:+.3f}" for i in dominant_indices
    )

    eigenvalue = getattr(mode, "EV", float("nan"))

    fig.suptitle(
        f"Alignment eigenmode {mode_index + 1} | EV = {eigenvalue:.4g}\n"
        f"Dominant coefficients: {dominant}",
        fontsize=11,
    )

    fig.tight_layout(rect=[0, 0.08, 1, 0.89])

    fig.savefig(
        f"alignment_mode_{mode_index + 1:03d}.png",
        dpi=180,
        bbox_inches="tight",
    )


# ============================================================
# 4. Plot all eigenmodes
# ============================================================


def plot_all_modes(eigenvectors):
    for i, mode in enumerate(eigenvectors):
        plot_mode(mode, i)
