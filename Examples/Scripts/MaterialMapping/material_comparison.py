import matplotlib.pyplot as plt
import numpy as np
import uproot
import math

# Make a file list, e.g.
files = [
    (
        "geant4_material_tracks.root:material-tracks",
        "Geant4",
        "blue",
        "blue",
        "Blues",
        3,
    ),
    (
        "acts_itk_validation_tracks.root:material-tracks",
        "Acts (try all)",
        "orange",
        "none",
        "Oranges",
        4,
    ),
    (
        "acts_itk_propagation_tracks.root:material-tracks",
        "Acts (nav)",
        "green",
        "none",
        "Greens",
        5,
    ),
]

fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)

x_range = (-4.0, 4.0)  # [-math.pi, math.pi ] #
x_bins = 80
x_val_name = "v_eta"  # v_phi' #eta'
y_val_name = "t_X0"

y_lim = 0

y_ratio_values = []
y_ratio_errors = [0.0 for i in range(x_bins)]

for ifile, (fname, label, color, mfc, cmap, msize) in enumerate(files):

    # labels += [ label+" mean" ]
    tree = uproot.open(fname)

    # fig, ax = plt.subplots()
    y_arr = np.array(tree[y_val_name].array())
    x_arr = np.array(tree[x_val_name].array())

    # axs[ifile].hist2d(x_arr, y_arr, bins=x_bins, cmap=cmap, vmin=0.01)
    y_max = y_arr.max()
    y_lim = y_max if y_max > y_lim else y_lim

    # Generate the central bin values
    x_step = (x_range[1] - x_range[0]) / x_bins
    x_vals = [x_range[0] + (ix + 0.5) * x_step for ix in range(x_bins)]

    # Prepare the min /max
    y_min_vals = [1000.0] * x_bins
    y_max_vals = [0.0] * x_bins
    y_vals_sorted = [np.array([])] * x_bins

    for iv in range(len(x_arr)):
        x_b = int((x_arr[iv] - x_range[0]) / x_step)
        y_v = y_arr[iv]
        # Take min / max
        y_min_vals[x_b] = y_v if y_v < y_min_vals[x_b] else y_min_vals[x_b]
        y_max_vals[x_b] = y_v if y_v > y_max_vals[x_b] else y_max_vals[x_b]
        # Regulate the x value
        y_vals_sorted[x_b] = np.append(y_vals_sorted[x_b], y_v)

    axs[0].fill_between(
        x=x_vals,
        y1=y_min_vals,
        y2=y_max_vals,
        alpha=0.1,
        label=label + " spread",
        color=color,
    )
    axs[0].grid(axis="x")
    y_vals_mean = [y_bin_vals.mean() for y_bin_vals in y_vals_sorted]

    y_ratio_values += [y_vals_mean]

    y_vals_std = [y_bin_vals.std() for y_bin_vals in y_vals_sorted]
    y_vals_mse = [
        y_bin_vals.std() ** 2 / len(y_bin_vals) for y_bin_vals in y_vals_sorted
    ]

    axs[0].errorbar(
        x=x_vals,
        y=y_vals_mean,
        yerr=y_vals_mse,
        markersize=msize,
        marker="o",
        mfc=mfc,
        linestyle="none",
        label=label + " mean",
        color=color,
    )

    if ifile > 0:
        y_ratios = [
            y_ratio_values[ifile][ib] / y_ratio_values[0][ib] for ib in range(x_bins)
        ]
        axs[1].errorbar(
            x=x_vals,
            y=y_ratios,
            yerr=y_ratio_errors,
            markersize=msize,
            marker="o",
            mfc=mfc,
            linestyle="none",
            color=color,
            label=label,
        )
        axs[1].set_ylabel("ratio to Geant4")


axs[0].set_ylim(0.0, y_lim)
axs[0].grid(axis="x", linestyle="dotted")

axs[1].set_ylim(0.9, 1.1)
axs[1].grid(axis="x", linestyle="dotted")
axs[1].axhline(y=1.0, color="black", linestyle="-")

axs[0].legend(loc="upper center")
axs[1].legend(loc="upper center")


# Set the range of x-axis
plt.xlabel("η")
plt.show()
