# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# python includes
from collections import namedtuple
import math
import numpy as np

# python based plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import ticker
import matplotlib.style as style
from mpl_toolkits.axes_grid1 import make_axes_locatable

# detray imports
from .plot_helpers import plt_data, axis_options, legend_options

style.use("tableau-colorblind10")
# style.use('seaborn-colorblind')

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.size": 25,
        "font.family": "serif",
    }
)


# See: https://stackoverflow.com/questions/42142144/displaying-first-decimal-digit-in-scientific-notation-in-matplotlib
class ScalarFormatterForceFormat(ticker.ScalarFormatter):
    def _set_format(self):
        self.format = "%3.1f"


# ------------------------------------------------------------------------------
# Global identifiers
# ------------------------------------------------------------------------------

""" Default color for graphs and histograms """
default_color = "tab:blue"

# ------------------------------------------------------------------------------
# Data Plotting
# ------------------------------------------------------------------------------

""" Plotter interface that uses pyplot/matplotlib. """


class pyplot_factory:

    def __init__(self, out_dir, logger, atlas_badge=""):
        self.name = ("Pyplot",)
        self.output_prefix = out_dir
        self.logger = logger
        self.atlas_badge = atlas_badge
        self.badge_scale = 1.1
        self.axis_formatter = ScalarFormatterForceFormat()
        self.axis_formatter.set_powerlimits((-2, 2))

    # Add legend to a plot. Labels must be defined.
    def __add_legend(self, ax, options=legend_options()):
        return ax.legend(
            title=options.title,
            loc=options.loc,
            bbox_to_anchor=(options.horiz_anchor, options.vert_anchor),
            ncol=options.ncol,
            borderpad=0.3,
            columnspacing=options.colspacing,
            handletextpad=options.handletextpad,
        )

    # Adjust label spacing in legend
    def __adjust_lgd_label_spacing(self, lgd):
        # Refine legend
        lgd.legend_handles[0].set_visible(False)
        for handle in lgd.legend_handles[1:]:
            handle.set_sizes([40])

        # Adjust spacing in box
        for vpack in lgd._legend_handle_box.get_children()[:1]:
            for hpack in vpack.get_children():
                hpack.get_children()[0].set_width(0)

    # Update after adding new entry to existing legend
    def __update_legend(self, lgd):
        handles, labels = lgd.axes.get_legend_handles_labels()
        lgd._legend_box = None
        lgd._init_legend_box(handles, labels)
        lgd._set_loc(lgd._loc)
        lgd.set_title(lgd.get_title().get_text())

    # Find the axis boundaries either from data or custom boundaries
    def __get_axis_boundaries(self, data, axis_opts):
        if axis_opts.min is not None and axis_opts.max is not None:
            return axis_opts.min, axis_opts.max
        else:
            return np.min(data), np.max(data)

    # Apply boundary to input data
    def __apply_boundary(self, data, min_v, max_v):
        if min_v is not None and max_v is not None:
            out = data[np.nonzero(data >= min_v)]
            out = out[np.nonzero(out <= max_v)]
            return out
        else:
            return data

    # Set axis tick label formatting
    def __set_label_format(self, label_format, axis):
        if label_format is None:
            return

        if label_format == "default":
            axis.set_major_formatter(self.axis_formatter)
        else:
            tick_formatter = ticker.StrMethodFormatter(label_format)
            axis.set_major_formatter(tick_formatter)
            axis.set_minor_formatter(tick_formatter)

    """ Create a graph from given input data. """

    def graph(
        self,
        x,
        y,
        y_errors=None,
        title="",
        label="",
        x_axis=axis_options(label="x"),
        y_axis=axis_options(label="y"),
        color=None,
        marker=".",
        lgd_ops=legend_options(),
        figsize=(8, 8),
        layout="constrained",
    ):
        # Create fresh plot
        fig = plt.figure(figsize=figsize, layout=layout)
        ax = fig.add_subplot(1, 1, 1)

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(x_axis.label)
        ax.set_ylabel(y_axis.label)
        ax.grid(True, alpha=0.25)

        # Plot log scale
        if x_axis.log_scale is not None:
            ax.set_xscale("log", base=x_axis.log_scale)
        if y_axis.log_scale is not None:
            ax.set_yscale("log", base=y_axis.log_scale)

        if x_axis.tick_positions is not None:
            ax.set_xticks(x_axis.tick_positions)
            ax.tick_params(axis="x", which="major", pad=7)

        if y_axis.tick_positions is not None:
            ax.set_yticks(y_axis.tick_positions)
            ax.tick_params(axis="y", which="major", pad=7)

        # Restrict x and y ranges
        x = self.__apply_boundary(x, x_axis.min, x_axis.max)
        y = self.__apply_boundary(y, y_axis.min, y_axis.max)

        # Format of tick labels
        self.__set_label_format(x_axis.label_format, ax.xaxis)
        self.__set_label_format(y_axis.label_format, ax.yaxis)

        # Nothing left to do
        if len(x) == 0:
            self.logger.debug(rf" create graph: empty data {label}")
            return plt_data(fig=fig, axes=ax)

        if len(x) != len(y):
            self.logger.debug(rf" create graph: x range does match y range {label}")
            return plt_data(fig=fig, axes=ax, errors=y_errors)

        data = ax.errorbar(
            x=x, y=y, label=label, yerr=y_errors, marker=marker, color=color
        )

        # Add legend
        lgd = self.__add_legend(ax, lgd_ops)

        return plt_data(fig=fig, axes=ax, lgd=lgd, data=data, errors=y_errors)

    """ Add new graph to an existing plot """

    def add_graph(
        self,
        plot,
        x,
        y,
        y_errors=None,
        label="",
        marker="+",
        color=None,
    ):
        # Nothing left to do
        if len(y) == 0 or plot.data is None:
            self.logger.debug(rf" add graph: empty data {label}")
            return plot

        # Add new data to old plot axis
        data = plot.axes.errorbar(
            x=x,
            y=y,
            label=label,
            yerr=y_errors,
            color=color,
            marker=marker,
        )

        self.__update_legend(plot.lgd)

        # Rescale the plot
        plot.axes.relim()
        plot.axes.autoscale_view()

        return plt_data(
            fig=plot.fig, axes=plot.axes, lgd=plot.lgd, data=data, errors=y_errors
        )

    """
    Create a histogram from given input data. The normalization is achieved by
    dividing the bin count by the total number of observations. The error is
    calculated as the square root of the bin content.
    """

    def hist1D(
        self,
        x,
        bins=1,
        errors=None,
        w=None,
        title="",
        label="",
        x_axis=axis_options(label="x"),
        y_axis=axis_options(label=""),
        color=default_color,
        alpha=0.75,
        normalize=False,
        show_error=False,
        show_stats=True,
        u_outlier=-1,
        o_outlier=-1,
        lgd_ops=legend_options(),
        figsize=(8, 8),
        layout="compressed",
    ):

        # Create fresh plot
        fig = plt.figure(figsize=figsize, layout=layout)
        ax = fig.add_subplot(1, 1, 1)

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(x_axis.label)
        ax.set_ylabel(y_axis.label)
        ax.grid(True, alpha=0.25)

        # Plot log scale
        if x_axis.log_scale is not None:
            ax.set_xscale("log", base=x_axis.log_scale)
        if y_axis.log_scale is not None:
            ax.set_yscale("log", base=y_axis.log_scale)

        # Leave x-axis with default formatter for 1D histograms
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        self.__set_label_format(y_axis.label_format, ax.yaxis)

        # Do calculations on data in the range of the histogram
        x_min, x_max = self.__get_axis_boundaries(x, x_axis)
        if x_axis.min is not None and x_axis.max is not None:
            x = self.__apply_boundary(x, x_min, x_max)

        # Display number of entries in under- and overflow bins
        underflow = len(np.argwhere(x < x_min))
        overflow = len(np.argwhere(x > x_max))
        if u_outlier >= 0 or o_outlier >= 0:
            underflow = underflow + u_outlier
            overflow = overflow + o_outlier

        # Nothing left to do
        if len(x) == 0:
            self.logger.debug(rf" create hist: empty data {label}")
            return plt_data(fig=fig, axes=ax)

        # Histogram normalization
        scale = 1.0 / len(x) if normalize else 1.0

        # Format the 'newline'
        newline = "\n"

        # Name of the data collection
        label_str = f"{label} ({len(x)} entries)"
        if u_outlier >= 0 or o_outlier >= 0:
            label_str = (
                label_str
                + f"{newline} underflow: {underflow}"
                + f"{newline} overflow:  {overflow}"
            )

        # Fill data
        data, bins, _ = ax.hist(
            x,
            weights=w,
            range=(x_min, x_max),
            bins=bins,
            label=label_str,
            histtype="stepfilled",
            density=normalize,
            alpha=alpha,
            facecolor=color,
            edgecolor=color,
        )

        # Add some additional information
        if show_stats:
            mean = np.mean(x, axis=0)
            # rms  = np.sqrt(np.mean(np.square(x)))
            stdev = np.std(x, axis=0)

            # Create empty plot with blank marker containing the extra label
            ax.plot(
                [],
                [],
                " ",
                label=rf"data:"
                rf"{newline}mean    = {mean:.2e}"
                rf"{newline}stddev  = {stdev:.2e}",
            )
        else:
            mean = None
            stdev = None

        # Calculate the bin error
        bin_centers = 0.5 * (bins[1:] + bins[:-1])
        err = np.sqrt(scale * data) if errors is None else errors
        if show_error or errors is not None:
            ax.errorbar(
                bin_centers,
                data,
                yerr=err,
                fmt=".",
                linestyle="",
                linewidth=0.4,
                color="black",
                capsize=2.5,
            )

        # Add legend
        lgd = self.__add_legend(ax, lgd_ops)

        # Adjust spacing in box
        lgd.legend_handles[0].set_visible(False)
        if show_stats:
            lgd.legend_handles[1].set_visible(False)
        for vpack in lgd._legend_handle_box.get_children():
            for hpack in vpack.get_children():
                hpack.get_children()[0].set_width(0)

        return plt_data(
            fig=fig,
            axes=ax,
            lgd=lgd,
            data=data,
            bins=bins,
            mu=mean,
            rms=stdev,
            errors=err,
        )

    """ Add new histogram to an existing plot """

    def add_hist(
        self,
        old_hist,
        x,
        errors=None,
        w=None,
        label="",
        color="tab:orange",
        alpha=0.75,
        normalize=False,
        show_error=False,
    ):

        # Do calculations on data in the range of the histogram
        x = self.__apply_boundary(x, np.min(old_hist.bins), np.max(old_hist.bins))

        # Nothing left to do
        if len(x) == 0 or old_hist.data is None:
            self.logger.debug(rf" add hist: empty data {label}")
            return old_hist

        # Add new data to old hist axis
        scale = 1.0 / len(x) if normalize else 1.0
        data, bins, _ = old_hist.axes.hist(
            x=x,
            bins=old_hist.bins,
            label=f"{label} ({len(x)} entries)",
            weights=w,
            histtype="stepfilled",
            facecolor=color,
            alpha=alpha,
            edgecolor=color,
        )

        # Calculate the bin error
        bin_centers = 0.5 * (bins[1:] + bins[:-1])
        err = np.sqrt(scale * data) if errors is None else errors
        if show_error or errors is not None:
            old_hist.axes.errorbar(
                bin_centers,
                data,
                yerr=err,
                fmt=".",
                linestyle="",
                linewidth=0.4,
                color="black",
                capsize=2.5,
            )

        # Update legend
        self.__update_legend(old_hist.lgd)

        return plt_data(
            fig=old_hist.fig,
            axes=old_hist.axes,
            lgd=old_hist.lgd,
            data=data,
            bins=bins,
            errors=err,
        )

    """
    Plot the ratio of two histograms. The data is assumed to be uncorrelated.
    """

    def add_ratio(
        self, nom, denom, label, color="tab:red", set_log=False, show_error=False
    ):

        # Resize figure
        nom.fig.set_figheight(7)
        nom.fig.set_figwidth(8)

        if nom.bins is None or denom.bins is None:
            return plt_data(fig=nom.fig, axes=nom.axes)

        if len(nom.bins) != len(denom.bins):
            return plt_data(fig=nom.fig, axes=nom.axes)

        # Remove ticks/labels that are already visible on the ratio plot
        x_label = nom.axes.xaxis.get_label().get_text()
        nom.axes.tick_params(
            axis="x", which="both", bottom=True, top=False, labelbottom=False
        )
        nom.axes.set_xlabel("")

        # Don't print a warning when dividing by zero
        with np.errstate(divide="ignore"), np.errstate(invalid="ignore"):
            # Filter out nan results from division by zero
            ratio = np.nan_to_num(nom.data / denom.data, nan=0, posinf=0)

            # Calculate errors by Gaussian propagation
            bin_centers = 0.5 * (nom.bins[1:] + nom.bins[:-1])
            n_data, d_data = (nom.data, denom.data)

            # Gaussian approximation for large number of events in bin
            # Note: Should be Clopper-Pearson
            n_err, d_err = (nom.errors, denom.errors)
            errors = np.nan_to_num(
                np.sqrt(
                    np.square(n_err / d_data)
                    + np.square(n_data * d_err / np.square(d_data))
                ),
                nan=0,
                posinf=0,
            )

        # Create new axes on the bottom of the current axes
        # The first argument of the new_vertical(new_horizontal) method is
        # the height (width) of the axes to be created in inches.
        divider = make_axes_locatable(nom.axes)
        ratio_plot = divider.append_axes("bottom", 1.2, pad=0.2, sharex=nom.axes)
        # Ratio should be around 1: Don't use scientific notation/offset
        ratio_plot.axes.yaxis.set_major_formatter(
            ticker.ScalarFormatter(useOffset=False)
        )
        if show_error:
            ratio_plot.errorbar(
                bin_centers, ratio, yerr=errors, label=label, color=color, fmt="."
            )
        else:
            ratio_plot.plot(
                bin_centers,
                ratio,
                label=label,
                color=color,
                marker=".",
                linestyle="",
            )

        # Refine plot
        ratio_plot.set_xlabel(x_label)
        ratio_plot.set_ylabel("ratio")
        ratio_plot.grid(True, alpha=0.25)

        # Plot log scale
        if set_log:
            ratio_plot.set_yscale("log")

        # Add a horizontal blue line at y = 1.
        ratio_plot.axline((nom.bins[0], 1), (nom.bins[-1], 1), linewidth=1, color="b")

        nom.fig.set_size_inches((9, 9))

        return plt_data(fig=nom.fig, axes=ratio_plot, errors=errors)

    """
    Create a 2D histogram from given input data. If z values are given they will
    be used as weights per bin.
    """

    def hist2D(
        self,
        x,
        y,
        z=None,
        x_bins=1,
        y_bins=1,
        x_axis=axis_options(label="x"),
        y_axis=axis_options(label="y"),
        z_axis=axis_options(label=""),
        title="",
        label="",
        color=default_color,
        alpha=0.75,
        show_stats=True,
        figsize=(8, 6),
    ):

        # Create fresh plot
        fig = plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(1, 1, 1)

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(x_axis.label)
        ax.set_ylabel(y_axis.label)

        # Do calculations on data in the range of the histogram
        x_min, x_max = self.__get_axis_boundaries(x, x_axis)
        if x_axis.min is not None and x_axis.max is not None:
            x = self.__apply_boundary(x, x_min, x_max)

        y_min, y_max = self.__get_axis_boundaries(y, y_axis)
        if y_axis.min is not None and y_axis.max is not None:
            y = self.__apply_boundary(y, y_min, y_max)

        # Nothing left to do
        if len(x) == 0 or len(y) == 0:
            self.logger.debug(rf" create hist: empty data {label}")
            return plt_data(fig=fig, axes=ax)

        # Fill data
        data, _, _, hist = ax.hist2d(
            x,
            y,
            weights=z,
            range=[(x_min, x_max), (y_min, y_max)],
            bins=(x_bins, y_bins),
            label=f"{label}  ({len(x)*len(y)} entries)",
            facecolor=mcolors.to_rgba(color, alpha),
            edgecolor=None,
            rasterized=True,
        )

        # Add some additional information
        if show_stats:
            x_mean = np.mean(x, axis=0)
            x_rms = np.sqrt(np.mean(np.square(x)))
            y_mean = np.mean(y, axis=0)
            y_rms = np.sqrt(np.mean(np.square(y)))

            # Create empty plot with blank marker containing the extra label
            newline = "\n"
            ax.plot(
                [],
                [],
                " ",
                label=rf"xMean = {x_mean:.2e}"
                rf"{newline}xRMS  = {x_rms:.2e}"
                rf"yMean = {y_mean:.2e}"
                rf"{newline}yRMS  = {y_rms:.2e}",
            )

        # Add the colorbar
        fig.colorbar(hist, label=z_axis.label)

        return plt_data(fig=fig, axes=ax, data=data)

    """ Create a 2D scatter plot """

    def scatter(
        self,
        x,
        y,
        x_axis=axis_options(label=""),
        y_axis=axis_options(label=""),
        title="",
        label="",
        color=default_color,
        alpha=1,
        show_stats=lambda x, _: f"{len(x)} entries",
        lgd_ops=legend_options(),
        figsize=(8, 6),
    ):

        fig = plt.figure(figsize=figsize, layout="constrained")
        ax = fig.add_subplot(1, 1, 1)

        # Refine plot
        ax.set_title(title)
        ax.set_xlabel(x_axis.label)
        ax.set_ylabel(y_axis.label)
        ax.grid(True, alpha=0.25)

        # Create empty plot with blank marker containing the extra label
        ax.plot([], [], " ", label=show_stats(x, y))
        scatter = ax.scatter(
            x, y, label=label, c=color, s=0.1, alpha=alpha, rasterized=True
        )

        # Add legend
        lgd = self.__add_legend(ax, lgd_ops)

        # Refine legend
        self.__adjust_lgd_label_spacing(lgd)

        return plt_data(fig=fig, axes=ax, lgd=lgd, data=scatter)

    """ Add new data in a different color to a scatter plot """

    def highlight_region(self, plot_data, x, y, color, label=""):

        if label == "":
            plot_data.axes.scatter(x, y, c=color, alpha=1, s=0.1, rasterized=True)
        else:
            plot_data.axes.scatter(
                x, y, c=color, alpha=1, s=0.1, label=label, rasterized=True
            )

            # Update legend
            self.__update_legend(plot_data.lgd)

        # Refine legend
        self.__adjust_lgd_label_spacing(plot_data.lgd)

    """ Fit a Gaussian to a 1D distribution and plot in the same figure. """

    def fit_gaussian(self, dist, color="tab:orange"):

        # Calculate bin centers from bin edges
        bins = dist.bins
        if bins is None:
            # If fit failed, return empty result
            return None, None

        bin_centers = [(b1 + b2) / 2 for b1, b2 in zip(bins, bins[1:])]

        # Gaussian distribution with all fit parameters
        def __gaussian(x, a, mean, sigma):
            return (
                a
                / (math.sqrt(2 * math.pi) * sigma)
                * np.exp(-((x - mean) ** 2 / (2 * sigma**2)))
            )

        # Gaussian fit
        try:
            from scipy.optimize import curve_fit
        except ImportError:
            print("WARNING: Could not find scipy: Skipping fit")
        else:
            try:
                # Initial estimators
                mean = np.mean(bin_centers, axis=0)
                sigma = np.std(bin_centers, axis=0)
                a = np.max(dist.data) * (math.sqrt(2 * math.pi) * sigma)

                popt, _ = curve_fit(
                    __gaussian, bin_centers, dist.data, p0=[a, mean, sigma]
                )
            except RuntimeError:
                # If fit failed, return empty result
                return None, None

            # If the fitting was successful, plot the curve
            mu = float(f"{popt[1]:.2e}")  # < formatting the sig. digits
            sig = float(f"{popt[2]:.2e}")
            newline = "\n"
            plot_label = (
                rf"gaussian fit:{newline}$\mu$ = {mu:.2e}"
                + rf"{newline}$\sigma$ = {abs(sig):.2e}"
            )

            # Generate points for the curve
            min_val = min(bin_centers)
            max_val = max(bin_centers)
            step = (max_val - min_val) / 1000
            x = [v for v in np.arange(min_val, max_val + step, step)]

            dist.axes.plot(
                x,
                __gaussian(x, *popt),
                label=plot_label,
                color=color,
            )

            # Update legend
            self.__update_legend(dist.lgd)

            # Adjust spacing in box
            dist.lgd.legend_handles[0].set_visible(False)
            for vpack in dist.lgd._legend_handle_box.get_children()[:-1]:
                for hpack in vpack.get_children():
                    hpack.get_children()[0].set_width(0)

            return popt[1], abs(popt[2])

        return None, None

    """ Draw a vertical line in a given plot"""

    def vertical_line(self, plot_data, x, y=None, color="b", label=""):

        plot_data.axes.axvline(x=x, color=color, linestyle="--")

        ymin, ymax = plot_data.axes.get_ylim()
        plot_data.axes.text(
            x,
            ymin + (ymax - ymin) / 2 if y is None else y,
            label,
            ha="center",
            va="center",
            backgroundcolor="white",
        )

    """ Write a plot to disk """

    def write_plot(
        self, plot_data, name="plot", file_format="svg", out_prefix="", dpi=450
    ):
        if out_prefix == "":
            file_name = self.output_prefix + "/" + name + "." + file_format
        else:
            file_name = out_prefix + name + "." + file_format

        plot_data.fig.savefig(file_name, dpi=dpi)
        plt.close(plot_data.fig)

    """ Write a plot as svg """

    def write_svg(self, plot_data, name, out_prefix=""):
        self.write_plot(plot_data, name, ".svg", out_prefix)

    """ Write a plot as pdf """

    def write_pdf(self, plot_data, name, out_prefix=""):
        self.write_plot(plot_data, name, ".pdf", out_prefix)

    """ Write a plot as png """

    def write_png(self, plot_data, name, out_prefix=""):
        self.write_plot(plot_data, name, ".png", out_prefix)
