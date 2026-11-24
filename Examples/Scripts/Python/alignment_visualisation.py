import numpy as np
import matplotlib.pyplot as plt
from enum import IntEnum
import os
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


def plotAlignImpacts(
    labels,
    centralValues,
    sigmas,
    formatDicts,
    fname="align_impacts",
    title=None,
    xLim=None,
):
    """
    "Impacts plot" style comparison of the alignment corrections.

    Each row of the Y axis is a parameter (layer, coordinate).
    The X axis is the correction in mm,
    centred at zero.

    Parameters
    ----------
    labels: 1D array of parameter labels - defines the number nLabels of vertical entries
    centralValues: 2D array of shape (nSeries, nLabels) - central parameter values for each series
    sigmas: 2D array of shape (nSeries, nLabels) - parameter uncertainty values for each series
    formatDicts: 1D list of dicts, length nSeries - format and label settings for each series
    fname: file name to save the plot
    title: a title to plot
    xLim: an x-axis range limit

    Sign convention: the correction is drawn exactly as the fits return it.
    If truth is the INJECTED displacement, the expected correction is its
    opposite; pass -truth if you want to superimpose them.
    """
    nParams = labels.shape[0]
    nSeries = len(formatDicts)
    assert centralValues.shape[0] == nSeries
    assert centralValues.shape[1] == nParams
    assert sigmas.shape[1] == nParams
    assert sigmas.shape[1] == nParams

    # ---- figure: height proportional to the number of rows
    height = max(4.5, 0.34 * nParams + 1.4)
    fig, ax = plt.subplots(figsize=(7.5, height))

    # ---- zero line
    yPos = np.asarray(range(nParams), int) + 1
    ax.axvline(0.0, color="k", lw=1.1, zorder=2)
    offsetY = 0.3 / nSeries
    offsetSteps = nSeries - 1
    offsetStep = offsetSteps / 2
    for thing in range(nSeries):
        yPositions = yPos - offsetStep * offsetY
        means = centralValues[thing, :]
        errors = sigmas[thing, :]
        format = formatDicts[thing]
        hasErrs = np.count_nonzero(errors)
        if not hasErrs:
            ax.plot(means, yPositions, **format)
        else:
            ax.errorbar(means, yPositions, xerr=errors, **format)
        offsetStep -= 1

    # ---- axes
    ax.set_yticks(yPos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_ylim(yPos[-1] + 0.6, yPos[0] - 0.6)  # row 0 at the top
    ax.set_xlabel("alignment correction [mm]")
    ax.grid(axis="x", ls=":", color="0.8", zorder=0)
    ax.set_axisbelow(True)

    if xLim is not None:
        ax.set_xlim(*xLim)

    if title:
        ax.set_title(title, fontsize=11)

    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(fname + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(fname + ".pdf", bbox_inches="tight")
    return fig, ax


if __name__ == "__main__":

    labels = np.array([1, 2, 3, 4], int)
    centralValues = np.zeros(dtype=float, shape=(3, len(labels)))
    sigmas = np.zeros(dtype=float, shape=(3, len(labels)))
    centralValues[0] = [-0.11, 0.1, 0.0, 0.09]
    centralValues[1] = [-0.12, 0.08, 0.0, 0.11]
    centralValues[2] = [-0.1, 0.1, 0.0, 0.1]
    sigmas[0] = [0.015, 0.015, 0.015, 0.015]
    sigmas[1] = [0.015, 0.015, 0.015, 0.015]
    formatDicts = [
        {"label": "Series 1", "fmt": "o", "color": "k", "markersize": 5},
        {"label": "Series 2", "fmt": "o", "color": "blue", "markersize": 5},
        {
            "label": "Truth",
            "marker": "s",
            "markeredgecolor": "#f31d1d",
            "markerfacecolor": "none",
            "linestyle": "none",
            "markersize": 11,
            "markeredgewidth": 1.6,
        },
    ]

    plotAlignImpacts(
        labels,
        centralValues,
        sigmas,
        formatDicts,
        fname="testImpacts",
        title=None,
        xLim=None,
    )
