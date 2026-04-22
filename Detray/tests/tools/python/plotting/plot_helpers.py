# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from collections import namedtuple
import numpy as np

# ------------------------------------------------------------------------------
# Common helper types to configure plots
# ------------------------------------------------------------------------------

""" Pass plotting data between functions """
plt_data = namedtuple(
    "plt_data",
    "fig axes lgd data bins mu rms errors",
    defaults=[None, None, None, None, 0, -1, -1, None],
)

""" Configuration for plot axes """
axis_options = namedtuple(
    "axis_options",
    "label min max log_scale tick_positions label_format",
    defaults=["x", None, None, None, None, None],
)

""" Configuration for plot legends """
legend_options = namedtuple(
    "legend_options",
    "title loc ncol colspacing handletextpad horiz_anchor vert_anchor",
    defaults=[None, "upper right", 1, 1, 1, 0.5, 0.5],
)

# ------------------------------------------------------------------------------
# Common helpers for plotting measurement data
# ------------------------------------------------------------------------------

""" Filter the data in a data frame by a given prescription """


def filter_data(data, filter=lambda df: [], variables=[]):
    data_coll = []

    # Get global data
    if len(filter(data)) == 0:
        for var in variables:
            data_coll.append(data[var].to_numpy(dtype=np.double))

    # Filtered data
    else:
        filtered = data.loc[filter]
        for var in variables:
            data_coll.append(filtered[var].to_numpy(dtype=np.double))

    return data_coll[0] if len(data_coll) == 1 else tuple(data_coll)
