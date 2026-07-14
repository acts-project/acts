# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0


def harmonic_sum(vals):
    """
    Return the harmonic addition of values, i.e. the reciprocal of the sum of
    the reciprocals.
    """
    return 1.0 / sum(1.0 / x for x in vals)
