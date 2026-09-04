# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import collections

GpuSpec = collections.namedtuple("GpuSpec", ["n_sm", "n_threads_per_sm"])
