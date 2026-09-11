# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0


def simplify_name(name):
    """
    Simplify a kernel name in a NCU profile, stripping the return type,
    arguments, and templates.
    """
    if name[:5] == "void ":
        name = name[5:]

    val = ""

    while name:
        if name[:2] == "::":
            val = ""
            name = name[2:]
        elif name[0] == "(" or name[0] == "<":
            return val
        else:
            val = val + name[0]
            name = name[1:]

    raise RuntimeError("An error occured in name simpliciation")


def map_name(name):
    """
    Map some non-traccc kernel names to human-readable versions.
    """
    if name in [
        "DeviceRadixSortUpsweepKernel",
        "RadixSortScanBinsKernel",
        "DeviceRadixSortDownsweepKernel",
        "DeviceRadixSortSingleTileKernel",
        "DeviceMergeSortBlockSortKernel",
        "DeviceMergeSortMergeKernel",
        "DeviceMergeSortPartitionKernel",
        "DeviceRadixSortHistogramKernel",
        "DeviceRadixSortOnesweepKernel",
        "DeviceRadixSortExclusiveSumKernel",
    ]:
        return "Thrust::sort"
    elif name in ["_kernel_agent", "static_kernel"]:
        return "unknown"
    else:
        return name
