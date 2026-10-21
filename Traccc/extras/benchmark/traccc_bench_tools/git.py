# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0


def is_parent_of(subj, parent_str):
    for p in subj.iter_parents():
        if str(subj) == parent_str or str(p) == parent_str:
            return True
    return False
