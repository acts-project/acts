/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <cassert>

namespace traccc::details {

TRACCC_HOST_DEVICE inline unsigned int find_root(
    const vecmem::device_vector<unsigned int>& labels, unsigned int e) {

    unsigned int r = e;
    assert(r < labels.size());
    while (labels[r] != r) {
        r = labels[r];
        assert(r < labels.size());
    }
    return r;
}

TRACCC_HOST_DEVICE inline unsigned int make_union(
    vecmem::device_vector<unsigned int>& labels, unsigned int e1,
    unsigned int e2) {

    unsigned int e;
    if (e1 < e2) {
        e = e1;
        assert(e2 < labels.size());
        labels[e2] = e;
    } else {
        e = e2;
        assert(e1 < labels.size());
        labels[e1] = e;
    }
    return e;
}

template <typename T1, typename T2>
TRACCC_HOST_DEVICE inline bool is_adjacent(const edm::silicon_cell<T1>& a,
                                           const edm::silicon_cell<T2>& b) {

    return (a.channel0() - b.channel0()) * (a.channel0() - b.channel0()) <= 1 &&
           (a.channel1() - b.channel1()) * (a.channel1() - b.channel1()) <= 1 &&
           a.module_index() == b.module_index();
}

template <typename T1, typename T2>
TRACCC_HOST_DEVICE inline bool is_far_enough(const edm::silicon_cell<T1>& a,
                                             const edm::silicon_cell<T2>& b) {

    assert((a.channel1() >= b.channel1()) ||
           (a.module_index() != b.module_index()));
    return (a.channel1() > (b.channel1() + 1)) ||
           (a.module_index() != b.module_index());
}

TRACCC_HOST_DEVICE inline unsigned int sparse_ccl(
    const edm::silicon_cell_collection::const_device& cells,
    vecmem::device_vector<unsigned int>& labels) {

    unsigned int nlabels = 0;

    // The number of cells.
    const unsigned int n_cells = cells.size();

    // first scan: pixel association
    unsigned int start_j = 0;
    for (unsigned int i = 0; i < n_cells; ++i) {

        labels[i] = i;
        unsigned int ai = i;
        for (unsigned int j = start_j; j < i; ++j) {

            if (is_adjacent(cells.at(i), cells.at(j))) {
                ai = make_union(labels, ai, find_root(labels, j));
            } else if (is_far_enough(cells.at(i), cells.at(j))) {
                ++start_j;
            }
        }
    }

    // second scan: transitive closure
    for (unsigned int i = 0; i < n_cells; ++i) {
        if (labels[i] == i) {
            labels[i] = nlabels++;
        } else {
            labels[i] = labels[labels[i]];
        }
    }
    return nlabels;
}

}  // namespace traccc::details
