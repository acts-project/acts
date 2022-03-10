// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

// TODO top-level description
// TODO finish descriptions here
// TODO better name

// namespace Acts {
// namespace HK {

enum Label { None = 0 };

// Simple wrapper around a generic cell type
// The Cell type should have the following accessor functions defined:
//   int get_cell_row(const Cell&)
//   int get_cell_column(const Cell&)
//   float get_cell_activation(const Cell&)
template <typename Cell>
struct LabeledCell {
    Cell const * ptr;
    mutable Label lbl;
    explicit LabeledCell(const Cell& cell); : ptr{&cell}, lbl{Label::None} {}
};

template <typename Cell>
using LabeledCellCollection = std::vector<LabeledCell<Cell>>;

/// @brief createClusters
/// TODO: descr
/// @param [in] begin input iterator to beginning of cell sequence
/// @param [in] end input iterator to end of cell sequence
/// @param [out] out output iterator to cluster sequence (e.g. std::back_inserter)
/// @param [in] commonCorner flag indicating if cells sharing a common corner should be merged into one cluster
/// @param [in] threshold possible activation threshold to be applied
/// @return nothing (clusters are inserted in `out')
template <typename Cell, typename ClusterT, typename InputIt, typename OutputIt>
void createClusters(InputIt begin, InputIt end, OutputIt out, bool commonCorner = true, float threshold = 0.);

/// @brief labelClusters
/// TODO: descr
/// @param [in] begin input interator to beginning of cell sequence
/// @param [in] end input iterator to end of cell sequence
/// @param [in] commonCorner flag indicating if cells sharing a common corner should be merged into one cluster
/// @param [in] threshold possible activation threshold to be applied
/// @return collection of cells wrapped with connected component labels. May not be in same order as the input collection.
template <typename Cell, typename InputIt>
LabeledCellCollection<Cell> labelClusters(InputIt begin, InputIt end, bool commonCorner = true, float threshold = 0.);

// } // HK namespace
// } // Acts namespace

#include "Acts/Clusterization/NewClusterization.ipp"
