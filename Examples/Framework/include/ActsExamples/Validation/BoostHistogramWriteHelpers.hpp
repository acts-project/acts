// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsExamples {

class BoostHistogram1D;
class BoostHistogram2D;
class BoostProfileHistogram;
class BoostEfficiency1D;
class BoostEfficiency2D;

/// Helper functions to write boost histograms to ROOT files
///
/// These functions encapsulate the conversion from boost::histogram to ROOT
/// and handle the ROOT file I/O. They convert the boost histogram to a ROOT
/// histogram using BoostHistogramToRoot::toRoot(), write it to the current
/// ROOT directory, and clean up the temporary ROOT object.
namespace BoostHistogramWriteHelpers {

/// Write BoostHistogram1D to current ROOT directory
///
/// @param hist The boost histogram to write
void write(const BoostHistogram1D& hist);

/// Write BoostHistogram2D to current ROOT directory
///
/// @param hist The boost histogram to write
void write(const BoostHistogram2D& hist);

/// Write BoostProfileHistogram to current ROOT directory
///
/// @param hist The boost profile histogram to write
void write(const BoostProfileHistogram& hist);

/// Write BoostEfficiency1D to current ROOT directory
///
/// @param hist The boost efficiency histogram to write
void write(const BoostEfficiency1D& hist);

/// Write BoostEfficiency2D to current ROOT directory
///
/// @param hist The boost efficiency histogram to write
void write(const BoostEfficiency2D& hist);

}  // namespace BoostHistogramWriteHelpers
}  // namespace ActsExamples
