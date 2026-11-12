// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>

namespace ActsExamples {

/// The projection of a cluster onto one axis of the segmentation.
///
/// The charges are calibrated charges, in the unit of the readout. A raw
/// quantity, for example the time over threshold, must be converted before.
struct ClusterAxisProjection {
  /// The edges of the first bin of the cluster on this axis
  double firstLower = 0;
  double firstUpper = 0;
  /// The edges of the last bin of the cluster on this axis
  double lastLower = 0;
  double lastUpper = 0;
  /// The summed charge of the fired cells in the first and the last bin
  double firstCharge = 0;
  double lastCharge = 0;
  /// The number of fired cells in the first and the last bin
  std::size_t firstCells = 1;
  std::size_t lastCells = 1;
  /// The widths of the bins before the first and after the last bin. Zero if
  /// the sensor ends there.
  double beforeWidth = 0;
  double afterWidth = 0;
};

/// The footprint of the predicted track in the sensor, along one axis.
///
/// The footprint comes from the predicted direction as it is. The error of
/// the direction is not used: a track fit can inflate it far above the
/// actual error, and the actual error is small after the first measurements.
struct ClusterAxisFootprint {
  /// The extent of the path through the sensor on this axis
  double length = 0;
  /// The length of the full path through the sensor
  double pathLength = 0;
};

/// The response of one cell.
struct ClusterAxisReadout {
  /// A cell fires if its charge is above this threshold
  double threshold = 0;
  /// The standard deviation of the additive Gaussian noise of one cell
  double noise = 0;
  /// The readout reports the charge of a fired cell, not only the hit
  bool analogue = false;
  /// The expected charge per path length, in the unit of the charges. The
  /// geometric digitization uses the path length as the charge, so it is 1.
  double chargePerLength = 1;
};

/// The mean and the variance of the track position on one axis.
struct ClusterAxisPosterior {
  double mean = 0;
  double variance = 0;
};

/// The log likelihood of the projected cluster for the track midpoint `mid`.
///
/// The track deposits its charge uniformly along a straight segment with the
/// length of the footprint, centred on `mid`. Each cell gets additive Gaussian
/// noise. The first and the last bin fire, and for analogue readout their
/// charge is Gaussian around the deposited charge. A touched neighbour bin
/// stays below the threshold. An untouched cell never fires.
///
/// @param projection The cluster on this axis
/// @param footprint The footprint of the track
/// @param readout The response of one cell
/// @param mid The midpoint of the track segment on this axis
///
/// @return The log likelihood up to a constant, which is zero for a perfect
///         match, or minus infinity if the segment does not touch both edge
///         bins
double clusterAxisLogLikelihood(const ClusterAxisProjection& projection,
                                const ClusterAxisFootprint& footprint,
                                const ClusterAxisReadout& readout, double mid);

/// Calibrate the track position on one axis from the projected cluster.
///
/// Computes the mean and the variance of the posterior of the midpoint, with
/// `clusterAxisLogLikelihood` and a uniform prior, on a grid. The grid is
/// refined to the support of the posterior.
///
/// @param projection The cluster on this axis
/// @param footprint The footprint of the track
/// @param readout The response of one cell
///
/// @return The posterior moments, or uniform over the cluster if the cluster
///         is incompatible with the footprint, which is mostly a merged
///         cluster
ClusterAxisPosterior calibrateClusterAxis(
    const ClusterAxisProjection& projection,
    const ClusterAxisFootprint& footprint, const ClusterAxisReadout& readout);

}  // namespace ActsExamples
