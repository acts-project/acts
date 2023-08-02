// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Result.hpp"

#include <unordered_map>

namespace Acts {

/// @class AdaptiveGridTrackDensity
/// @brief Implements a 1-dim density grid to be filled with
/// track Gaussian distributions. Each single track is modelled
/// as a 2-dim Gaussian distribution grid in the d0-z0 plane,
/// but only the overlap with the z-axis (i.e. a 1-dim density
/// vector) needs to be calculated.
/// The position of the highest track density (of either a single
/// bin or the sum of a certain region) can be determined.
/// Single tracks can be cached and removed from the overall density.
/// Unlike the GaussianGridTrackDensity, the overall density vector
/// grows adaptively with the tracks densities being added to the grid.
///
/// @tparam trkGridSize The 2-dim grid size of a single track, i.e.
/// a single track is modelled as a (trkGridSize x trkGridSize) grid
/// in the d0-z0 plane. Note: trkGridSize has to be an odd value.
template <int trkGridSize = 15>
class AdaptiveGridTrackDensity {
  // Assert odd trkGridSize
  static_assert(trkGridSize % 2);

 public:
  using DensityMap = std::unordered_map<int, float>;

  /// The configuration struct
  struct Config {
    /// @param binSize_ The binSize in mm
    Config(float binSize_ = 0.1) : binSize(binSize_) {}

    // Z size of one single bin in grid
    float binSize;  // mm

    // Do NOT use just the z-bin with the highest
    // track density, but instead check the (up to)
    // first three density maxima (only those that have
    // a maximum relative deviation of 'relativeDensityDev'
    // from the main maximum) and take the z-bin of the
    // maximum with the highest surrounding density sum
    bool useHighestSumZPosition = false;

    // The maximum relative density deviation from the main
    // maximum to consider the second and third maximum for
    // the highest-sum approach from above
    float maxRelativeDensityDev = 0.01;
  };

  AdaptiveGridTrackDensity(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Finds the maximum density of a DensityMap
  /// @param densityMap Map between z bins and corresponding density value
  /// @return Iterator of the map entry with the highest density value
  DensityMap::const_iterator highestDensityEntry(
      const DensityMap& densityMap) const;

  /// @brief Returns the z position of maximum (surrounding) track density
  ///
  /// @param densityMap Map from z bins to corresponding track density
  ///
  /// @return The z position of maximum track density
  Result<float> getMaxZPosition(DensityMap& densityMap) const;

  /// @brief Returns the z position of maximum track density and
  /// the estimated width
  ///
  /// @param densityMap Map from z bins to corresponding track density
  ///
  /// @return The z position of maximum track density and width
  Result<std::pair<float, float>> getMaxZPositionAndWidth(
      DensityMap& densityMap) const;

  /// @brief Adds a single track to the overall grid density
  ///
  /// @param trk The track to be added.
  /// @param mainDensityMap Map from z bins to corresponding track density.
  ///
  /// @return The density map of the track that was added
  DensityMap addTrack(const BoundTrackParameters& trk,
                      DensityMap& mainDensityMap) const;

  /// @brief Removes a track from the overall grid density
  ///
  /// @param trackDensityMap Map from z bins to corresponding track density. The track density comes from a single track.
  /// @param mainDensityMap Map from z bins to corresponding track density. The track density comes an arbitrary number of tracks.
  void subtractTrack(const DensityMap& trackDensityMap,
                     DensityMap& mainDensityMap) const;

 private:
  /// @brief Function that creates a track density map, i.e., a map of z bins to corresponding density values coming from a single track.
  ///
  /// @param offset Offset in d0 direction, to account for the 2-dim part
  /// of the Gaussian track distribution
  /// @param cov The track covariance matrix
  /// @param distCtrD The distance in d0 from the track position to its
  /// bin center in the 2-dim grid
  /// @param centralZBin Central z bin of the track (where its density is the highest)
  /// @param distCtrZ The distance in z0 from the track position to its
  /// bin center in the 2-dim grid
  DensityMap createTrackGrid(int offset, const SymMatrix2& cov, float distCtrD,
                             int centralZBin, float distCtrZ) const;

  /// @brief Function that estimates the seed width based on the full width
  /// at half maximum (FWHM) of the maximum density peak
  ///
  /// @param densityMap Map from z bins to corresponding track density
  /// @param maxZ z-position of the maximum density value
  ///
  /// @return The width
  Result<float> estimateSeedWidth(const DensityMap& densityMap,
                                  float maxZ) const;

  /// @brief Helper to retrieve values according to a 2-dim normal distribution
  float normal2D(float d, float z, const SymMatrix2& cov) const;

  /// @brief Checks the (up to) first three density maxima (only those that have
  /// a maximum relative deviation of 'relativeDensityDev' from the main
  /// maximum) and take the z-bin of the maximum with the highest surrounding
  /// density
  ///
  /// @param densityMap Map between z bins and corresponding density value
  ///
  /// @return The z-bin position
  int highestDensitySumBin(DensityMap& densityMap) const;

  /// @brief Calculates the density sum of a z-bin and its two neighboring bins
  /// as needed for 'highestDensitySumBin'
  ///
  /// @param densityMap Map between z bins and corresponding density value
  /// @param zBin The center z-bin whose neighbors we want to sum up
  ///
  /// @return The sum
  float getDensitySum(const DensityMap& densityMap, int zBin) const;

  Config m_cfg;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveGridTrackDensity.ipp"
