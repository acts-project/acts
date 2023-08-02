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
 // TODO: go back to unordered map
  using DensityMap = std::map<int, float>;
  using TrackGridVector = Eigen::Matrix<float, trkGridSize, 1>;

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

  /// @brief TODO
  /// @param densityMap Map between z bins and corresponding density value
  /// @return 
  DensityMap::const_iterator highestDensityEntry(const DensityMap& densityMap) const;

  /// @brief Returns the z position of maximum (surrounding) track density
  ///
  /// @param mainGridDensity The main 1-dim density grid along the z-axis
  /// @param mainGridZValues The corresponding z-bin values of the track
  /// densities along the z-axis
  ///
  /// @return The z position of maximum track density
  Result<float> getMaxZPosition(std::vector<float>& mainGridDensity,
                                const std::vector<int>& mainGridZValues) const;

  Result<float> getMaxZPosition(DensityMap& densityMap) const;

  /// @brief Returns the z position of maximum track density and
  /// the estimated width
  ///
  /// @param mainGridDensity The main 1-dim density grid along the z-axis
  /// @param mainGridZValues The corresponding z-bin values of the track
  /// densities along the z-axis
  ///
  /// @return The z position of maximum track density and width
  Result<std::pair<float, float>> getMaxZPositionAndWidth(
      std::vector<float>& mainGridDensity,
      const std::vector<int>& mainGridZValues) const;

  Result<std::pair<float, float>> getMaxZPositionAndWidth(
      DensityMap& densityMap) const;

  /// @brief Adds a single track to the overall grid density
  ///
  /// @param trk The track to be added
  /// @param mainGridDensity The main 1-dim density grid along the z-axis
  /// @param mainGridZValues The corresponding z-bin values of the track
  /// densities along the z-axis
  ///
  /// @return A pair storing information about the z-bin position
  /// the track was added (int) and the 1-dim density contribution
  /// of the track itself
  std::pair<int, TrackGridVector> addTrack(
      const BoundTrackParameters& trk, std::vector<float>& mainGridDensity,
      std::vector<int>& mainGridZValues) const;

  DensityMap addTrack(
      const BoundTrackParameters& trk, DensityMap& mainDensityMap) const;

  /// @brief Removes a track from the overall grid density
  ///
  /// @param zBin The center z-bin position the track needs to be
  /// removed from
  /// @param trkGrid The 1-dim density contribution of the track
  /// @param mainGridDensity The main 1-dim density grid along the z-axis
  /// @param mainGridZValues The corresponding z-bin values of the track
  /// densities along the z-axis
  void removeTrackGridFromMainGrid(
      int zBin, const TrackGridVector& trkGrid,
      std::vector<float>& mainGridDensity,
      const std::vector<int>& mainGridZValues) const;

  void subtractTrack(const DensityMap& trackDensityMap, DensityMap& mainDensityMap) const;

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
  DensityMap createTrackGrid(int offset, const SymMatrix2& cov,
                                  float distCtrD, int centralZBin, float distCtrZ) const;

  /// @brief Function that estimates the seed width based on the full width
  /// at half maximum (FWHM) of the maximum density peak
  ///
  /// @param mainGridDensity The main 1-dim density grid along the z-axis
  /// @param mainGridZValues The corresponding z-bin values of the track
  ///                        densities along the z-axis
  /// @param maxZ z-position of the maximum density value
  ///
  /// @return The width
  Result<float> estimateSeedWidth(const std::vector<float>& mainGridDensity,
                                  const std::vector<int>& mainGridZValues,
                                  float maxZ) const;

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
