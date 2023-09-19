// This file is part of the Acts project.
//
// Copyright (C) 2020-2023 CERN for the benefit of the Acts project
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
/// TODO update comment
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
/// TODO once we switch to c++ 20 we should use a single template
/// parameter of type std::pair<int, int>
template <int spatialTrkGridSize = 15, int temporalTrkGridSize = 1>
class AdaptiveGridTrackDensity {
  // Assert odd spatial and temporal track grid size
  static_assert(spatialTrkGridSize % 2);
  static_assert(temporalTrkGridSize % 2);

 public:
  using Bins = int;  // std::pair<int, int>;
  using DensityMap = std::unordered_map<Bins, float>;

  /// The configuration struct
  struct Config {
    /// @param spatialBinExtent_ The spatial extent of a bin in mm
    // TODO remove default value?
    Config(float spatialBinExtent_ = 0.1)
        : spatialBinExtent(spatialBinExtent_) {
      if (temporalTrkGridSize > 1) {
        throw std::invalid_argument(
            "temporalBinExtent must be provided if temporalTrkGridSize > 1 "
            "(i.e., if time vertex seeding is enabled).");
      }
    }

    /// @param spatialBinExtent_ The spatial extent of a bin in mm
    /// @param temporalBinExtent_ The temporal extent of a bin in TODO: unit
    Config(float spatialBinExtent_, float temporalBinExtent_)
        : spatialBinExtent(spatialBinExtent_),
          temporalBinExtent(temporalBinExtent_) {}

    // Spatial extent of a bin in d0 and z0 direction
    float spatialBinExtent;  // mm

    // Temporal extent of a bin
    std::optional<float> temporalBinExtent = std::nullopt;  // TODO: unit?

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

  /// @brief Calculates the bin center from the bin number
  /// @param bin Bin number
  /// @param binExtent Bin extent
  /// @return Bin center
  float getBinCenter(int bin, float binExtent) const;

  /// @brief Calculates the bin number corresponding to a d, z, or time value
  /// @param value d, z, or time value
  /// @param binExtent Bin extent
  /// @return Bin number
  int getBin(float value, float binExtent) const;

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
  /// @param trackDensityMap Map from z bins to corresponding track density.
  /// @note The track density comes from a single track.
  /// @param mainDensityMap Map from z bins to corresponding track density.
  /// @note The track density comes an arbitrary number of tracks.
  void subtractTrack(const DensityMap& trackDensityMap,
                     DensityMap& mainDensityMap) const;

 private:
  /// @brief Function that creates a track density map, i.e., a map of z bins
  /// to corresponding density values coming from a single track.
  ///
  /// @param d0 Transverse impact parameter
  /// @param z0 Longitudinal impact parameter
  /// @param centralZBin Central z bin of the track (where its density is the highest)
  /// @param cov 2x2 impact parameter covariance matrix
  DensityMap createTrackGrid(ActsScalar d0, ActsScalar z0, int centralZBin,
                             const Acts::SquareMatrix2& cov) const;

  /// @brief Function that estimates the seed width based on the full width
  /// at half maximum (FWHM) of the maximum density peak
  /// @note This only works if the maximum is sufficiently isolated since
  /// overlapping neighboring peaks might lead to an overestimation of the
  /// seed width.
  ///
  /// @param densityMap Map from z bins to corresponding track density
  /// @param maxZ z-position of the maximum density value
  ///
  /// @return The width
  Result<float> estimateSeedWidth(const DensityMap& densityMap,
                                  float maxZ) const;

  /// @brief Helper to retrieve values according of a nDim-dimensional normal distribution
  /// @note The constant prefactor (2 * pi)^(- nDim / 2) is discarded
  ///
  /// @param args Coordinates of a bin with respect to the track center
  /// @param cov Track covariance
  ///
  /// @return A value
  template <unsigned int nDim>
  float multivariateGaussian(const Acts::ActsVector<nDim>& args,
                             const Acts::ActsSquareMatrix<nDim>& cov) const;

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
