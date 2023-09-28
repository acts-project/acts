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

#include <boost/functional/hash.hpp>

namespace Acts {

/// @class AdaptiveGridTrackDensity
/// @brief Implements a 1D (no time seeding) / 2D (time seeding)
/// grid that is filled with track densities.
/// Each track is modelled by a 2D / 3D Gaussian distribution in the
/// d0-z0 / d0-z0-t0 plane, which is evaluated at d0=0. Therefore,
/// each track effectively lives in 1D / 2D.
/// The position of the highest track density (of either a single
/// bin or the sum of a certain region) can be determined.
/// Single tracks can be cached and removed from the overall density.
/// Unlike in the GaussianGridTrackDensity, the overall density map
/// grows adaptively when tracks densities are added to the grid.
///
/// @tparam spatialTrkGridSize Number of bins per track in z direction
/// @tparam temporalTrkGridSize Number of bins per track in t direction
/// @note In total, a track is represented by a grid of size
/// spatialTrkGridSize * temporalTrkGridSize
template <int spatialTrkGridSize = 15, int temporalTrkGridSize = 1>
class AdaptiveGridTrackDensity {
  // Assert odd spatial and temporal track grid size
  static_assert(spatialTrkGridSize % 2);
  static_assert(temporalTrkGridSize % 2);

 public:
  // The first (second) integer indicates the bin's z (t) position
  using Bin = std::pair<int, int>;
  // Mapping between bins and track densities
  using DensityMap = std::unordered_map<Bin, float, boost::hash<Bin>>;
  // Coordinates in the z-t plane; the t value will be set to 0 if time
  // vertex seeding is disabled
  using ZTPosition = std::pair<float, float>;
  // z-t position of a maximum and its width
  using ZTPositionAndWidth = std::pair<ZTPosition, float>;

  /// The configuration struct
  struct Config {
    /// @param spatialBinExtent_ The spatial extent of a bin in mm
    Config(float spatialBinExtent_) : spatialBinExtent(spatialBinExtent_) {
      if constexpr (temporalTrkGridSize > 1) {
        throw std::invalid_argument(
            "temporalBinExtent must be provided if temporalTrkGridSize > 1 "
            "(i.e., if time vertex seeding is enabled).");
      }
    }

    /// @param spatialBinExtent_ The spatial extent of a bin in mm
    /// @param temporalBinExtent_ The temporal extent of a bin in mm
    /// @note The speed of light is set to 1, hence the unit.
    Config(float spatialBinExtent_, float temporalBinExtent_)
        : spatialBinExtent(spatialBinExtent_),
          temporalBinExtent(temporalBinExtent_) {
      if constexpr (temporalTrkGridSize == 1) {
        throw std::invalid_argument(
            "temporalBinExtent must not be provided if temporalTrkGridSize == "
            "1 (i.e., if time vertex seeding is disabled).");
      }
    }

    // Spatial extent of a bin in d0 and z0 direction, should always be set to a
    // positive value
    float spatialBinExtent = 0.;  // mm

    // Temporal extent of a bin, should be set to 0 if time vertex seeding is
    // disabled (i.e., if temporalTrkGridSize = 1)
    float temporalBinExtent = 0.;  // mm

    // Do NOT use just the z-bin with the highest
    // track density, but instead check (up to)
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
  static float getBinCenter(int bin, float binExtent);

  /// @brief Calculates the bin number corresponding to a d, z, or time value
  /// @param value d, z, or time value
  /// @param binExtent Bin extent
  /// @return Bin number
  static int getBin(float value, float binExtent);

  /// @brief Finds the maximum density of a DensityMap
  /// @param densityMap Map between bins and corresponding density
  /// values
  /// @return Iterator of the map entry with the highest density
  DensityMap::const_iterator highestDensityEntry(
      const DensityMap& densityMap) const;

  /// @brief Returns the z and t coordinate of maximum (surrounding)
  /// track density
  /// @note if time vertex seeding is not enabled, the t coordinate
  /// will be set to 0.
  ///
  /// @param densityMap Map between bins and corresponding density
  /// values
  ///
  /// @return The z and t coordinates of maximum track density
  Result<ZTPosition> getMaxZTPosition(DensityMap& densityMap) const;

  /// @brief Returns the z-t position of maximum track density
  /// and the estimated z-width of the maximum
  ///
  /// @param densityMap Map between bins and corresponding density
  /// values
  ///
  /// @return The z-t position of the maximum track density and
  /// its width
  Result<ZTPositionAndWidth> getMaxZTPositionAndWidth(
      DensityMap& densityMap) const;

  /// @brief Adds a single track to the overall grid density
  ///
  /// @param trk The track to be added
  /// @param mainDensityMap Map between bins and corresponding density
  ///
  /// @return The density map of the track that was added
  DensityMap addTrack(const BoundTrackParameters& trk,
                      DensityMap& mainDensityMap) const;

  /// @brief Removes a track from the overall grid density.
  ///
  /// @param trackDensityMap Map between bins and corresponding density
  /// @note The track density comes from a single track
  /// @param mainDensityMap Map between bins and corresponding density
  /// @note The track density comes from an arbitrary number of tracks
  void subtractTrack(const DensityMap& trackDensityMap,
                     DensityMap& mainDensityMap) const;

 private:
  /// @brief Function that creates a track density map, i.e., a map from bins
  /// to the corresponding density values for a single track.
  ///
  /// @param impactParams vector containing d0, z0, and t0 of the track
  /// @param centralBin Central z and t bin of the track (where its
  /// density is the highest)
  /// @param cov 3x3 impact parameter covariance matrix
  DensityMap createTrackGrid(const Acts::Vector3& impactParams,
                             const Bin& centralBin,
                             const Acts::SquareMatrix3& cov) const;

  /// @brief Function that estimates the seed width in z direction based
  /// on the full width at half maximum (FWHM) of the maximum density peak
  /// @note This only works if the maximum is sufficiently isolated since
  /// overlapping neighboring peaks might lead to an overestimation of the
  /// seed width.
  ///
  /// @param densityMap Map from bins to corresponding track density
  /// @param maxZT z-t position of the maximum density value
  ///
  /// @return The width
  Result<float> estimateSeedWidth(const DensityMap& densityMap,
                                  const ZTPosition& maxZT) const;

  /// @brief Helper to retrieve values of an nDim-dimensional normal
  /// distribution
  /// @note The constant prefactor (2 * pi)^(- nDim / 2) is discarded
  ///
  /// @param args Coordinates where the Gaussian should be evaluated
  /// @note args must be in a coordinate system with origin at the mean
  /// values of the Gaussian
  /// @param cov Covariance matrix
  ///
  /// @return Multivariate Gaussian evaluated at args
  template <unsigned int nDim>
  static float multivariateGaussian(const Acts::ActsVector<nDim>& args,
                                    const Acts::ActsSquareMatrix<nDim>& cov);

  /// @brief Checks (up to) first three density maxima that have a
  /// maximum relative deviation of 'relativeDensityDev' from the
  /// global maximum. Returns the bin of the maximum that has the
  /// highest surrounding density in z direction.
  ///
  /// @param densityMap Map between bins and corresponding density values
  ///
  /// @return The bin corresponding to the highest surrounding density
  Bin highestDensitySumBin(DensityMap& densityMap) const;

  /// @brief Calculates the density sum of a bin and its two neighboring bins
  /// in z direction
  ///
  /// @param densityMap Map between bins and corresponding density values
  /// @param bin Bin whose neighbors in z we want to sum up
  ///
  /// @return The density sum
  float getDensitySum(const DensityMap& densityMap, const Bin& bin) const;

  Config m_cfg;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveGridTrackDensity.ipp"
