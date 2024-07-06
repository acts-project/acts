// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Result.hpp"

#include <Eigen/Dense>

namespace Acts {

/// @brief Implements a 1D (no time seeding) grid that is filled with
/// track densities.
///
/// Each track is modelled by a 2D Gaussian distribution in the
/// d0-z0 plane, which is evaluated at d0=0. Therefore, each track
/// effectively lives in 1D.
/// The position of the highest track density (of either a single
/// bin or the sum of a certain region) can be determined.
/// Single tracks can be cached and removed from the overall density.
/// Unlike in the GaussianGridTrackDensity, the overall density vector
/// grows adaptively when tracks densities are added to the grid.
///
/// @tparam trkGridSize The 1-dim grid size of a single track, i.e.
/// a single track is modelled as a (trkGridSize) grid in the z0 axis.
/// Note: trkGridSize has to be an odd value.
template <int trkGridSize = 15>
class AdaptiveGridTrackDensity {
 public:
  static_assert(trkGridSize % 2 == 1, "Assert odd trkGridSize");

  /// z-t position of a maximum and its width
  using ZPositionAndWidth = std::pair<double, double>;

  /// Mapping between bins and single track density
  struct TrackDensityMap {
    using Vector = Eigen::Matrix<double, trkGridSize, 1>;

    /// The central z-bin position of the track density
    std::int32_t centralZBin = 0;
    /// The 1-dim density contribution of the track
    Vector density = Vector::Zero();
  };

  /// Mapping between bins and track densities
  struct MainDensityMap {
    /// The main 1-dim density grid along the z-axis
    std::vector<double> density;
    /// The corresponding z-bin values of the track densities along the z-axis
    std::vector<std::int32_t> zBin;

    std::size_t size() const { return density.size(); }

    bool empty() const { return density.empty(); }

    void clear() {
      density.clear();
      zBin.clear();
    }
  };

  /// The configuration struct
  struct Config {
    Config() = default;

    /// @param binSize_ The binSize in mm
    explicit Config(double binSize_) : binSize(binSize_) {}

    /// Z size of one single bin in grid
    double binSize = 0.1 * UnitConstants::mm;

    /// Do NOT use just the z-bin with the highest
    /// track density, but instead check (up to)
    /// first three density maxima (only those that have
    /// a maximum relative deviation of 'relativeDensityDev'
    /// from the main maximum) and take the z-bin of the
    /// maximum with the highest surrounding density sum
    bool useHighestSumZPosition = false;

    /// The maximum relative density deviation from the main
    /// maximum to consider the second and third maximum for
    /// the highest-sum approach from above
    double maxRelativeDensityDev = 0.01;
  };

  explicit AdaptiveGridTrackDensity(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Returns the z coordinate of maximum (surrounding)
  /// track density
  ///
  /// @param mainDensityMap Map between bins and corresponding density
  /// values
  ///
  /// @return The z coordinate of maximum track density
  Result<double> getMaxZPosition(MainDensityMap& mainDensityMap) const;

  /// @brief Returns the z position of maximum track density
  /// and the estimated width of the maximum
  ///
  /// @param mainDensityMap Map between bins and corresponding density
  /// values
  ///
  /// @return The z position of the maximum track density and
  /// its width
  Result<ZPositionAndWidth> getMaxZPositionAndWidth(
      MainDensityMap& mainDensityMap) const;

  /// @brief Adds a single track to the overall grid density
  ///
  /// @param trk The track to be added
  /// @param mainDensityMap Map between bins and corresponding density
  ///
  /// @return The density map of the track that was added
  TrackDensityMap addTrack(const BoundTrackParameters& trk,
                           MainDensityMap& mainDensityMap) const;

  /// @brief Removes a track from the overall grid density.
  ///
  /// @param trackDensityMap Map between bins and corresponding density
  /// @note The track density comes from a single track
  /// @param mainDensityMap Map between bins and corresponding density
  /// @note The track density comes from an arbitrary number of tracks
  void subtractTrack(const TrackDensityMap& trackDensityMap,
                     MainDensityMap& mainDensityMap) const;

  // TODO this should not be public
  /// @brief Calculates the bin center from the bin number
  /// @param bin Bin number
  /// @return Bin center
  double getBinCenter(std::int32_t bin) const;

 private:
  Config m_cfg;

  /// @brief Calculates the bin number corresponding to a d, z, or time value
  /// @param value d, z, or time value
  /// @return Bin number
  std::int32_t getBin(double value) const;

  /// @brief Function that creates a track density map, i.e., a map from bins
  /// to the corresponding density values for a single track.
  ///
  /// @param centralZBin Central z bin of the track (where its
  /// density is the highest)
  /// @param impactParams vector containing d0, z0 of the track
  /// @param cov 2x2 impact parameter covariance matrix
  ///
  /// @return The track density map
  TrackDensityMap createTrackGrid(std::int32_t centralZBin,
                                  const Vector2& impactParams,
                                  const SquareMatrix2& cov) const;

  /// @brief Function that estimates the seed width in z direction based
  /// on the full width at half maximum (FWHM) of the maximum density peak
  /// @note This only works if the maximum is sufficiently isolated since
  /// overlapping neighboring peaks might lead to an overestimation of the
  /// seed width.
  ///
  /// @param mainDensityMap Map from bins to corresponding track density
  /// @param maxZ z position of the maximum density value
  ///
  /// @return The width
  Result<double> estimateSeedWidth(const MainDensityMap& mainDensityMap,
                                   double maxZ) const;

  /// @brief Helper to retrieve values according to a 2-dim normal distribution
  double normal2D(double d, double z, const SquareMatrix2& cov) const;

  /// @brief Checks (up to) first three density maxima that have a
  /// maximum relative deviation of 'relativeDensityDev' from the
  /// global maximum. Returns the bin of the maximum that has the
  /// highest surrounding density in z direction.
  ///
  /// @param mainDensityMap Map between bins and corresponding density values
  ///
  /// @return The bin index corresponding to the highest surrounding density
  std::size_t highestDensitySumBinIndex(MainDensityMap& mainDensityMap) const;

  /// @brief Calculates the density sum of a bin and its two neighboring bins
  /// in z direction
  ///
  /// @param mainDensityMap Map between bins and corresponding density values
  /// @param zBinIndex Bin index whose neighbors in z we want to sum up
  ///
  /// @return The density sum
  double getDensitySum(const MainDensityMap& mainDensityMap,
                       std::size_t zBinIndex) const;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveGridTrackDensity.ipp"
