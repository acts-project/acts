// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <string>
#include <vector>

namespace Acts {

/// @class HoughVertexFinder
///
/// @brief Implements the vertex finder based on the spacepoints using Hough transform
/// For more information, see arXiv:2410.14494
/// 0. Assumes there is only 1 vertex and that it has a high multiplicity
/// 1. Estimates what eta range is really necessary
/// 2. Creates Hough space (z_vtx - cot(theta)) from spacepoints within that eta
/// range
/// 3. Subtracts the coincidentally crossed lines in the Hough space
/// 4. Makes a projection to the Z axis and finds a peak - that is the vertex
/// position
/// 5. Repeats 2-4 if necessary
///
template <typename spacepoint_t>
class HoughVertexFinder {
 public:
  /// Configuration struct
  struct Config {
    /// Ideal amount of spacepoints; |eta| range will be limited to
    /// contain approximately this amount of SPs
    std::uint32_t targetSPs = 10000;

    /// Minimum and maximum ranges in |eta|; the |eta| will not be
    /// set outside these bounds even if targetSPs is not reached
    double minAbsEta = 0.3f;
    double maxAbsEta = 4.0f;

    /// Minimum number of hits in Hough plane to consider
    /// the cell to contain a track
    std::uint32_t minHits = 4;

    /// Number of neighbouring bins in Hough plane
    /// to fill in the cot(theta) direction
    std::uint32_t fillNeighbours = 0;

    /// The algorithm dynamically choose |eta| range necessary to
    /// have a reasonable precision, based on the distribution
    /// of the measurements.
    /// The first |eta| range starts at 0., the others start at
    /// the endpoint of the previous range.
    std::vector<double> absEtaRanges{2., 4.};
    /// The amount of measurements in the |eta| range expressed
    /// as a fraction of all measurements within the whole |eta| range.
    /// Measurements are assumed to be distributed uniformly with
    /// the |eta| range.
    std::vector<double> absEtaFractions{0.4f, 0.6f};

    /// The algorithm starts peak searching considering wide range
    /// in Z and then iteratively constrain itself to narrower ranges
    /// around previously found peak. At the same time it may adapt
    /// other Hough-image parameters in each step.
    std::vector<double> rangeIterZ{200. * UnitConstants::mm,
                                   30. * UnitConstants::mm,
                                   16. * UnitConstants::mm};
    /// Number of bins along z-axis of the Hough image for each iteration
    std::vector<std::uint32_t> nBinsZIterZ{800, 180, 80};
    /// Number of bins along cot(theta)-axis of the Hough image for
    /// each iteration
    std::vector<std::uint32_t> nBinsCotThetaIterZ{8000, 8000, 8000};

    /// If the actual number of measurements is below "targetSPs", then
    /// the number of bins along cot(theta)-axis decreases. Thus, the actual
    /// number of bins can be smaller than stated in "nBinsCotThetaIterZ".
    /// For every magnitude (in natural logarithm) below targetSPs, the number
    /// of bins in cot(theta) will decrease by this factor.
    double binsCotThetaDecrease = 1.35f;

    /// Width of the peak when estimating vertex position
    std::uint32_t peakWidth = 3;

    /// Default position of the vertex in X, Y, and Z coordinates
    Vector3 defVtxPosition{0. * UnitConstants::mm, 0. * UnitConstants::mm,
                           0. * UnitConstants::mm};
  };

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  /// @brief Constructor
  /// @param cfg Configuration object
  /// @param lgr Logging instance
  explicit HoughVertexFinder(
      Config cfg,
      std::unique_ptr<const Logger> lgr = getDefaultLogger("HoughVertexFinder",
                                                           Logging::INFO));

  using HoughCount_t = std::uint16_t;
  using HoughAxis = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  using HoughHist = Grid<HoughCount_t, HoughAxis, HoughAxis>;

  /// @brief Finds the vertex based on the provided spacepoints
  /// @param spacepoints Vector of the input spacepoints; they do not need to be sorted anyhow
  /// @return Position of the vertex
  Acts::Result<Acts::Vector3> find(
      const std::vector<spacepoint_t>& spacepoints) const;

 private:
  /// Configuration instance
  const Config m_cfg;

  /// @brief Returns the positions of the peak along Z axis in the pojection of the Hough plane
  /// @param spacepoints Set of all spacepoints within the event
  /// @param vtxOld Previous position of the vertex
  /// @param rangeZ Range in along Z around vtxOld_z to consider when looking for the new vertex
  /// @param numZBins Number of bins along Z axis
  /// @param minCotTheta Minimum theta to consider for the spacepoint
  /// @param maxCotTheta Maximum theta to consider for the spacepoint
  /// @param numCotThetaBins Number of bins along cot(theta) axis
  /// @return Position of the vertex in (X,Y,Z)
  Acts::Result<Acts::Vector3> findHoughVertex(
      const std::vector<spacepoint_t>& spacepoints, const Acts::Vector3& vtxOld,
      double rangeZ, std::uint32_t numZBins, double minCotTheta,
      double maxCotTheta, std::uint32_t numCotThetaBins) const;

  /// @brief Finds the peak in the Z axis projection of the Hough space
  /// @param houghZProjection Hough space projection after the cleaning procedure
  /// @param vtxZPositions Bins position in the Hough space projection
  /// @return Position of the peak
  Acts::Result<double> findHoughPeak(
      const std::vector<std::uint32_t>& houghZProjection,
      const std::vector<double>& vtxZPositions) const;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "HoughVertexFinder.ipp"
