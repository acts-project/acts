// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MuonHoughMaximum.hpp"
#include "ActsExamples/EventData/MuonSegment.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/// @brief Example implementation of a muon hough transform seeder
/// Uses the hough tools from the ACTS Core repo
/// Reads CSV files with muon sim hits (= true trajectories)
/// and drift circles (= measurements), performs
/// a hough transform to the drift circles in each station,
/// and compares to the true parameters of the sim hit in the
/// given station.
class MuonHoughSeeder final : public IAlgorithm {
 public:
  /// @brief Abbrivation of the HoughPlane_t
  using HoughPlane_t =
      Acts::HoughTransformUtils::HoughPlane<const MuonSpacePoint*>;
  /// @brief Abbrivation of the PeakFinder
  using PeakFinder_t = Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
      const MuonSpacePoint*>;
  /// @brief Abbrivation of the PeakFinder configuration object
  using PeakFinderCfg_t =
      Acts::HoughTransformUtils::PeakFinders::IslandsAroundMaxConfig;
  /// @brief Abbrivation of the HoughMaximum type returned by the PeakFinder
  using Maximum_t = PeakFinder_t::Maximum;
  /// @brief Abbrivation of the HoughMaximum vector
  using MaximumVec_t = std::vector<Maximum_t>;
  /// @brief Abbrivation of the HoughTransform axis utils
  using AxisRange_t = Acts::HoughTransformUtils::HoughAxisRanges;
  /// @brief Abbrivation of the space point id
  using MuonId = MuonSpacePoint::MuonId;

  /// @brief Configuration object of the Hough seeder
  struct Config {
    /// @brief Container name of the truth segments (used for validation)
    std::string inTruthSegments{};
    /// @brief Container name of the space point collection
    std::string inSpacePoints{};
    /// @brief Container name of the output hough seed collection
    std::string outHoughMax{};
    /// @brief Extra margin added to both y-sides of the eta-hough accumulator plane
    double etaPlaneMarginIcept{10. * Acts::UnitConstants::cm};
    /// @brief Extra margin added to both y-sides of the phi-hough accumulator plane
    double phiPlaneMarginIcept{10. * Acts::UnitConstants::cm};
    /// @brief Number of bins to scan tan (theta)
    unsigned nBinsTanTheta = 25;
    /// @brief Number of bins in y0 space (complementary to tan (theta))
    unsigned nBinsY0 = 25;
    /// @brief Number of bins to scan tan (phi)
    unsigned nBinsTanPhi = 10;
    /// @brief Number of bins in x0 space (omplementary to tan (phi))
    unsigned nBinsX0 = 10;
    /// @brief Visualize the Hough plane maxima
    bool dumpVisualization{false};
    /// @brief Visualization function (optional, e.g., for ROOT-based visualization)
    /// Takes: outputPath, bucketId, maxima, plane, axis, truthSegments, logger
    std::function<void(const std::string&, const MuonId&, const MaximumVec_t&,
                       const HoughPlane_t&, const AxisRange_t&,
                       const MuonSegmentContainer&, const Acts::Logger&)>
        visualizationFunction{};
  };

  MuonHoughSeeder(const Config& cfg, Acts::Logging::Level lvl);
  ~MuonHoughSeeder() override;

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief Find eta maxima from the space point bucket and fills them into a new
  ///        maximum container
  /// @param ctx: Algorithm context needed for the display of the truth-parameters
  /// @param bucket: Spacepoint bucket of interest
  /// @param plane: Allocated hough plane to be recycled for all hough searches in the event
  MuonHoughMaxContainer constructEtaMaxima(const AlgorithmContext& ctx,
                                           const MuonSpacePointBucket& bucket,
                                           HoughPlane_t& plane) const;
  /// @brief Extends the obtained eta maxima and tries to attach straight line parameters
  ///        in the non-precision plane (phi).
  MuonHoughMaxContainer extendMaximaWithPhi(const AlgorithmContext& ctx,
                                            MuonHoughMaxContainer&& etaMaxima,
                                            HoughPlane_t& plane) const;

  /// @brief Dispatches visualization of the found maxima via the configured
  ///        visualization function (no-op if not configured)
  /// @param ctx: Algorithm context to fetch the truth segment parameters
  /// @param bucketId: identifier of the bucket and whether it's an eta / phi maximum
  /// @param maxima: List of maxima from the PeakFinder
  /// @param plane: Filled hough plane
  /// @param axis: Axis range needed to interpret the hough binning
  void displayMaxima(const AlgorithmContext& ctx, const MuonId& bucketId,
                     const MaximumVec_t& maxima, const HoughPlane_t& plane,
                     const AxisRange_t& axis) const;

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }

  ReadDataHandle<MuonSegmentContainer> m_inputTruthSegs{this,
                                                        "InputTruthSegments"};
  ReadDataHandle<MuonSpacePointContainer> m_inputSpacePoints{
      this, "InputSpacePoints"};
  WriteDataHandle<MuonHoughMaxContainer> m_outputMaxima{this, "OutputHoughMax"};
};

}  // namespace ActsExamples
