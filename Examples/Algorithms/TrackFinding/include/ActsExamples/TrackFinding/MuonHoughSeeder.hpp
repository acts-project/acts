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
#include <memory>
#include <string>
#include <vector>

class TCanvas;

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
  /// config
  struct Config {
    std::string inTruthSegments{};
    std::string inSpacePoints{};
    std::string outHoughMax{};

    /** @brief Extra margin added to both y-sides of the eta-hough accumulator plane */
    double etaPlaneMarginIcept{10. * Acts::UnitConstants::cm};
    /** @brief Extra margin added to both y-sides of the phi-hough accumulator plane */
    double phiPlaneMarginIcept{10. * Acts::UnitConstants::cm};
  };

  MuonHoughSeeder(Config cfg, Acts::Logging::Level lvl);
  ~MuonHoughSeeder() override;

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;
  ProcessCode initialize() final;
  ProcessCode finalize() final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /** @brief Abbrivation of the HoughPlane_t */
  using HoughPlane_t =
      Acts::HoughTransformUtils::HoughPlane<const MuonSpacePoint*>;
  /** @brief Abbrivation of the PeakFinder */
  using PeakFinder_t = Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
      const MuonSpacePoint*>;
  /** @brief Abbrivation of the PeakFinder configuration object  */
  using PeakFinderCfg_t =
      Acts::HoughTransformUtils::PeakFinders::IslandsAroundMaxConfig;
  /** @brief Abbrivation of the HoughMaximum type returned by the PeakFinder */
  using Maximum_t = PeakFinder_t::Maximum;
  /** @brief Abbrivation of the HoughMaximum vector */
  using MaximumVec_t = std::vector<Maximum_t>;
  /** @brief Abbrivation of the HoughTransform axis utils */
  using AxisRange_t = Acts::HoughTransformUtils::HoughAxisRanges;
  /** @brief Abbrivation of the space point id */
  using MuonId = MuonSpacePoint::MuonId;

  /** @brief Find eta maxima from the space point bucket and fills them into a new
   *         maximum container
   *  @param ctx: Algorithm context needed for the display of the truth-parameters
   *  @param bucket: Spacepoint bucket of interest
   *  @param plane: Allocated hough plane to be recycled for all hough searches in the event*/
  MuonHoughMaxContainer constructEtaMaxima(const AlgorithmContext& ctx,
                                           const MuonSpacePointBucket& bucket,
                                           HoughPlane_t& plane) const;
  /** @brief Extends the obtained eta maxima and tries to attach straight line parameters
   *         in the non-precision plane (phi).
   */
  MuonHoughMaxContainer extendMaximaWithPhi(const AlgorithmContext& ctx,
                                            MuonHoughMaxContainer&& etaMaxima,
                                            HoughPlane_t& plane) const;

  /** @brief Displays the found maxima onto a TCanvas
   *  @param ctx: Algorithm context to fetch the truth segment parameters
   *  @param bucketId: identifier of the bucket to display on the Canvas and also
   *                   to determine whether it's an eta / phi maximum
   *  @param maxima: List of maxima from the PeakFinder
   *  @param plane: Filled hough plane
   *  @param axis: Axis range needed to interpret the hough binning */
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
  /// use ROOT for visualisation
  std::unique_ptr<TCanvas> m_outCanvas;
};

}  // namespace ActsExamples
