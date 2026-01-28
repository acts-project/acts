// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

/// Tools to make duplication ratio (formerly called duplication rate) and
/// duplication number plots to show tracking duplication.
///
/// The duplication ratio is evaluated on truth-matched reco tracks. If there
/// is more than one track matched to the same truth particle, the reco track
/// with the highest matching probability is tagged as 'real' and the others are
/// 'duplicated'.
class DuplicationPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
        {"Eta", PlotHelpers::Binning::Uniform("#eta", 40, -4, 4)},
        {"Phi", PlotHelpers::Binning::Uniform("#phi", 100, -3.15, 3.15)},
        {"Pt", PlotHelpers::Binning::Uniform("pT [GeV/c]", 40, 0, 100)},
        {"Num", PlotHelpers::Binning::Uniform("N", 30, -0.5, 29.5)}};
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  DuplicationPlotTool(const Config& cfg, Acts::Logging::Level lvl);
  DuplicationPlotTool(const DuplicationPlotTool&) = delete;
  DuplicationPlotTool& operator=(const DuplicationPlotTool&) = delete;
  DuplicationPlotTool(DuplicationPlotTool&&) noexcept = default;
  DuplicationPlotTool& operator=(DuplicationPlotTool&&) noexcept = default;
  ~DuplicationPlotTool();

  /// @brief book the duplication plots
  ///
  void book();

  /// @brief fill duplication ratio w.r.t. fitted track parameters
  ///
  /// @param fittedParameters fitted track parameters of this track
  /// @param status the (truth-matched) reconstructed track is duplicated or not
  void fill(const Acts::BoundTrackParameters& fittedParameters, bool status);

  /// @brief fill number of duplicated tracks for a truth particle seed
  ///
  /// @param truthParticle the truth Particle
  /// @param nMatchedTracks the number of matched tracks
  void fill(const SimParticleState& truthParticle, std::size_t nMatchedTracks);

  /// @brief write the duplication plots to file
  ///
  void write();

 private:
  struct Impl;

  /// The Config class
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Impl> m_impl;

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
