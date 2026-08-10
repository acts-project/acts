// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/GaussianHistogramFit.hpp"
#include "ActsExamples/Validation/HistogramFit.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace ActsExamples {

/// Collects track-fitter performance histograms without any file I/O.
///
/// Collects residual/pull histograms, efficiency plots, and track summary
/// information for track fitting performance evaluation. The Gaussian fit
/// backend is supplied by the caller via @c Config::fitFunction, so this
/// collector is agnostic to whether it runs Core's own fit, ROOT's `TH1::Fit`,
/// or a Python callable.
///
/// @note The caller must ensure exclusive access (e.g. hold a mutex) when
///       calling fill(). This class applies no locking of its own.
class TrackFitterPerformanceCollector {
 public:
  struct Config {
    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// The Gaussian fit backend used by @c fitProfiles. Defaults to Core's
    /// own ROOT-free @c gaussianHistogramFit; pass e.g.
    /// @c ActsPlugins::RootHistogramFit or a Python callable instead to use a
    /// different backend. If explicitly cleared, @c fitProfiles logs a
    /// warning and returns no profiles instead of fitting.
    HistogramFitFunction fitFunction = &gaussianHistogramFit;

    /// Minimum number of entries in a bin for it to be included in the
    /// mean/width fit.
    int fitMinEntries = 10;
    /// The range in sigma for the iterative Gaussian fit
    double fitSigmaRange = 3.0;
    /// The maximum number of iterations for the iterative Gaussian fit
    int fitIterations = 3;
    /// Threshold for warning about fit failure fraction in profile
    /// extraction.
    double warningThresholdFitFailureFraction = 0.55;
  };

  /// @param cfg The configuration
  /// @param logger Logger, also used for diagnostics from @c fitProfiles
  TrackFitterPerformanceCollector(Config cfg,
                                  std::unique_ptr<const Acts::Logger> logger);

  /// Fill histograms for one event.
  ///
  /// @note The caller must ensure exclusive access (e.g. hold a mutex).
  void fill(const Acts::GeometryContext& geoContext,
            const ConstTrackContainer& tracks,
            const SimParticleContainer& particles,
            const TrackParticleMatching& trackParticleMatching);

  /// Summary count statistics accumulated across all filled events.
  struct Stats {
    std::size_t nTotalTracks = 0;
    std::size_t nTotalMatchedTracks = 0;
    std::size_t nTotalFakeTracks = 0;
    std::size_t nTotalParticles = 0;
    std::size_t nTotalMatchedParticles = 0;
  };

  /// Return accumulated event counts.
  const Stats& stats() const { return m_stats; }

  /// Emit summary statistics via the internal logger.
  void logSummary() const;

  /// @name Accessors for the underlying plot tools
  /// @{
  const ResPlotTool& resPlotTool() const { return m_resPlotTool; }
  const EffPlotTool& effPlotTool() const { return m_effPlotTool; }
  const TrackSummaryPlotTool& trackSummaryPlotTool() const {
    return m_trackSummaryPlotTool;
  }
  /// @}

  /// Mean/width profiles fitted from every residual and pull histogram.
  ///
  /// @c profiles1 holds the outputs of the 2D (vs. eta, vs. pT) inputs;
  /// @c profiles2 the 3D (vs. eta-phi, vs. eta-pT) ones. Each profile carries
  /// its own name, e.g. `"resmean_d0_vs_eta"` / `"reswidth_d0_vs_eta"`.
  struct FittedProfiles {
    std::vector<Acts::Experimental::ValueHistogram1> profiles1;
    std::vector<Acts::Experimental::ValueHistogram2> profiles2;
  };

  /// Fit every residual/pull profile histogram with @c Config::fitFunction.
  ///
  /// Emits a warning via the internal logger for any input histogram whose
  /// fit failure fraction reaches @c Config::warningThresholdFitFailureFraction.
  /// If @c Config::fitFunction is unset, logs a warning and returns no
  /// profiles instead.
  FittedProfiles fitProfiles() const;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  /// Fit every histogram in @p histMap and append the resulting mean/width
  /// profiles to @p out, warning on excessive fit failures.
  template <std::size_t Dim>
  void addFittedProfiles(
      const std::map<std::string, Acts::Experimental::Histogram<Dim>>& histMap,
      const std::string& meanPrefix, const std::string& widthPrefix,
      std::vector<Acts::Experimental::ValueHistogram<Dim - 1>>& out) const;

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ResPlotTool m_resPlotTool;
  EffPlotTool m_effPlotTool;
  TrackSummaryPlotTool m_trackSummaryPlotTool;

  Stats m_stats;
};

}  // namespace ActsExamples
