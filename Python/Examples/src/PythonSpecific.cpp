// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/PatternRecognitionPerformanceCollector.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackParameterPerformanceCollector.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <mutex>
#include <stdexcept>
#include <string>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace {

/// Insert @p obj into @p d under @p name, throwing if the name is already
/// present. Histogram names are expected to be unique; a collision indicates
/// a bug in the plot-tool configuration rather than an intentional overwrite.
template <typename T>
void insertUniqueHistogram(py::dict& d, const std::string& name, const T& obj) {
  py::str key(name);
  if (d.contains(key)) {
    throw std::runtime_error("Duplicate histogram name: " + name);
  }
  d[key] = py::cast(obj, py::return_value_policy::copy);
}

/// A ROOT-free writer that collects pattern-recognition performance histograms
/// and exposes them to Python via histograms() after s.run().
class PythonPatternRecognitionPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (found) tracks collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input track-particle matching.
    std::string inputParticleTrackMatching;
    /// Input particle measurements map.
    std::string inputParticleMeasurementsMap;
    /// Label for histogram titles and names.
    std::string label = "track";
    /// Plot tool configurations (inlined from
    /// PatternRecognitionPerformanceCollector).
    EffPlotTool::Config effPlotToolConfig;
    FakePlotTool::Config fakePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    TrackQualityPlotTool::Config trackQualityPlotToolConfig;
    /// Optional per-subdetector track summary plots.
    std::map<std::string, std::set<int>> subDetectorTrackSummaryVolumes;
  };

  PythonPatternRecognitionPerformanceWriter(Config cfg,
                                            Acts::Logging::Level lvl)
      : WriterT(cfg.inputTracks, "PythonPatternRecognitionPerformanceWriter",
                lvl),
        m_cfg(std::move(cfg)),
        m_collector(
            PatternRecognitionPerformanceCollector::Config{
                m_cfg.label, m_cfg.effPlotToolConfig, m_cfg.fakePlotToolConfig,
                m_cfg.duplicationPlotToolConfig,
                m_cfg.trackSummaryPlotToolConfig,
                m_cfg.trackQualityPlotToolConfig,
                m_cfg.subDetectorTrackSummaryVolumes},
            logger().clone()) {
    if (m_cfg.inputParticles.empty()) {
      throw std::invalid_argument("Missing particles input collection");
    }
    if (m_cfg.inputTrackParticleMatching.empty()) {
      throw std::invalid_argument("Missing input track particles matching");
    }
    if (m_cfg.inputParticleTrackMatching.empty()) {
      throw std::invalid_argument("Missing input particle track matching");
    }
    if (m_cfg.inputParticleMeasurementsMap.empty()) {
      throw std::invalid_argument("Missing input measurement particles map");
    }

    m_inputParticles.initialize(m_cfg.inputParticles);
    m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
    m_inputParticleTrackMatching.initialize(m_cfg.inputParticleTrackMatching);
    m_inputParticleMeasurementsMap.initialize(
        m_cfg.inputParticleMeasurementsMap);
  }

  ProcessCode finalize() override {
    m_collector.logSummary();
    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

  /// Return all filled histograms as a Python dict keyed by histogram name.
  py::dict histograms() const {
    py::dict d;
    const auto& coll = m_collector;

    for (const auto& [name, eff] : coll.effPlotTool().efficiencies1D()) {
      insertUniqueHistogram(d, name, eff);
    }
    for (const auto& [name, eff] : coll.effPlotTool().efficiencies2D()) {
      insertUniqueHistogram(d, name, eff);
    }
    for (const auto& eff : coll.effPlotTool().trackEffVsEtaInPtRanges()) {
      insertUniqueHistogram(d, eff.name(), eff);
    }
    for (const auto& eff : coll.effPlotTool().trackEffVsPtInAbsEtaRanges()) {
      insertUniqueHistogram(d, eff.name(), eff);
    }

    for (const auto& [name, hist] : coll.fakePlotTool().histograms()) {
      insertUniqueHistogram(d, name, hist);
    }
    for (const auto& [name, eff] : coll.fakePlotTool().efficiencies()) {
      insertUniqueHistogram(d, name, eff);
    }

    for (const auto& [name, prof] : coll.duplicationPlotTool().profiles()) {
      insertUniqueHistogram(d, name, prof);
    }
    for (const auto& [name, eff] : coll.duplicationPlotTool().efficiencies()) {
      insertUniqueHistogram(d, name, eff);
    }

    for (const auto& [name, prof] : coll.trackSummaryPlotTool().profiles()) {
      insertUniqueHistogram(d, name, prof);
    }
    for (const auto& [key, tool] : coll.subDetectorSummaryTools()) {
      for (const auto& [name, prof] : tool.profiles()) {
        insertUniqueHistogram(d, name, prof);
      }
    }

    for (const auto& [name, prof] : coll.trackQualityPlotTool().profiles()) {
      insertUniqueHistogram(d, name, prof);
    }

    return d;
  }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override {
    const auto& particles = m_inputParticles(ctx);
    const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);
    const auto& particleTrackMatching = m_inputParticleTrackMatching(ctx);
    const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);

    std::lock_guard<std::mutex> lock(m_writeMutex);
    m_collector.fill(ctx.recoGeoContext, tracks, particles,
                     trackParticleMatching, particleTrackMatching,
                     particleMeasurementsMap);
    return ProcessCode::SUCCESS;
  }

  Config m_cfg;
  std::mutex m_writeMutex;
  PatternRecognitionPerformanceCollector m_collector;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<ParticleTrackMatching> m_inputParticleTrackMatching{
      this, "InputParticleTrackMatching"};
  ReadDataHandle<InverseMultimap<SimBarcode>> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};
};

}  // namespace

/// A ROOT-free writer that collects track-parameter performance histograms and
/// exposes them to Python via histograms() after s.run().
class PythonTrackParameterPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input tracks collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Output filename (optional).
    std::string filePath = "performance_track_parameters.root";
    /// Plot tool configurations.
    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    /// The Gaussian fit backend, a Python callable with the matching
    /// signature. Required: there is no sensible default.
    HistogramFitFunction fitFunction;
    /// Fit parameters.
    int fitMinEntries = 10;
    double fitSigmaRange = 3.0;
    int fitIterations = 3;
    double warningThresholdFitFailureFraction = 0.55;
  };

  /// Translate the writer configuration into the collector configuration.
  static TrackParameterPerformanceCollector::Config collectorConfig(
      const Config& cfg) {
    TrackParameterPerformanceCollector::Config collectorCfg;
    collectorCfg.resPlotToolConfig = cfg.resPlotToolConfig;
    collectorCfg.effPlotToolConfig = cfg.effPlotToolConfig;
    collectorCfg.trackSummaryPlotToolConfig = cfg.trackSummaryPlotToolConfig;
    collectorCfg.fitFunction = cfg.fitFunction;
    collectorCfg.fitMinEntries = cfg.fitMinEntries;
    collectorCfg.fitSigmaRange = cfg.fitSigmaRange;
    collectorCfg.fitIterations = cfg.fitIterations;
    collectorCfg.warningThresholdFitFailureFraction =
        cfg.warningThresholdFitFailureFraction;
    return collectorCfg;
  }

  PythonTrackParameterPerformanceWriter(Config cfg, Acts::Logging::Level lvl)
      : WriterT(cfg.inputTracks, "PythonTrackParameterPerformanceWriter", lvl),
        m_cfg(std::move(cfg)),
        m_collector(collectorConfig(m_cfg), logger().clone()) {
    if (m_cfg.inputParticles.empty()) {
      throw std::invalid_argument("Missing particles input collection");
    }
    if (m_cfg.inputTrackParticleMatching.empty()) {
      throw std::invalid_argument("Missing input track particles matching");
    }

    m_inputParticles.initialize(m_cfg.inputParticles);
    m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  }

  ProcessCode finalize() override {
    m_collector.logSummary();
    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

  /// Return all filled histograms as a Python dict keyed by histogram name.
  py::dict histograms() const {
    py::dict d;
    const auto& coll = m_collector;

    // Residual histograms
    for (const auto& [name, hist] : coll.resPlotTool().res()) {
      insertUniqueHistogram(d, "res_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().resVsEta()) {
      insertUniqueHistogram(d, "resVsEta_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().resVsPt()) {
      insertUniqueHistogram(d, "resVsPt_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().resVsEtaPhi()) {
      insertUniqueHistogram(d, "resVsEtaPhi_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().resVsEtaPt()) {
      insertUniqueHistogram(d, "resVsEtaPt_" + name, hist);
    }

    // Pull histograms
    for (const auto& [name, hist] : coll.resPlotTool().pull()) {
      insertUniqueHistogram(d, "pull_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().pullVsEta()) {
      insertUniqueHistogram(d, "pullVsEta_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().pullVsPt()) {
      insertUniqueHistogram(d, "pullVsPt_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().pullVsEtaPhi()) {
      insertUniqueHistogram(d, "pullVsEtaPhi_" + name, hist);
    }
    for (const auto& [name, hist] : coll.resPlotTool().pullVsEtaPt()) {
      insertUniqueHistogram(d, "pullVsEtaPt_" + name, hist);
    }

    // Efficiency histograms
    for (const auto& [name, eff] : coll.effPlotTool().efficiencies1D()) {
      insertUniqueHistogram(d, name, eff);
    }
    for (const auto& [name, eff] : coll.effPlotTool().efficiencies2D()) {
      insertUniqueHistogram(d, name, eff);
    }
    for (const auto& eff : coll.effPlotTool().trackEffVsEtaInPtRanges()) {
      insertUniqueHistogram(d, eff.name(), eff);
    }
    for (const auto& eff : coll.effPlotTool().trackEffVsPtInAbsEtaRanges()) {
      insertUniqueHistogram(d, eff.name(), eff);
    }

    // Track summary histograms
    for (const auto& [name, prof] : coll.trackSummaryPlotTool().profiles()) {
      insertUniqueHistogram(d, name, prof);
    }

    // Fitted mean/width profiles
    const auto fittedProfiles = coll.fitProfiles();
    for (const auto& profile : fittedProfiles.profiles1) {
      insertUniqueHistogram(d, profile.name(), profile);
    }
    for (const auto& profile : fittedProfiles.profiles2) {
      insertUniqueHistogram(d, profile.name(), profile);
    }

    return d;
  }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override {
    const auto& particles = m_inputParticles(ctx);
    const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);

    std::lock_guard<std::mutex> lock(m_writeMutex);
    m_collector.fill(ctx.recoGeoContext, tracks, particles,
                     trackParticleMatching);
    return ProcessCode::SUCCESS;
  }

  Config m_cfg;
  std::mutex m_writeMutex;
  TrackParameterPerformanceCollector m_collector;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
};

namespace ActsPython {

void addPythonSpecific(py::module_& mex) {
  {
    using Writer = PythonPatternRecognitionPerformanceWriter;
    using Config = Writer::Config;

    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "PythonPatternRecognitionPerformanceWriter")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly("config", &Writer::config)
                 .def("histograms", &Writer::histograms);

    auto c = py::class_<Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, inputTracks, inputParticles,
                       inputTrackParticleMatching, inputParticleTrackMatching,
                       inputParticleMeasurementsMap, label, effPlotToolConfig,
                       fakePlotToolConfig, duplicationPlotToolConfig,
                       trackSummaryPlotToolConfig, trackQualityPlotToolConfig,
                       subDetectorTrackSummaryVolumes);
  }

  {
    using Writer = PythonTrackParameterPerformanceWriter;
    using Config = Writer::Config;

    auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
                 mex, "PythonTrackParameterPerformanceWriter")
                 .def(py::init<const Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def_property_readonly("config", &Writer::config)
                 .def("histograms", &Writer::histograms);

    auto c = py::class_<Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, inputTracks, inputParticles,
                       inputTrackParticleMatching, filePath, resPlotToolConfig,
                       effPlotToolConfig, trackSummaryPlotToolConfig,
                       fitFunction, fitMinEntries, fitSigmaRange, fitIterations,
                       warningThresholdFitFailureFraction);
  }
}

}  // namespace ActsPython
