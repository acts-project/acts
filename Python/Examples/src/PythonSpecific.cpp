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
#include "ActsExamples/Validation/TrackFinderPerformanceCollector.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <mutex>
#include <stdexcept>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace {

/// A ROOT-free writer that collects track-finder performance histograms and
/// exposes them to Python via histograms() after s.run().
class PythonTrackFinderPerformanceWriter final
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
    /// Plot tool configurations (inlined from TrackFinderPerformanceCollector).
    EffPlotTool::Config effPlotToolConfig;
    FakePlotTool::Config fakePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    TrackQualityPlotTool::Config trackQualityPlotToolConfig;
    /// Optional per-subdetector track summary plots.
    std::map<std::string, std::set<int>> subDetectorTrackSummaryVolumes;
  };

  PythonTrackFinderPerformanceWriter(Config cfg, Acts::Logging::Level lvl)
      : WriterT(cfg.inputTracks, "PythonTrackFinderPerformanceWriter", lvl),
        m_cfg(std::move(cfg)),
        m_collector(
            TrackFinderPerformanceCollector::Config{
                m_cfg.effPlotToolConfig, m_cfg.fakePlotToolConfig,
                m_cfg.duplicationPlotToolConfig,
                m_cfg.trackSummaryPlotToolConfig,
                m_cfg.trackQualityPlotToolConfig,
                m_cfg.subDetectorTrackSummaryVolumes},
            lvl) {
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
    m_collector.logSummary(logger());
    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

  /// Return all filled histograms as a Python dict keyed by histogram name.
  py::dict histograms() const {
    py::dict d;
    const auto& coll = m_collector;

    for (const auto& [name, eff] : coll.effPlotTool().efficiencies1D()) {
      d[py::str(name)] = py::cast(eff, py::return_value_policy::copy);
    }
    for (const auto& [name, eff] : coll.effPlotTool().efficiencies2D()) {
      d[py::str(name)] = py::cast(eff, py::return_value_policy::copy);
    }
    for (const auto& eff : coll.effPlotTool().trackEffVsEtaInPtRanges()) {
      d[py::str(eff.name())] = py::cast(eff, py::return_value_policy::copy);
    }
    for (const auto& eff : coll.effPlotTool().trackEffVsPtInAbsEtaRanges()) {
      d[py::str(eff.name())] = py::cast(eff, py::return_value_policy::copy);
    }

    for (const auto& [name, hist] : coll.fakePlotTool().histograms()) {
      d[py::str(name)] = py::cast(hist, py::return_value_policy::copy);
    }
    for (const auto& [name, eff] : coll.fakePlotTool().efficiencies()) {
      d[py::str(name)] = py::cast(eff, py::return_value_policy::copy);
    }

    for (const auto& [name, prof] : coll.duplicationPlotTool().profiles()) {
      d[py::str(name)] = py::cast(prof, py::return_value_policy::copy);
    }
    for (const auto& [name, eff] : coll.duplicationPlotTool().efficiencies()) {
      d[py::str(name)] = py::cast(eff, py::return_value_policy::copy);
    }

    for (const auto& [name, prof] : coll.trackSummaryPlotTool().profiles()) {
      d[py::str(name)] = py::cast(prof, py::return_value_policy::copy);
    }
    for (const auto& [key, tool] : coll.subDetectorSummaryTools()) {
      for (const auto& [name, prof] : tool.profiles()) {
        d[py::str(name)] = py::cast(prof, py::return_value_policy::copy);
      }
    }

    for (const auto& [name, prof] : coll.trackQualityPlotTool().profiles()) {
      d[py::str(name)] = py::cast(prof, py::return_value_policy::copy);
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
    return m_collector.fill(ctx, tracks, particles, trackParticleMatching,
                            particleTrackMatching, particleMeasurementsMap);
  }

  Config m_cfg;
  std::mutex m_writeMutex;
  TrackFinderPerformanceCollector m_collector;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<ParticleTrackMatching> m_inputParticleTrackMatching{
      this, "InputParticleTrackMatching"};
  ReadDataHandle<InverseMultimap<SimBarcode>> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};
};

}  // namespace

namespace ActsPython {

void addPythonSpecific(py::module_& mex) {
  using Writer = PythonTrackFinderPerformanceWriter;
  using Config = Writer::Config;

  auto w = py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
               mex, "PythonTrackFinderPerformanceWriter")
               .def(py::init<const Config&, Acts::Logging::Level>(),
                    py::arg("config"), py::arg("level"))
               .def_property_readonly("config", &Writer::config)
               .def("histograms", &Writer::histograms);

  auto c = py::class_<Config>(w, "Config").def(py::init<>());
  ACTS_PYTHON_STRUCT(c, inputTracks, inputParticles, inputTrackParticleMatching,
                     inputParticleTrackMatching, inputParticleMeasurementsMap,
                     effPlotToolConfig, fakePlotToolConfig,
                     duplicationPlotToolConfig, trackSummaryPlotToolConfig,
                     trackQualityPlotToolConfig,
                     subDetectorTrackSummaryVolumes);
}

}  // namespace ActsPython
