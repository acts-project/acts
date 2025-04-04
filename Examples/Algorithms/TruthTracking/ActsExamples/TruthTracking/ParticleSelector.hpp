// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <limits>
#include <string>

namespace ActsExamples {
struct AlgorithmContext;

/// Select particles by applying some selection cuts.
class ParticleSelector final : public IAlgorithm {
 public:
  struct MeasurementCounter {
    // Combination of a geometry hierarchy map and a minimum hit count
    using CounterElement =
        std::pair<Acts::GeometryHierarchyMap<unsigned int>, unsigned int>;

    boost::container::small_vector<CounterElement, 4> counters;

    bool isValidParticle(
        const SimParticle& particle,
        const InverseMultimap<SimBarcode>& particleMeasurementsMap,
        const MeasurementContainer& measurements) const;

    void addCounter(const std::vector<Acts::GeometryIdentifier>& identifiers,
                    unsigned int threshold) {
      std::vector<Acts::GeometryHierarchyMap<unsigned int>::InputElement>
          elements;
      for (const auto& id : identifiers) {
        elements.emplace_back(id, 0);
      }
      counters.emplace_back(std::move(elements), threshold);
    }
  };

  struct Config {
    /// The input particles collection.
    std::string inputParticles;
    /// (Optionally) The input particle measurements map. Only required for
    /// measurement-based cuts.
    std::string inputParticleMeasurementsMap;
    /// (Optionally) The input measurements collection. Only required for
    /// measurement-based cuts.
    std::string inputMeasurements;
    /// The output particles collection.
    std::string outputParticles;

    // Minimum/maximum distance from the origin in the transverse plane.
    double rhoMin = 0;
    double rhoMax = std::numeric_limits<double>::infinity();
    // Minimum/maximum absolute distance from the origin along z.
    double absZMin = 0;
    double absZMax = std::numeric_limits<double>::infinity();
    // Minimum/maximum particle time.
    double timeMin = -std::numeric_limits<double>::infinity();
    double timeMax = std::numeric_limits<double>::infinity();
    // Direction cuts.
    double phiMin = -std::numeric_limits<double>::infinity();
    double phiMax = std::numeric_limits<double>::infinity();
    double etaMin = -std::numeric_limits<double>::infinity();
    double etaMax = std::numeric_limits<double>::infinity();
    double absEtaMin = 0;
    double absEtaMax = std::numeric_limits<double>::infinity();
    // Momentum cuts.
    double ptMin = 0;
    double ptMax = std::numeric_limits<double>::infinity();
    // Rest mass cuts
    double mMin = 0;
    double mMax = std::numeric_limits<double>::infinity();
    // Hit count cuts
    std::size_t hitsMin = 0;
    std::size_t hitsMax = std::numeric_limits<std::size_t>::max();
    // Measurement number cuts
    std::size_t measurementsMin = 0;
    std::size_t measurementsMax = std::numeric_limits<std::size_t>::max();
    /// Remove charged particles.
    bool removeCharged = false;
    /// Remove neutral particles.
    bool removeNeutral = false;
    /// Remove secondaries.
    bool removeSecondaries = false;
    /// Exclude particles depending on absolute pdg value
    std::vector<int> excludeAbsPdgs;

    /// Min primary vertex ID cut
    std::uint64_t minPrimaryVertexId = 0;
    /// Max primary vertex ID cut
    std::uint64_t maxPrimaryVertexId =
        std::numeric_limits<std::uint64_t>::max();

    /// The measurement counter to be used for the measurement cuts.
    MeasurementCounter measurementCounter;
  };

  ParticleSelector(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<InverseMultimap<SimBarcode>> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
};

}  // namespace ActsExamples
