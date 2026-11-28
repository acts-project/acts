// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticleDecay.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace ActsExamples::Geant4 {

/// A persistent Geant4 master and native worker pool for one ACTS sequence.
class NativeSimulationSession final {
 public:
  NativeSimulationSession(const Geant4Simulation::Config& cfg,
                          const Acts::Logger& logger);
  ~NativeSimulationSession();

  void beginRun();
  void endRun();
  EventResult simulate(const SimParticleContainer& particles,
                       const Acts::GeometryContext& geometry,
                       std::uint64_t seed, std::size_t eventNumber,
                       const SimParticleDecayGraph* decayGraph);

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples::Geant4
