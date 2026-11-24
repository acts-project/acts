// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Manager.hpp"
#include "ActsExamples/Geant4/RunManager.hpp"

#include <cstdint>
#include <functional>
#include <memory>

class G4VPhysicalVolume;

namespace ActsExamples::Geant4 {

/// A single, serial transport session sharing the application's Geant4 handle.
class SimulationSession {
 public:
  SimulationSession(std::shared_ptr<Geant4Handle> handle,
                    const std::string& physicsList, double largestStep);

  SimulationSession(const SimulationSession&) = delete;
  SimulationSession& operator=(const SimulationSession&) = delete;

  void install();
  void setMagneticField(
      G4VPhysicalVolume& world,
      std::shared_ptr<const Acts::MagneticFieldProvider> field);
  void beginRun();
  EventResult simulate(const SimParticleContainer& particles,
                       const Acts::GeometryContext& geometry,
                       std::uint64_t seed);
  void endRun();

  const std::shared_ptr<Geant4Handle> handle;
  const std::shared_ptr<EventContext> eventContext =
      std::make_shared<EventContext>();
  /// Install optional actions once before the first run. A failed callback
  /// keeps the session available for another initialization attempt.
  std::function<void()> configureActions;

 private:
  RunManager& m_manager;
  double m_largestStep;
  bool m_runStarted = false;
};

}  // namespace ActsExamples::Geant4
