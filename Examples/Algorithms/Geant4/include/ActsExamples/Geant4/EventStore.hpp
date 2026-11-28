// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Utilities/detail/ScopeExit.hpp"
#include "ActsExamples/EventData/PropagationSummary.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticleDecay.hpp"

#include <exception>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "G4Types.hh"

namespace ActsExamples::Geant4 {

/// Owned collections and diagnostics returned by a transport event.
struct EventResult {
  /// Detect duplicate particle identities during transport.
  using ParticleContainer = std::set<SimParticle, detail::CompareParticleId>;

  /// Final particle states that pass the configured output selection.
  ParticleContainer particlesSimulated;
  /// Hits in transport order, before the framework sorts them by sensor.
  SimHitContainer::sequence_type hits;
  /// Per-track material recorded by the material stepping action.
  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack> materialTracks;
  /// Optional propagation summaries by Geant4 track ID.
  std::unordered_map<G4int, PropagationSummary> propagationRecords;

  /// Number of steps contributing to sensitive hits across the event.
  std::size_t numberGeantSteps = 0;
  /// Largest number of merged Geant4 steps in one hit.
  std::size_t maxStepsForHit = 0;
  /// Counts of truth bookkeeping failures for reporting after transport.
  std::size_t particleIdCollisionsInitial = 0;
  std::size_t particleIdCollisionsFinal = 0;
  std::size_t parentIdNotFound = 0;
};

/// Inputs, outputs and scratch storage for one synchronous transport call.
struct EventStore {
  /// @param geometry Alignment context for the duration of transport.
  /// @param particles Read-only input for the duration of transport.
  /// @param graph Optional generator-prescribed decays for the selected roots.
  EventStore(const Acts::GeometryContext& geometry,
             const SimParticleContainer& particles,
             const SimParticleDecayGraph* graph = nullptr)
      : geoContext(geometry), inputParticles(particles), decayGraph(graph) {}

  /// Only the user actions see alignment here. The Geant4 geometry stays fixed.
  const Acts::GeometryContext& geoContext;
  /// Flat collection to translate to Geant4 primaries.
  const SimParticleContainer& inputParticles;
  /// The event owns the optional graph for the duration of transport.
  const SimParticleDecayGraph* decayGraph;
  /// A rejected primary input leaves an empty G4 event. Report its error after
  /// Geant4 has taken ownership and completed event cleanup.
  std::exception_ptr primaryError;
  /// The result moves to the caller when transport completes.
  EventResult result;

  /// Initial states used by the tracking action to assemble output particles.
  EventResult::ParticleContainer particlesInitial;

  /// Hit buffer for step merging (multiple steps in sensitive volume)
  std::vector<ActsFatras::Hit> hitBuffer;

  /// Particle hit count (for hit indexing)
  std::unordered_map<SimBarcode, std::size_t> particleHitCount;
  /// Particle status
  std::unordered_map<SimBarcode, ActsFatras::SimulationOutcome> particleOutcome;
  /// Geant4 Track ID to Barcode mapping
  std::unordered_map<G4int, SimBarcode> trackIdMapping;
  /// Geant4 Track ID subparticle counter (for subparticle indexing)
  std::unordered_map<G4int, std::size_t> trackIdSubparticleCount;

  /// Store subparticle count for {primVertex, secVertex, part, gen}
  /// This is done using a pseudo-barcode that contains all fields but not the
  /// subparticle counter. This can be used as key in a map to store the
  /// subparticle information
  std::unordered_map<SimBarcode, std::size_t> subparticleMap;
};

/// Connects run-long Geant4 actions to the event being transported.
/// Access is confined to callbacks on the session thread.
class EventContext {
 public:
  EventContext() = default;
  EventContext(const EventContext&) = delete;
  EventContext& operator=(const EventContext&) = delete;

  /// Bind an event for the lifetime of the returned guard.
  /// @param event Storage that outlives the guard.
  /// @return Scope guard that clears the active event on destruction.
  [[nodiscard]] auto bind(EventStore& event) {
    if (m_event != nullptr) {
      throw std::logic_error("A Geant4 event is already bound");
    }
    m_event = &event;
    return Acts::detail::ScopeExit([this] { m_event = nullptr; });
  }

  /// Access event state while Geant4 is processing an event.
  EventStore& current() { return activeEvent(); }
  /// Read event state while Geant4 is processing an event.
  const EventStore& current() const { return activeEvent(); }

 private:
  EventStore& activeEvent() const {
    if (m_event == nullptr) {
      throw std::logic_error("No active Geant4 event");
    }
    return *m_event;
  }

  EventStore* m_event = nullptr;
};

}  // namespace ActsExamples::Geant4
