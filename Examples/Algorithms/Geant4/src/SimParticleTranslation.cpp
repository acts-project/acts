// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/SimParticleTranslation.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/UnitConversion.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

#include <G4ChargedGeantino.hh>
#include <G4Decay.hh>
#include <G4DecayProducts.hh>
#include <G4DynamicParticle.hh>
#include <G4Event.hh>
#include <G4Geantino.hh>
#include <G4LogicalVolume.hh>
#include <G4Navigator.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessType.hh>
#include <G4Track.hh>
#include <G4TransportationManager.hh>
#include <G4VSolid.hh>
#include <G4VUserPrimaryParticleInformation.hh>

namespace ActsExamples::Geant4 {

namespace {

class ParticleSource final : public G4VUserPrimaryParticleInformation {
 public:
  explicit ParticleSource(const SimParticle& particle) : source(particle) {}
  void Print() const override {}
  const SimParticle& source;
};

// Match G4PrimaryTransformer::IsGoodForTrack. Otherwise Geant4 silently
// promotes a short-lived particle's daughters to its own parent vertex.
bool canTrack(const G4ParticleDefinition* definition) {
  return definition != nullptr && (!definition->IsShortLived() ||
                                   definition->GetDecayTable() != nullptr);
}

bool hasActiveDecay(const G4ParticleDefinition& definition) {
  auto* manager = definition.GetProcessManager();
  if (!canTrack(&definition) || definition.GetPDGStable() ||
      manager == nullptr) {
    return false;
  }
  for (int i = 0; i < manager->GetProcessListLength(); ++i) {
    auto* decay = dynamic_cast<G4Decay*>((*manager->GetProcessList())[i]);
    if (decay != nullptr && decay->IsApplicable(definition) &&
        manager->GetProcessActivation(i)) {
      return true;
    }
  }
  return false;
}

void attachDecay(const SimParticle& parent, G4PrimaryParticle& primary,
                 EventStore& store, std::unordered_set<SimBarcode>& attached) {
  const auto it = store.decayGraph->find(parent.particleId());
  if (it == store.decayGraph->end()) {
    return;
  }
  const auto& decay = it->second;
  if (!hasActiveDecay(*primary.GetG4code()) ||
      !std::isfinite(decay.properTime) || decay.properTime < 0. ||
      decay.daughters.size() < 2) {
    throw std::invalid_argument(
        "Invalid prescribed decay or missing active Geant4 decay for PDG " +
        std::to_string(parent.pdg()));
  }
  primary.SetProperTime(decay.properTime * convertTimeToGeant4);
  for (const SimParticle& child : decay.daughters) {
    const auto id = child.particleId();
    if (id.generation() == 0 || store.inputParticles.contains(id) ||
        !attached.insert(id).second ||
        child.parentParticleId() != parent.particleId() ||
        id != parent.particleId().makeDescendant(id.subParticle())) {
      throw std::invalid_argument(
          "Prescribed daughter has overlapping or inconsistent ancestry");
    }
    // Transport-generated children must start beyond the assigned generator
    // IDs.
    auto& last = store.subparticleMap[id.withoutSubparticle()];
    last = std::max(last, static_cast<std::size_t>(id.subParticle()));
    auto* def = G4ParticleTable::GetParticleTable()->FindParticle(child.pdg());
    if (!canTrack(def)) {
      throw std::invalid_argument(
          "Geant4 cannot track prescribed daughter PDG " +
          std::to_string(child.pdg()));
    }
    auto daughter = std::make_unique<G4PrimaryParticle>(def);
    daughter->SetMass(child.mass() * convertEnergyToGeant4);
    daughter->SetCharge(child.charge());
    const Acts::Vector3 mom = child.momentum() * convertEnergyToGeant4;
    daughter->SetMomentum(mom[0], mom[1], mom[2]);
    daughter->SetUserInformation(new ParticleSource(child));
    attachDecay(child, *daughter, store, attached);
    primary.SetDaughter(daughter.release());
  }
}

}  // namespace

const SimParticle* SimParticleTranslation::sourceParticle(
    const G4Track& track) {
  const auto* primary = track.GetDynamicParticle()->GetPrimaryParticle();
  if (primary == nullptr ||
      (track.GetParentID() != 0 &&
       (track.GetCreatorProcess() == nullptr ||
        track.GetCreatorProcess()->GetProcessType() != fDecay))) {
    return nullptr;
  }
  const auto* info =
      dynamic_cast<const ParticleSource*>(primary->GetUserInformation());
  return info == nullptr ? nullptr : &info->source;
}

void SimParticleTranslation::prepareDecay(const G4Track& track) {
  const auto* dynamic = track.GetDynamicParticle();
  const auto* primary = dynamic->GetPrimaryParticle();
  if (primary == nullptr || primary->GetDaughter() == nullptr) {
    return;
  }
  if (const auto* existing = dynamic->GetPreAssignedDecayProducts()) {
    // G4PrimaryTransformer derives each daughter's kinetic energy from the
    // Geant4 table mass, then applies its generator mass. Reapply the lab
    // momentum after that assignment to preserve the input four-vector.
    for (G4int i = 0; i < existing->entries(); ++i) {
      auto* child = (*existing)[i];
      child->SetMomentum(child->GetPrimaryParticle()->GetMomentum());
    }
    return;
  }
  // G4Decay copies its daughters, keeping their primary pointer but dropping
  // their own prescribed products and lifetime. Use the original generator
  // parent as the reference for G4DecayProducts::Boost at the actual decay.
  G4DynamicParticle reference(primary->GetG4code(), primary->GetMomentum());
  reference.SetMass(primary->GetMass());
  reference.SetMomentum(primary->GetMomentum());
  auto products = std::make_unique<G4DecayProducts>(reference);
  for (auto* child = primary->GetDaughter(); child != nullptr;
       child = child->GetNext()) {
    auto item = std::make_unique<G4DynamicParticle>(child->GetG4code(),
                                                    child->GetMomentum());
    item->SetMass(child->GetMass());
    item->SetMomentum(child->GetMomentum());
    item->SetCharge(child->GetCharge());
    item->SetPrimaryParticle(child);
    products->PushProducts(item.release());
  }
  // G4Track exposes only a const view of its owned, mutable dynamic particle.
  auto& particle = *const_cast<G4DynamicParticle*>(dynamic);
  particle.SetPreAssignedDecayProperTime(primary->GetProperTime());
  particle.SetPreAssignedDecayProducts(products.release());
}

SimParticleTranslation::SimParticleTranslation(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4VUserPrimaryGeneratorAction(),
      m_cfg(cfg),
      m_logger(std::move(logger)) {}

void SimParticleTranslation::GeneratePrimaries(G4Event* eventPtr) {
  assert(eventPtr != nullptr);
  G4Event& event = *eventPtr;

  if (!m_cfg.preserveGeant4EventId) {
    event.SetEventID(m_eventNr++);
  }
  std::uint32_t eventID = event.GetEventID();

  ACTS_DEBUG("Primary Generator Action for Event: " << eventID);

  // G4RunManager::GenerateEvent still owns only a raw local event here.
  // Unwinding through it leaks the event and terminates the serial run.
  try {
    addPrimaries(event);
  } catch (...) {
    eventStore().primaryError = std::current_exception();
  }
}

void SimParticleTranslation::addPrimaries(G4Event& event) {
  const auto& inputParticles = eventStore().inputParticles;

  // Reserve hopefully enough hit space
  eventStore().result.hits.reserve(inputParticles.size() *
                                   m_cfg.reserveHitsPerParticle);

  if (inputParticles.empty()) {
    return;
  }
  auto* world = G4TransportationManager::GetTransportationManager()
                    ->GetNavigatorForTracking()
                    ->GetWorldVolume();
  const auto* solid = world->GetLogicalVolume()->GetSolid();
  EInside region = kOutside;
  // Default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4PrimaryVertex* pVertex = nullptr;
  std::unordered_set<SimBarcode> attached;
  std::vector<std::unique_ptr<G4PrimaryVertex>> pending;

  // We are looping through the particles and flush per vertex
  std::optional<Acts::Vector4> lastVertex;

  G4int nextTrackId = 1;
  // Group adjacent roots by their production vertex.
  for (const auto& part : inputParticles) {
    const Acts::Vector4 currentVertex = part.fourPosition();
    // Add a new primary to the vertex

    const Acts::Vector4 mom4 = part.fourMomentum() * convertEnergyToGeant4;

    // Particle properties, may be forced to specific value
    const G4int particlePdgCode = m_cfg.forcedPdgCode.value_or(part.pdg());
    const G4double particleCharge = m_cfg.forcedCharge.value_or(part.charge());
    const G4double particleMass =
        m_cfg.forcedMass.value_or(part.mass() * convertEnergyToGeant4);

    // Check if it is a Geantino / ChargedGeantino
    const G4ParticleDefinition* particleDefinition =
        particleTable->FindParticle(particlePdgCode);
    if (particleDefinition == nullptr) {
      if (particlePdgCode == 0 && particleMass == 0 && particleCharge == 0) {
        particleDefinition = G4Geantino::Definition();
      }
      if (particlePdgCode == 0 && particleMass == 0 && particleCharge != 0) {
        if (particleCharge != 1) {
          ACTS_ERROR("invalid charged geantino charge " << particleCharge
                                                        << ". should be 1");
        }
        particleDefinition = G4ChargedGeantino::Definition();
      }
    }

    // Skip if translation failed
    if (particleDefinition == nullptr) {
      if (eventStore().decayGraph != nullptr &&
          eventStore().decayGraph->contains(part.particleId())) {
        throw std::invalid_argument(
            "Geant4 cannot track prescribed parent PDG " +
            std::to_string(part.pdg()));
      }
      ACTS_DEBUG(
          "Could not translate particle with PDG code : " << particlePdgCode);
      continue;
    }

    const bool newVertex = !lastVertex || !currentVertex.isApprox(*lastVertex);
    if (newVertex) {
      lastVertex = currentVertex;
      const G4ThreeVector position(currentVertex[0] * convertLengthToGeant4,
                                   currentVertex[1] * convertLengthToGeant4,
                                   currentVertex[2] * convertLengthToGeant4);
      region = currentVertex.allFinite() ? solid->Inside(position) : kOutside;
      auto vertex = std::make_unique<G4PrimaryVertex>(
          position, currentVertex[3] * convertTimeToGeant4);
      pVertex = vertex.get();
      pending.push_back(std::move(vertex));
    }
    // Geant4 dereferences the starting touchable in SetInitialStep. An
    // outside root has none. A root exactly on the boundary can enter or
    // leave the world depending on its direction, so consult a private
    // navigator without altering the tracking navigator's state.
    bool outside = region == kOutside;
    if (region == kSurface) {
      G4Navigator navigator;
      navigator.SetWorldVolume(world);
      const auto& dir = part.direction();
      const G4ThreeVector direction(dir[0], dir[1], dir[2]);
      outside =
          navigator.LocateGlobalPointAndSetup(
              pVertex->GetPosition(), &direction, false, false) == nullptr;
    }
    if (outside) {
      std::ostringstream message;
      message << "Geant4 primary " << part.particleId() << " (PDG "
              << part.pdg() << ") starts outside the world at "
              << pVertex->GetPosition()
              << " mm. Select input roots within the detector world.";
      throw std::invalid_argument(message.str());
    }

    ACTS_VERBOSE("Adding particle with name '"
                 << particleDefinition->GetParticleName()
                 << "' and properties:");
    ACTS_VERBOSE(" -> mass: " << particleMass);
    ACTS_VERBOSE(" -> charge: " << particleCharge);
    ACTS_VERBOSE(" -> momentum: " << mom4.transpose());

    auto particle = std::make_unique<G4PrimaryParticle>(particleDefinition);

    particle->SetMass(particleMass);
    particle->SetCharge(particleCharge);
    if (eventStore().decayGraph != nullptr) {
      // Preserve the specified mass, including for highly boosted parents.
      particle->SetMomentum(mom4[0], mom4[1], mom4[2]);
      particle->SetUserInformation(new ParticleSource(part));
      attachDecay(part, *particle, eventStore(), attached);
    } else {
      particle->Set4Momentum(mom4[0], mom4[1], mom4[2], mom4[3]);
      // The flat-input path also supplies mappings to material-only actions.
      particle->SetTrackID(nextTrackId++);
      eventStore().particlesInitial.insert(part);
      eventStore().trackIdMapping.emplace(particle->GetTrackID(),
                                          part.particleId());
    }
    pVertex->SetPrimary(particle.release());
  }
  // Only commit primaries once every root and its optional decay tree is valid.
  for (auto& vertex : pending) {
    event.AddPrimaryVertex(vertex.release());
  }
}

}  // namespace ActsExamples::Geant4
