// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/EventConfig.hpp"

#include "Acts/Definitions/ParticleData.hpp"

#include <optional>
#include <stdexcept>

using namespace Acts::UnitLiterals;

namespace ActsFatras::Synthetic {

float EventConfig::particleMass() const {
  const std::optional<float> mass =
      Acts::findMass(Acts::makeAbsolutePdgParticle(particlePdg));
  if (!mass.has_value()) {
    throw std::invalid_argument(
        "EventConfig: the core library has no mass for this particlePdg");
  }
  return *mass;
}

EventConfig EventConfig::itkPixelTtbarPu200() {
  EventConfig config;
  GenerationConfig& generation = config.generation;
  generation.pileup = 200;
  // Fitted against 250 events of a GNN4ITk ttbar dump. 8.9k primaries per event
  // leave an ITk pixel cluster inside the comparison's acceptance; the
  // generator makes more, reaching below it
  generation.chargedPerUnitRapidity = 6.15f;
  generation.minPt = 0.02_GeV;
  generation.ptScale = 0.672f;
  generation.ptExponent = 5.44f;
  generation.minRapidity = -4.3f;
  generation.maxRapidity = 4.3f;
  generation.rapidityEdge = 5.30f;
  generation.rapidityEdgeWidth = 1.30f;
  generation.beamspotSigmaZ = 50.f;
  generation.d0Sigma = 0.0120f;

  SimulationConfig& simulation = config.simulation;
  MaterialConfig& material = simulation.material;
  material.maxDiscPathLength = 4.00f;
  material.maxCylinderPathLength = 100.00f;
  material.scale = 1.000f;
  material.multipleScattering = true;
  material.energyLoss = true;
  material.energyLossModel = EnergyLossModel::Mode;
  material.maxEnergyLossFraction = 0.5f;
  // Only a curler gets this far. The reference's soft primaries put five
  // clusters on their busiest layer, which is two and a half turns.
  simulation.propagation.maxTurns = 3.00f;
  // 50 um pitch read out digitally, so pitch / sqrt(12)
  simulation.measurement.positionSmearing = 15_um;

  SecondaryConfig& secondaries = simulation.secondaries;
  secondaries.electronRate = 3.202f;
  secondaries.nuclearRate = 7.483f;
  secondaries.decayYield = 0.185f;
  secondaries.decayLength = 60_mm;
  // taken from the ODD and scaled to this layout's rate: a reference that
  // records a secondary only above a few hundred MeV cannot measure stubs
  secondaries.stubRate = 1.267f;
  secondaries.stubClusters = 2.1f;
  secondaries.stubReach = 4_mm;

  SecondarySamplingConfig& sampling = secondaries.sampling;
  sampling.minPt = 5_MeV;
  // measured on ColliderML, which links every cluster and so needs no
  // threshold correction; see the ODD preset below
  sampling.electronScale = 0.137f;
  sampling.electronExponent = 0.155f;
  sampling.electronSpread = 2.800f;
  sampling.electronKt = 0.014f;
  // the longitudinal law, and the kick beside it rather than out of it; both
  // measured on the dump's parent links over the hadronic channel alone
  sampling.momentumScale = 0.448f;
  sampling.momentumExponent = 0.429f;
  sampling.momentumSpread = 1.670f;
  sampling.kt = 0.319f;
  sampling.evaporationFraction = 0.300f;
  sampling.evaporationScale = 0.300f;

  config.particlePdg = Acts::PdgParticle::ePionPlus;
  config.bFieldZ = 2_T;
  config.seed = 12345;
  return config;
}

EventConfig EventConfig::openDataDetectorTtbarPu200() {
  EventConfig config;
  GenerationConfig& generation = config.generation;
  generation.pileup = 200;
  // 8.1k primaries per event leave an ODD pixel cluster inside the
  // comparison's acceptance; the generator makes more, reaching below it
  generation.chargedPerUnitRapidity = 5.58f;
  generation.minPt = 0.02_GeV;
  generation.ptScale = 0.663f;
  generation.ptExponent = 5.30f;
  generation.minRapidity = -4.3f;
  generation.maxRapidity = 4.3f;
  generation.rapidityEdge = 4.90f;
  generation.rapidityEdgeWidth = 1.30f;
  generation.beamspotSigmaZ = 57.f;
  generation.d0Sigma = 0.0142f;

  SimulationConfig& simulation = config.simulation;
  MaterialConfig& material = simulation.material;
  material.maxDiscPathLength = 4.00f;
  material.maxCylinderPathLength = 100.00f;
  material.scale = 1.000f;
  material.multipleScattering = true;
  material.energyLoss = true;
  material.energyLossModel = EnergyLossModel::Mode;
  material.maxEnergyLossFraction = 0.5f;
  // @copydoc itkPixelTtbarPu200
  simulation.propagation.maxTurns = 3.00f;
  // core of the cluster-to-truth residual at normal incidence, below the 14 um
  // of a digital 50 um pitch because clusters are charge-weighted
  simulation.measurement.positionSmearing = 8_um;

  SecondaryConfig& secondaries = simulation.secondaries;
  secondaries.electronRate = 2.744f;
  secondaries.nuclearRate = 6.410f;
  // an order below the ITk's: ColliderML counts far fewer of its secondaries as
  // produced inside the beam pipe than the Athena barcode convention does, and
  // this is the whole of the missing |d0| > 100 mm tail
  secondaries.decayYield = 0.015f;
  secondaries.decayLength = 60_mm;
  secondaries.stubRate = 1.800f;
  secondaries.stubClusters = 2.1f;
  secondaries.stubReach = 4_mm;

  SecondarySamplingConfig& sampling = secondaries.sampling;
  sampling.minPt = 5_MeV;
  // where the law was measured. How an interaction shares momentum out does not
  // depend on which detector is watching, so the ITk carries the same numbers.
  sampling.electronScale = 0.137f;
  sampling.electronExponent = 0.155f;
  sampling.electronSpread = 2.800f;
  sampling.electronKt = 0.014f;
  // the longitudinal law, and the kick beside it rather than out of it; both
  // measured on the dump's parent links over the hadronic channel alone
  sampling.momentumScale = 0.448f;
  sampling.momentumExponent = 0.429f;
  sampling.momentumSpread = 1.670f;
  sampling.kt = 0.319f;
  sampling.evaporationFraction = 0.300f;
  sampling.evaporationScale = 0.300f;

  config.particlePdg = Acts::PdgParticle::ePionPlus;
  config.bFieldZ = 2_T;
  config.seed = 12345;
  return config;
}

}  // namespace ActsFatras::Synthetic
