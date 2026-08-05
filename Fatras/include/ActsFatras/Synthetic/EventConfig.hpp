// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Steering for a synthetic pile-up event. See @ref fatras_synthetic_events.

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"

#include <cstddef>
#include <cstdint>

namespace ActsFatras::Synthetic {

/// Which point of the energy loss distribution a surface crossing costs.
enum class EnergyLossModel {
  /// The mean loss, `Acts::computeEnergyLossMean`
  Mean,
  /// The most probable loss, `Acts::computeEnergyLossMode`
  Mode,
};

/// What the collisions make: how many charged particles and with which
/// kinematics. It knows nothing of the detector they are then walked through.
struct GenerationConfig {
  /// Number of minimum-bias interactions overlaid in the event
  std::size_t pileup = 200;
  /// Charged particles per interaction and unit of rapidity above `minPt`,
  /// averaged over the generated range rather than central
  float chargedPerUnitRapidity = 5.12f;

  /// Lowest generated transverse momentum; below it a particle is a stub
  float minPt = 0.02 * Acts::UnitConstants::GeV;
  /// Momentum scale of the spectrum, see `detail::samplePt`
  float ptScale = 1.338 * Acts::UnitConstants::GeV;
  /// Falloff exponent of the spectrum, see `detail::samplePt`
  float ptExponent = 4.89f;

  /// Lower end of the generated *rapidity* range: not pseudorapidity, because
  /// `dN/deta = (p/E) dN/dy` below the pion mass.
  ///
  /// @note Wider than the presets are validated over, so that the tracks just
  ///       past it leave their clusters.
  float minRapidity = -4.3f;
  /// @copydoc minRapidity
  float maxRapidity = 4.3f;
  /// Rapidity at which the plateau has fallen to half, the shape being a Fermi
  /// factor; see `detail::sampleRapidity`.
  float rapidityEdge = 4.3f;
  /// Width of that edge; zero leaves the plateau flat over the whole range.
  /// @copydoc rapidityEdge
  float rapidityEdgeWidth = 0.f;

  /// Longitudinal spread of the luminous region
  float beamspotSigmaZ = 50 * Acts::UnitConstants::mm;
  /// Transverse impact parameter spread of the primaries.
  /// @note The core only; the tail a full simulation calls primary is heavy
  ///       flavour, which this generator makes as secondaries.
  float d0Sigma = 11.5 * Acts::UnitConstants::um;

  /// Number of primary particles the pile-up settings amount to.
  /// @return the number of primaries to generate
  std::size_t numPrimaries() const {
    return static_cast<std::size_t>(static_cast<float>(pileup) *
                                    chargedPerUnitRapidity *
                                    (maxRapidity - minRapidity));
  }
};

/// How far a track is walked, and how much room the walk is given.
struct PropagationConfig {
  /// Turning angle to propagate a track through, in turns; past half a turn a
  /// soft one curls back through the layers it has already crossed
  float maxTurns = 0.5f;

  /// @name Storage
  ///
  /// What one primary is expected to amount to, i.e. how much room an event is
  /// given before it is generated: too little costs a reallocation, too much
  /// memory.
  ///
  /// @{

  /// Tracks a primary leads to, itself and the secondaries it makes
  float tracksPerPrimary = 2.f;
  /// Space points all of those leave
  float hitsPerPrimary = 12.f;

  /// @}

  /// Ceiling on the tracks one event may hold, as a multiple of the primaries;
  /// past it no further secondaries are produced. Two orders above
  /// `tracksPerPrimary`, being a guard rather than an expectation.
  float maxTracksPerPrimary = 100.f;
};

/// What the material a track crosses does to it. Both effects displace the
/// *measured* position of a hit and leave the trajectory alone: a scatter
/// slides a hit by the angle times the lever arm, energy loss by half the
/// curvature change times that lever arm squared.
struct MaterialConfig {
  /// @name Path length
  ///
  /// Path through a surface, in units of its normal thickness, above which the
  /// material it is worth stops growing: the aspect ratio of a module.
  ///
  /// @{

  /// For a disc, bounded by the radial size of a module
  float maxDiscPathLength = 4.f;
  /// For a cylinder, effectively unbounded, a beam pipe being no module
  float maxCylinderPathLength = 100.f;

  /// @}

  /// Multiplier on the material the layout carries; zero turns it off
  float scale = 1.f;
  /// Whether a track is deflected by the material it crosses
  bool multipleScattering = true;
  /// Whether a track loses energy to it
  bool energyLoss = true;
  /// Whether a crossing costs the mean loss or the most probable one.
  /// @note The mode by default; the delta rays the mean includes are made as
  ///       secondaries here, so it would count them twice.
  EnergyLossModel energyLossModel = EnergyLossModel::Mode;
  /// Fractional momentum loss past which the *displacement* stops growing,
  /// being first order in it. The momentum itself keeps falling.
  float maxEnergyLossFraction = 0.5f;
};

/// What a sensitive crossing is read out as.
struct MeasurementConfig {
  /// Single space point resolution, along the two directions a sensor measures:
  /// `r*phi` and z on a cylinder, `r*phi` and r on a disc
  float positionSmearing = 15 * Acts::UnitConstants::um;
  /// Multiplier on the module overlap the layout carries; zero turns it off.
  /// See `DetectorSurface::overlapProbability`.
  float overlapScale = 1.f;
};

/// Where one secondary goes and how hard it is. The momentum along the parent
/// and the kick across it are drawn independently, the opening angle being
/// their ratio. Three channels: a converted or knocked-out **electron**, a
/// forward **cascade** product and an isotropic **evaporation** product. How
/// many are made is `SecondaryConfig`'s.
struct SecondarySamplingConfig {
  /// Transverse momentum below which a secondary becomes a *stub* rather than
  /// leaving the surface it was made on; see `SecondaryConfig::stubClusters`.
  float minPt = 5 * Acts::UnitConstants::MeV;

  /// Median momentum of such an electron, at a parent of one GeV
  float electronScale = 0.137 * Acts::UnitConstants::GeV;
  /// Power of the parent momentum its median follows
  float electronExponent = 0.155f;
  /// Width of the log-normal about that median, in e-folds
  float electronSpread = 2.80f;
  /// Scale of its transverse kick, twenty times below the nuclear one.
  ///
  /// @note Always about the *radial* axis, a collinear daughter of a charged
  ///       parent otherwise overflowing the innermost decade of |d0|.
  float electronKt = 0.014 * Acts::UnitConstants::GeV;

  /// Median momentum of a cascade product *along its parent*, at a parent of
  /// one GeV
  float momentumScale = 0.448 * Acts::UnitConstants::GeV;
  /// Power of the parent momentum that median follows
  float momentumExponent = 0.429f;
  /// Width of the log-normal about that median, in e-folds.
  float momentumSpread = 1.67f;
  /// Rayleigh scale of a cascade product's transverse kick, i.e. the limited
  /// transverse momentum of soft hadroproduction.
  ///
  /// @note The *scale*, neither the median nor the mean.
  float kt = 0.319 * Acts::UnitConstants::GeV;

  /// Share of the nuclear channel emitted isotropically, i.e. the struck
  /// nucleus breaking up
  float evaporationFraction = 0.30f;
  /// Rayleigh scale of its momentum, set by the Fermi momentum
  float evaporationScale = 0.30 * Acts::UnitConstants::GeV;
};

/// How much a crossing produces, and how deep the cascade it starts may go.
struct SecondaryConfig {
  /// Secondaries a crossing yields per radiation length, the electron channel.
  /// Effective rather than an interaction probability, so it belongs to a
  /// layout as much as to a detector.
  float electronRate = 4.941f;
  /// The same per nuclear interaction length, the nuclear channel: two rates
  /// and not one, the channels counting in different lengths.
  float nuclearRate = 11.545f;
  /// Mean number of K0S and Lambda per primary, the only source of secondaries
  /// produced away from a surface
  float decayYield = 0.261f;
  /// Mean flight length of such a parent *at a GeV*; it scales with momentum
  float decayLength = 60 * Acts::UnitConstants::mm;

  /// Mean number of stubs a crossing yields per radiation length, on top of the
  /// daughters that fall below `SecondarySamplingConfig::minPt`; its own rate,
  /// being a different process rather than the soft tail of that one.
  float stubRate = 0.f;
  /// Mean number of space points a stub leaves, which is at least one.
  float stubClusters = 2.1f;
  /// How far a stub's space points sit from where it was made, i.e. the curl
  /// radius of a particle this soft
  float stubReach = 4 * Acts::UnitConstants::mm;

  /// How many interactions deep a cascade may go: one lets only a primary
  /// produce, two lets its daughters produce in turn
  std::uint32_t maxGenerations = 1;
  /// Ceiling on the mean number one crossing may yield, a guard against
  /// exhausting memory
  float maxPerCrossing = 20.f;

  /// What one of them comes out as
  SecondarySamplingConfig sampling;
};

/// What the detector does to a track: how far it is followed, what the material
/// does to it, what is measured of it, and what it makes on the way.
struct SimulationConfig {
  /// How far a track is followed
  PropagationConfig propagation;
  /// What it crosses on the way
  MaterialConfig material;
  /// What a sensitive crossing is read out as
  MeasurementConfig measurement;
  /// What a crossing makes
  SecondaryConfig secondaries;
};

/// Steering for the synthetic event: what is made, what the detector does to
/// it, and the few settings both halves share. The detector enters several of
/// these, so prefer a preset over the defaults.
struct EventConfig {
  /// The particles of the event, before any detector
  GenerationConfig generation;
  /// What that detector then does to them
  SimulationConfig simulation;

  /// The species every track is taken to be. Both halves need it: a rapidity is
  /// an angle only given a mass, and so is an energy loss.
  Acts::PdgParticle particlePdg = Acts::PdgParticle::ePionPlus;
  /// Solenoid field along z
  float bFieldZ = 2 * Acts::UnitConstants::T;
  /// Random seed, so that events are reproducible across runs and seeders
  std::uint32_t seed = 12345;

  /// Mass of `particlePdg`.
  /// @return the mass
  /// @throws std::invalid_argument if the core library has no mass for it
  float particleMass() const;

  /// @name Tuned configurations
  ///
  /// One per layout, fitted against a full simulation of that detector. Each
  /// spells out every field, so that a change to a default cannot retune one.
  ///
  /// @{

  /// Configuration for an ITk-like ttbar event at a pile-up of 200.
  /// @return the configuration
  static EventConfig itkPixelTtbarPu200();

  /// Configuration for an Open Data Detector ttbar event at a pile-up of 200.
  /// @return the configuration
  static EventConfig openDataDetectorTtbarPu200();

  /// @}
};

}  // namespace ActsFatras::Synthetic
