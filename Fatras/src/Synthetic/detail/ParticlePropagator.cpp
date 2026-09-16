// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/detail/ParticlePropagator.hpp"

#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"
#include "ActsFatras/Synthetic/detail/Propagation.hpp"
#include "ActsFatras/Synthetic/detail/Sampling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

using namespace Acts::UnitLiterals;

namespace ActsFatras::Synthetic::detail {

namespace {

constexpr float pi = std::numbers::pi_v<float>;

/// Turn a count given per primary into one for the whole event.
std::size_t totalForEvent(const float perPrimary,
                          const std::size_t numPrimaries) {
  return static_cast<std::size_t>(std::max(0.f, perPrimary) *
                                  static_cast<float>(numPrimaries));
}

/// A direction, split the way the scattering frame and the projections onto a
/// surface both want it.
struct TrackDirection {
  /// Path along the track per unit of transverse path, i.e. cosh(eta)
  float invSinTheta{};
  /// Sine of the polar angle
  float sinTheta{};
  /// Cosine of the polar angle
  float cosTheta{};
};

TrackDirection trackDirection(const float cotTheta) {
  const float invSinTheta = std::hypot(1.f, cotTheta);
  const float sinTheta = 1.f / invSinTheta;
  return TrackDirection{invSinTheta, sinTheta, cotTheta * sinTheta};
}

/// What a crossing of a surface is worth, its material in the length each
/// effect counts in.
struct CrossingMaterial {
  /// Multiplier from the surface's normal thickness to the path taken through
  /// it
  float scale{};
  /// Radiation lengths crossed
  float xOverX0{};
  /// Mean number of secondaries produced
  float rate{};
  /// Share of those that are electrons rather than nuclear products
  float electronShare{};
};

CrossingMaterial crossingMaterial(const SimulationConfig& cfg,
                                  const DetectorSurface& surface,
                                  const SurfaceCrossing& crossing) {
  const float maxPath = surface.shape == SurfaceShape::Cylinder
                            ? cfg.material.maxCylinderPathLength
                            : cfg.material.maxDiscPathLength;
  // the band the crossing landed in, which is how a ring is told from the gap
  // beside it and the near end of a service cylinder from its far end
  const Acts::MaterialSlab& slab = *crossing.material;

  CrossingMaterial material;
  material.scale =
      cfg.material.scale * std::clamp(crossing.pathLength, 1.f, maxPath);
  material.xOverX0 = slab.thicknessInX0() * material.scale;
  // Each channel in the length it counts in, so the surface's composition
  // decides the mix and nothing is fitted on top.
  const float electronPart = cfg.secondaries.electronRate * material.xOverX0;
  const float nuclearPart =
      cfg.secondaries.nuclearRate * slab.thicknessInL0() * material.scale;
  material.rate = electronPart + nuclearPart;
  material.electronShare =
      material.rate > 0.f ? electronPart / material.rate : 1.f;
  return material;
}

/// A secondary that goes nowhere, left where the crossing made it.
PendingStub makeStub(std::mt19937& rng, const SecondarySamplingConfig& cfg,
                     const SurfaceCrossing& crossing,
                     const std::uint8_t generation) {
  const float pt = cfg.minPt * sampleUniform01(rng);
  return PendingStub{
      crossing.r,        crossing.phi,   crossing.z,       pt,
      sampleCharge(rng), crossing.layer, crossing.surface, generation};
}

/// Record what a sensitive crossing measures: one space point where the track
/// met the module, and a second where an overlapping module measured it too.
/// `displaceU` and `displaceV` are how far the material crossed so far has
/// moved the track across itself, in the transverse plane and out of it.
/// Returns how many space points were appended to `hits`.
std::uint32_t recordHits(std::mt19937& rng, const DetectorLayout& layout,
                         const MeasurementConfig& cfg, const Helix& helix,
                         const TrackDirection& direction,
                         const SurfaceCrossing& crossing, const float displaceU,
                         const float displaceV, const std::uint32_t particle,
                         std::vector<SmearedHit>& hits) {
  const DetectorSurface& surface = layout.surfaces[crossing.surface];
  const bool cylinder =
      layout.layers[crossing.layer].shape == SurfaceShape::Cylinder;
  // Where the track points relative to the surface it is crossing, which both
  // the displacement and the stagger below are projected onto.
  const float deltaPhi =
      helix.phi0 - helix.charge * crossing.gamma - crossing.phi;
  const float cosDeltaPhi = std::cos(deltaPhi);
  const float sinDeltaPhi = std::sin(deltaPhi);

  // smear along the measured directions of the sensor, leaving the normal to
  // it at the nominal surface position
  auto [smearRZ, smearRPhi] = sampleNormalPair(rng);
  smearRZ *= cfg.positionSmearing;
  smearRPhi *= cfg.positionSmearing;
  // kept apart from what the displacement adds below, an overlapping module
  // measuring the same track with an error of its own
  const float firstSmearRZ = smearRZ;
  const float firstSmearRPhi = smearRPhi;

  // Project the displacement onto the surface; the second term is a track
  // displaced across itself meeting the surface elsewhere.
  const float alongRPhi =
      displaceU * cosDeltaPhi - displaceV * direction.cosTheta * sinDeltaPhi;
  const float alongR =
      -displaceU * sinDeltaPhi - displaceV * direction.cosTheta * cosDeltaPhi;
  const float alongZ = displaceV * direction.sinTheta;
  const float slide = (cylinder ? alongR : alongZ) * crossing.pathLength;
  smearRPhi += alongRPhi - slide * direction.sinTheta * sinDeltaPhi;
  smearRZ += cylinder ? alongZ - slide * direction.cosTheta
                      : alongR - slide * direction.sinTheta * cosDeltaPhi;

  std::uint32_t numHits = 0;
  // Displaced, it may land off the module that made it, or in the gap between
  // two rings, where there is support rather than silicon.
  const float hitR = crossing.r + (cylinder ? 0.f : smearRZ);
  const float hitZ = crossing.z + (cylinder ? smearRZ : 0.f);
  if (const std::optional<std::uint32_t> layer =
          findLayer(layout, surface, hitR, hitZ);
      layer.has_value()) {
    hits.push_back(SmearedHit{hitR, crossing.phi + smearRPhi / crossing.r, hitZ,
                              *layer, particle});
    ++numHits;
  }

  // A crossing near the common edge of two modules is measured by both. They
  // alternate in depth, so the second sits a stagger along the normal away at
  // almost the same azimuth.
  const float overlap = surface.overlapProbability * cfg.overlapScale;
  if (overlap <= 0.f || sampleUniform01(rng) >= overlap) {
    return numHits;
  }
  const float dirR = direction.sinTheta * cosDeltaPhi;
  const float dirRPhi = direction.sinTheta * sinDeltaPhi;
  // along the outward normal, which for a disc is away from the interaction
  // point on either side of it
  const float normal =
      cylinder ? dirR : direction.cosTheta * (crossing.z < 0.f ? -1.f : 1.f);
  // a track this close to the tangent does not resolve two modules
  constexpr float minNormal = 0.25f;
  if (std::abs(normal) <= minNormal) {
    return numHits;
  }

  const float step = surface.overlapOffset / normal;
  auto [overlapRZ, overlapRPhi] = sampleNormalPair(rng);
  overlapRZ = smearRZ - firstSmearRZ + overlapRZ * cfg.positionSmearing;
  overlapRPhi = smearRPhi - firstSmearRPhi + overlapRPhi * cfg.positionSmearing;
  const float overlapR =
      crossing.r + step * dirR + (cylinder ? 0.f : overlapRZ);
  const float overlapZ =
      crossing.z + step * direction.cosTheta + (cylinder ? overlapRZ : 0.f);
  if (const std::optional<std::uint32_t> layer =
          findLayer(layout, surface, overlapR, overlapZ);
      layer.has_value()) {
    hits.push_back(SmearedHit{
        overlapR, crossing.phi + (step * dirRPhi + overlapRPhi) / crossing.r,
        overlapZ, *layer, particle});
    ++numHits;
  }
  return numHits;
}

/// Everything a crossing makes: the daughters it sends off, appended to
/// `tracks` and `particles` together, and the stubs it leaves where they were
/// born. `trackP` is what the parent has left at the crossing.
void produceSecondaries(std::mt19937& rng, const SecondaryConfig& cfg,
                        const float bFieldZ, const Helix& helix,
                        const std::uint8_t generation,
                        const SurfaceCrossing& crossing,
                        const CrossingMaterial& material, const float trackP,
                        std::vector<Helix>& tracks,
                        std::vector<GeneratedParticle>& particles,
                        std::vector<PendingStub>& stubs) {
  // Stubs of their own, which only a sensitive surface can record.
  if (crossing.sensitive() && cfg.stubRate > 0.f) {
    const std::uint32_t numStubs = samplePoisson(
        rng, std::min(cfg.stubRate * material.xOverX0, cfg.maxPerCrossing));
    for (std::uint32_t s = 0; s < numStubs; ++s) {
      stubs.push_back(makeStub(rng, cfg.sampling, crossing,
                               static_cast<std::uint8_t>(generation + 1)));
    }
  }

  // Both axes a daughter can come off: the track that crossed, and the line
  // from the beam line a neutral parent came in along.
  const float trackPhi = helix.phi0 - helix.charge * crossing.gamma;
  const float radialCotTheta = (crossing.z - helix.z0) / crossing.r;

  const std::uint32_t numSecondaries =
      samplePoisson(rng, std::min(material.rate, cfg.maxPerCrossing));
  for (std::uint32_t s = 0; s < numSecondaries; ++s) {
    const std::optional<SecondaryKinematics> kinematics = sampleSecondary(
        rng, cfg.sampling, trackPhi, helix.cotTheta, crossing.phi,
        radialCotTheta, trackP, material.electronShare);
    if (!kinematics.has_value()) {
      // Too soft to go anywhere: on a sensitive surface it still leaves the
      // cluster or two it made where it was born, on a passive one nothing.
      if (crossing.sensitive()) {
        stubs.push_back(makeStub(rng, cfg.sampling, crossing,
                                 static_cast<std::uint8_t>(generation + 1)));
      }
      continue;
    }
    // the same species as its parent, which is all a model with one of them
    // can say
    const float charge = sampleCharge(rng);
    const Helix secondary = makeHelixFromPoint(
        crossing.r, crossing.phi, crossing.z, kinematics->phi,
        kinematics->cotTheta, kinematics->pt, charge, bFieldZ);
    tracks.push_back(secondary);

    GeneratedParticle particle = makeParticle(secondary, kinematics->pt);
    particle.productionRadius = crossing.r;
    particle.productionZ = crossing.z;
    particle.productionSurface = crossing.surface;
    particle.generation = static_cast<std::uint8_t>(generation + 1);
    particles.push_back(particle);
  }
}

/// The neutral parents that come with the primaries: a K0S or a Lambda flies
/// straight out of the luminous region along the primary it was made with and
/// decays to charged hadrons. Nothing of the parent is recorded, so only where
/// it decayed is drawn. The primaries are the `numPrimaries` leading entries of
/// `tracks`, and the daughters are appended to it.
void produceDecays(std::mt19937& rng, const SecondaryConfig& cfg,
                   const float bFieldZ, std::vector<Helix>& tracks,
                   std::vector<GeneratedParticle>& particles,
                   const std::size_t numPrimaries) {
  // the primaries are the leading entries of both vectors
  const std::size_t base = particles.size() - numPrimaries;
  for (std::size_t t = 0; t < numPrimaries; ++t) {
    // the vectors grow while daughters are queued, so work on copies
    const Helix primary = tracks[t];
    const TrackDirection direction = trackDirection(primary.cotTheta);
    const float parentP = particles[base + t].pt * direction.invSinTheta;

    const std::uint32_t numDecays = samplePoisson(rng, cfg.decayYield);
    for (std::uint32_t s = 0; s < numDecays; ++s) {
      // a parent flies `beta * gamma * c * tau`, so held at a fixed length a
      // hard V0 decays where a soft one does and the |d0| tail never fills
      const float flight = -cfg.decayLength * (parentP / 1_GeV) *
                           std::log(sampleUniform01Positive(rng));
      const float r = flight * direction.sinTheta;
      const float z = primary.z0 + flight * direction.cosTheta;
      // A K0S or a Lambda decays to charged hadrons, so this is the nuclear
      // channel. It crossed no material, so both axes and the momentum are the
      // undegraded primary's.
      const std::optional<SecondaryKinematics> kinematics =
          sampleSecondary(rng, cfg.sampling, primary.phi0, primary.cotTheta,
                          primary.phi0, primary.cotTheta, parentP,
                          /*electronShare=*/0.f);
      if (!kinematics.has_value()) {
        continue;
      }
      const float charge = sampleCharge(rng);
      const Helix daughter = makeHelixFromPoint(
          std::max(r, 1e-3f), primary.phi0, z, kinematics->phi,
          kinematics->cotTheta, kinematics->pt, charge, bFieldZ);
      tracks.push_back(daughter);

      // a daughter of a primary, and it decayed away from any surface, so
      // the production surface stays at `kNoSurface`
      GeneratedParticle particle = makeParticle(daughter, kinematics->pt);
      particle.productionRadius = r;
      particle.productionZ = z;
      particle.generation = 1;
      particles.push_back(particle);
    }
  }
}

}  // namespace

ParticlePropagator::ParticlePropagator(const DetectorLayout& layout,
                                       const EventConfig& config)
    : m_layout(&layout),
      m_cfg(config),
      m_absPdg(Acts::makeAbsolutePdgParticle(config.particlePdg)),
      m_mass(config.particleMass()) {
  // a track cannot lose more than it carries, and the momentum it keeps
  // divides the loss of every crossing after
  if (config.simulation.material.maxEnergyLossFraction < 0.f ||
      config.simulation.material.maxEnergyLossFraction >= 1.f) {
    throw std::invalid_argument(
        "ParticlePropagator: maxEnergyLossFraction has to be in [0, 1)");
  }
}

void ParticlePropagator::propagate(
    std::mt19937& rng, GeneratorScratch& scratch,
    std::vector<GeneratedParticle>& particles) const {
  const DetectorLayout& layout = *m_layout;
  const SimulationConfig& cfg = m_cfg.simulation;
  const PropagationConfig& propagation = cfg.propagation;
  const SecondaryConfig& secondaries = cfg.secondaries;

  // Secondaries are appended while their parent is walked and picked up by the
  // same loop.
  std::vector<Helix>& tracks = scratch.tracks;
  std::vector<SmearedHit>& hits = scratch.hits;
  std::vector<SurfaceCrossing>& crossings = scratch.crossings;
  std::vector<PendingStub>& stubs = scratch.stubs;
  const std::size_t numPrimaries = tracks.size();
  if (particles.size() < numPrimaries) {
    throw std::invalid_argument(
        "ParticlePropagator: every track needs the particle it is recorded as");
  }
  // Track `t` is `particles[base + t]`, the two growing together from here on.
  const std::size_t base = particles.size() - numPrimaries;

  // room for what the primaries are expected to amount to, so that the first
  // event does not grow the vectors as it walks them
  const std::size_t expectedTracks =
      totalForEvent(propagation.tracksPerPrimary, numPrimaries);
  tracks.reserve(expectedTracks);
  particles.reserve(base + expectedTracks);
  hits.clear();
  hits.reserve(totalForEvent(propagation.hitsPerPrimary, numPrimaries));
  stubs.clear();

  // The neutral parents first, so that a decay is a straight line drawn before
  // anything is bent and its daughter is walked like any other track.
  if (secondaries.decayYield > 0.f) {
    produceDecays(rng, secondaries, m_cfg.bFieldZ, tracks, particles,
                  numPrimaries);
  }

  const float maxGamma = 2.f * pi * propagation.maxTurns;
  const std::size_t maxTracks =
      totalForEvent(propagation.maxTracksPerPrimary, numPrimaries);
  const bool hasMaterial = cfg.material.scale > 0.f;
  const bool scattering = hasMaterial && cfg.material.multipleScattering;
  const bool energyLoss = hasMaterial && cfg.material.energyLoss;
  const bool meanLoss = cfg.material.energyLossModel == EnergyLossModel::Mean;
  // Either effect accumulates along the path, so the crossings have to be
  // walked in the order the track meets them.
  const bool needsSortedCrossings = scattering || energyLoss;

  for (std::size_t t = 0; t < tracks.size(); ++t) {
    // the vector grows while secondaries are queued, so work on a copy
    const Helix track = tracks[t];
    const TrackDirection direction = trackDirection(track.cotTheta);
    // what it was recorded as when it was made; only its hits are still to come
    const auto particle = static_cast<std::uint32_t>(base + t);
    // what a daughter of it has to share; it falls as the track loses energy
    float trackP = particles[particle].pt * direction.invSinTheta;

    crossings.clear();
    // A stiff track is stopped where it leaves the tracker, a curling one only
    // by `maxTurns`.
    propagateHelix(layout, track, crossings,
                   particles[particle].productionSurface,
                   helixEscapeGamma(track, layout.escapeRadius,
                                    layout.escapeHalfZ, maxGamma));
    if (needsSortedCrossings) {
      // they come out in surface order, not in the order the track meets them
      std::ranges::sort(crossings, {}, &SurfaceCrossing::gamma);
    }

    // Running sums of the scattering angles and their moments along the path,
    // which is all the displacement `sum(theta_i * (s - s_i))` needs.
    float sumAngleU = 0.f;
    float sumAngleV = 0.f;
    float sumMomentU = 0.f;
    float sumMomentV = 0.f;
    // Radiation lengths crossed so far, which is what the logarithmic term of
    // Highland's formula is a function of.
    float sumXOverX0 = 0.f;
    // The same for the fractional momentum the track has lost, whose
    // displacement `sum(f_i * (gamma - gamma_i)^2)` needs one moment more.
    float sumLoss = 0.f;
    float sumLossGamma = 0.f;
    float sumLossGamma2 = 0.f;

    for (const SurfaceCrossing& crossing : crossings) {
      const float gamma = crossing.gamma;
      const float pathAlongTrack = track.radius * gamma * direction.invSinTheta;

      // How far the material crossed so far has moved the track across itself.
      float displaceU = 0.f;
      float displaceV = 0.f;
      if (scattering) {
        displaceU = sumAngleU * pathAlongTrack - sumMomentU;
        displaceV = sumAngleV * pathAlongTrack - sumMomentV;
      }
      if (energyLoss) {
        // a slower track turns faster, and the hit slides by the double
        // integral of that over the transverse arc `radius * gamma`
        displaceU -= 0.5f * track.charge * track.radius *
                     (sumLoss * gamma * gamma - 2.f * sumLossGamma * gamma +
                      sumLossGamma2);
      }

      // what this crossing is worth, in the length each effect counts in
      const CrossingMaterial material =
          crossingMaterial(cfg, layout.surfaces[crossing.surface], crossing);

      if (crossing.sensitive()) {
        particles[particle].numHits +=
            recordHits(rng, layout, cfg.measurement, track, direction, crossing,
                       displaceU, displaceV, particle, hits);
      }

      if (scattering && material.xOverX0 > 0.f) {
        // Highland describes a *small* deflection, and past this the
        // displacement model says nothing anyway.
        constexpr float maxScatteringAngle = 0.2f;
        constexpr float highland = 13.6_MeV;
        // the log correction is a function of the whole path so far, not of
        // the slab in front; per crossing it under-scatters by a fifth
        sumXOverX0 += material.xOverX0;
        const float theta0 =
            std::min(highland / trackP * std::sqrt(material.xOverX0) *
                         std::max(0.f, 1.f + 0.038f * std::log(sumXOverX0)),
                     maxScatteringAngle);
        auto [angleU, angleV] = sampleNormalPair(rng);
        angleU *= theta0;
        angleV *= theta0;
        sumAngleU += angleU;
        sumAngleV += angleV;
        sumMomentU += angleU * pathAlongTrack;
        sumMomentV += angleV * pathAlongTrack;
      }

      // Whether the track paid for this crossing with everything it had.
      bool stopped = false;
      if (energyLoss && material.xOverX0 > 0.f) {
        Acts::MaterialSlab slab = *crossing.material;
        slab.scaleThickness(material.scale);
        const float qOverP = track.charge / trackP;
        const float loss = meanLoss ? Acts::computeEnergyLossMean(
                                          slab, m_absPdg, m_mass, qOverP, 1.f)
                                    : Acts::computeEnergyLossMode(
                                          slab, m_absPdg, m_mass, qOverP, 1.f);
        // a track that cannot pay for the surface in front of it stops inside
        // it, which is how a soft particle ranges out
        const float energy = std::hypot(trackP, m_mass);
        stopped = loss >= energy - m_mass;
        // dp = dE / beta, and the direction is untouched, so the transverse
        // momentum and the whole of it lose the same fraction
        const float fraction =
            stopped ? 1.f : loss * energy / (trackP * trackP);
        // the slide is first order in what was lost, so past
        // `maxEnergyLossFraction` it stops growing while the momentum falls on
        const float counted = std::min(
            fraction,
            std::max(0.f, cfg.material.maxEnergyLossFraction - sumLoss));
        sumLoss += counted;
        sumLossGamma += counted * gamma;
        sumLossGamma2 += counted * gamma * gamma;
        trackP *= 1.f - fraction;
      }

      // deeper generations are off by default; see
      // `SecondaryConfig::maxGenerations`
      if (particles[particle].generation < secondaries.maxGenerations &&
          tracks.size() < maxTracks) {
        produceSecondaries(rng, secondaries, m_cfg.bFieldZ, track,
                           particles[particle].generation, crossing, material,
                           trackP, tracks, particles, stubs);
      }

      if (stopped) {
        // it made whatever this crossing was worth and goes no further
        break;
      }
    }
  }

  // The stubs, once every propagated track has its particle index. Each leaves
  // its space points on the layer that made it and nowhere else.
  for (const PendingStub& stub : stubs) {
    const auto particle = static_cast<std::uint32_t>(particles.size());
    // it goes nowhere, so its direction is where it was made and its impact
    // parameter the radius it can no longer come back from
    const float stubR = std::max(stub.r, 1e-3f);
    GeneratedParticle info;
    info.pt = stub.pt;
    info.eta = std::asinh(stub.z / stubR);
    info.phi = std::remainder(stub.phi, 2.f * pi);
    info.d0 = stub.r;
    info.z0 = stub.z;
    info.charge = stub.charge;
    info.productionRadius = stub.r;
    info.productionZ = stub.z;
    info.productionSurface = stub.surface;
    info.generation = stub.generation;
    particles.push_back(info);

    const bool cylinder =
        layout.layers[stub.layer].shape == SurfaceShape::Cylinder;
    const std::uint32_t numClusters =
        1 + samplePoisson(rng, std::max(0.f, secondaries.stubClusters - 1.f));
    for (std::uint32_t h = 0; h < numClusters; ++h) {
      // the curl swamps the sensor resolution by two orders of magnitude, so
      // the reach is the whole of the offset
      auto [alongRZ, alongRPhi] = sampleNormalPair(rng);
      alongRZ *= secondaries.stubReach;
      alongRPhi *= secondaries.stubReach;
      const float hitR = stub.r + (cylinder ? 0.f : alongRZ);
      const float hitZ = stub.z + (cylinder ? alongRZ : 0.f);
      const std::optional<std::uint32_t> layer =
          findLayer(layout, layout.surfaces[stub.surface], hitR, hitZ);
      if (!layer.has_value()) {
        continue;
      }
      hits.push_back(SmearedHit{hitR, stub.phi + alongRPhi / stubR, hitZ,
                                *layer, particle});
      ++particles[particle].numHits;
    }
  }
}

}  // namespace ActsFatras::Synthetic::detail
