// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfUtils.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <cstddef>
#include <cstdint>
#include <span>

namespace Acts {

double detail::Gsf::calculateDeterminant(
    const double *fullCalibratedCovariance,
    TrackStateTraits<kMeasurementSizeMax, true>::Covariance predictedCovariance,
    BoundSubspaceIndices projector, unsigned int calibratedSize) {
  return visit_measurement(calibratedSize, [&](auto N) {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;
    std::span<const std::uint8_t, kMeasurementSize> validSubspaceIndices(
        projector.begin(), projector.begin() + kMeasurementSize);
    FixedBoundSubspaceHelper<kMeasurementSize> subspaceHelper(
        validSubspaceIndices);

    typename Acts::TrackStateTraits<
        kMeasurementSize, true>::CalibratedCovariance calibratedCovariance{
        fullCalibratedCovariance};

    const auto H = subspaceHelper.projector();

    return (H * predictedCovariance * H.transpose() + calibratedCovariance)
        .determinant();
  });
}

void detail::Gsf::removeLowWeightComponents(std::vector<GsfComponent> &cmps,
                                            double weightCutoff) {
  auto proj = [](auto &cmp) -> double & { return cmp.weight; };

  normalizeWeights(cmps, proj);

  auto newEnd = std::remove_if(cmps.begin(), cmps.end(), [&](auto &cmp) {
    return proj(cmp) < weightCutoff;
  });

  // In case we would remove all components, keep only the largest
  if (std::distance(cmps.begin(), newEnd) == 0) {
    cmps = {*std::max_element(cmps.begin(), cmps.end(), [&](auto &a, auto &b) {
      return proj(a) < proj(b);
    })};
    cmps.front().weight = 1.0;
  } else {
    cmps.erase(newEnd, cmps.end());
    normalizeWeights(cmps, proj);
  }
}

double detail::Gsf::applyBetheHeitler(
    const GeometryContext &geoContext, const Surface &surface,
    Direction direction, const BoundTrackParameters &initialParameters,
    double initialWeight, const BetheHeitlerApprox &betheHeitlerApprox,
    std::vector<BetheHeitlerApprox::Component> &betheHeitlerCache,
    double weightCutoff, double transverseMomentumCut,
    std::vector<GsfComponent> &componentCache,
    Updatable<std::size_t> &nInvalidBetheHeitler,
    Updatable<double> &maxPathXOverX0, const Logger &logger) {
  const double initialMomentum = initialParameters.absoluteMomentum();
  const ParticleHypothesis &particleHypothesis =
      initialParameters.particleHypothesis();

  // Evaluate material slab
  MaterialSlab slab = surface.surfaceMaterial()->materialSlab(
      initialParameters.position(geoContext), direction,
      MaterialUpdateMode::FullUpdate);

  const double pathCorrection =
      surface.pathCorrection(geoContext, initialParameters.position(geoContext),
                             initialParameters.direction());
  slab.scaleThickness(pathCorrection);

  const double pathXOverX0 = slab.thicknessInX0();
  maxPathXOverX0.tmp() = std::max(maxPathXOverX0.tmp(), pathXOverX0);

  // Emit a warning if the approximation is not valid for this x/x0
  if (!betheHeitlerApprox.validXOverX0(pathXOverX0)) {
    ++nInvalidBetheHeitler.tmp();
    ACTS_DEBUG("Bethe-Heitler approximation encountered invalid value for x/x0="
               << pathXOverX0 << " at surface " << surface.geometryId());
  }

  // Get the mixture
  betheHeitlerCache.resize(betheHeitlerApprox.maxComponents());
  const auto mixture =
      betheHeitlerApprox.mixture(pathXOverX0, betheHeitlerCache);

  // Create all possible new components
  for (const GaussianComponent &gaussian : mixture) {
    // Here we combine the new child weight with the parent weight.
    // However, this must be later re-adjusted
    const double newWeight = gaussian.weight * initialWeight;

    if (newWeight < weightCutoff) {
      ACTS_VERBOSE("Skip component with weight " << newWeight);
      continue;
    }

    if (gaussian.mean < 1.e-8) {
      ACTS_WARNING("Skip component with gaussian " << gaussian.mean << " +- "
                                                   << gaussian.var);
      continue;
    }

    // compute delta p from mixture and update parameters
    BoundVector newPars = initialParameters.parameters();

    const auto delta_p = [&]() {
      if (direction == Direction::Forward()) {
        return initialMomentum * (gaussian.mean - 1.);
      } else {
        return initialMomentum * (1. / gaussian.mean - 1.);
      }
    }();

    assert(initialMomentum + delta_p > 0. && "new momentum must be > 0");

    // Apply pT cut here to avoid const of expansion and merging later
    const auto pT =
        (initialMomentum + delta_p) * std::sin(newPars[eBoundTheta]);
    if (pT < transverseMomentumCut) {
      ACTS_VERBOSE("Skip new component with pT=" << pT << " GeV");
      continue;
    }

    newPars[eBoundQOverP] = particleHypothesis.qOverP(
        initialMomentum + delta_p, initialParameters.charge());

    // compute inverse variance of p from mixture and update covariance
    BoundMatrix newCov = initialParameters.covariance().value();

    const auto varInvP = [&]() {
      if (direction == Direction::Forward()) {
        const double f = 1. / (initialMomentum * gaussian.mean);
        return f * f * gaussian.var;
      } else {
        return gaussian.var / (initialMomentum * initialMomentum);
      }
    }();

    newCov(eBoundQOverP, eBoundQOverP) += varInvP;
    assert(std::isfinite(newCov(eBoundQOverP, eBoundQOverP)) &&
           "new cov not finite");

    // Set the remaining things and push to vector
    componentCache.push_back({newWeight, newPars, newCov});
  }

  return pathXOverX0;
}

}  // namespace Acts
