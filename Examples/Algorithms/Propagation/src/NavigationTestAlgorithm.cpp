// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Propagation/NavigationTestAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <random>
#include <stdexcept>

using namespace Acts::UnitLiterals;

namespace ActsExamples {

ProcessCode NavigationTestAlgorithm::execute(
    const AlgorithmContext& context) const {
  // Create a random number generator
  ActsExamples::RandomEngine rng =
      m_cfg.randomNumberSvc->spawnGenerator(context);

  std::uniform_real_distribution<double> rDist(m_rMin, m_rMax);
  std::uniform_real_distribution<double> zDist(-m_halfLengthZ, m_halfLengthZ);
  std::uniform_real_distribution<double> phiDist(-M_PI, M_PI);
  std::uniform_real_distribution<double> etaDist(-5, 5);

  Acts::Navigator navigator{{m_cfg.trackingGeometry},
                            logger().cloneWithSuffix("Nav")};
  Acts::EigenStepper<> stepper{m_cfg.magneticField};
  Acts::Propagator propagator{std::move(stepper), std::move(navigator),
                              logger().cloneWithSuffix("Prop")};

  // loop over number of particles
  for (size_t it = 0; it < m_cfg.ntests; ++it) {
    double z = zDist(rng);
    double r = rDist(rng);
    double phi = phiDist(rng);
    ACTS_VERBOSE("Start location: " << it << " at r = " << r << ", z = " << z
                                    << ", phi = " << phi);
    std::shared_ptr<const Acts::PerigeeSurface> startSurface =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3(r * std::cos(phi), r * std::sin(phi), z));
    auto t = std::tie(*startSurface, context.geoContext);
    ACTS_VERBOSE("Start surface: " << t);

    z = zDist(rng);
    r = rDist(rng);
    phi = phiDist(rng);
    ACTS_VERBOSE("Target location " << it << " at r = " << r << ", z = " << z
                                    << ", phi = " << phi);
    std::shared_ptr<const Acts::PerigeeSurface> targetSurface =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3(r * std::cos(phi), r * std::sin(phi), z));
    ACTS_VERBOSE(
        "Target surface: " << std::tie(*targetSurface, context.geoContext));

    double dirPhi = phiDist(rng);
    double dirEta = etaDist(rng);
    double dirTheta = 2 * std::atan(std::exp(-dirEta));

    Acts::Vector3 dir{std::cos(dirPhi) * std::sin(dirTheta),
                      std::sin(dirPhi) * std::sin(dirTheta),
                      std::cos(dirTheta)};

    ACTS_VERBOSE("Direction: " << it << " phi = " << dirPhi << ", eta = "
                               << dirEta << ", theta = " << dirTheta << " -> "
                               << dir.transpose());

    Acts::Vector4 pos4;
    pos4.template head<3>() = startSurface->center(context.geoContext);
    pos4[3] = 0;

    auto intersections = targetSurface->intersect(
        context.geoContext, pos4.template head<3>(), dir, true);

    bool reachable = intersections.closest().status() ==
                         Acts::Intersection3D::Status::reachable &&
                     intersections.size() > 0;

    ACTS_VERBOSE("Surface is reachable a priori? "
                 << (reachable ? "yes" : "no"));
    // if (intersections.size() != 1) {
    // ACTS_ERROR("Target surface does not intersect with start parameters");
    // return ProcessCode::ABORT;
    // }
    // if (intersections.closest().status() !=
    // Acts::Intersection3D::Status::reachable) {
    // ACTS_ERROR("Target surface is not reachable");
    // return ProcessCode::ABORT;
    // }

    auto boundParams =
        Acts::BoundTrackParameters::create(startSurface, context.geoContext,
                                           pos4, dir, 1 / 2_GeV, std::nullopt,
                                           Acts::ParticleHypothesis::pion())
            .value();

    Acts::PropagatorOptions<> pOptions{context.geoContext,
                                       context.magFieldContext};

    auto pRes =
        propagator.template propagate(boundParams, *targetSurface, pOptions);

    if (!pRes.ok() && reachable) {
      ACTS_ERROR("Propagation failed: " << pRes.error());
      return ProcessCode::ABORT;
    }
  }

  return ProcessCode::SUCCESS;
}

NavigationTestAlgorithm::NavigationTestAlgorithm(
    const NavigationTestAlgorithm::Config& config, Acts::Logging::Level level)
    : IAlgorithm("NavigationTestAlgorithm", level), m_cfg(config) {
  if (!m_cfg.randomNumberSvc) {
    throw std::invalid_argument("No random number generator given");
  }

  const auto* topVolume = m_cfg.trackingGeometry->highestTrackingVolume();
  const auto& volumeBounds = topVolume->volumeBounds();

  const auto* cylVolBounds =
      dynamic_cast<const Acts::CylinderVolumeBounds*>(&volumeBounds);

  if (cylVolBounds == nullptr) {
    throw std::runtime_error{"Volume bounds are not cylindrical"};
  }

  const auto values = cylVolBounds->values();
  m_rMin = values[Acts::CylinderVolumeBounds::eMinR];
  m_rMax = values[Acts::CylinderVolumeBounds::eMaxR];
  m_halfLengthZ = values[Acts::CylinderVolumeBounds::eHalfLengthZ];

  ACTS_INFO("Cylinder bounds extracted: rMin = "
            << m_rMin << ", rMax = " << m_rMax
            << ", halfLengthZ = " << m_halfLengthZ);
}

}  // namespace ActsExamples
