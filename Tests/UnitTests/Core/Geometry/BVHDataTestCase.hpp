// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"

namespace bdata = boost::unit_test::data;
using Box = Acts::Volume::BoundingBox;
using Ray = Acts::Ray<double, 3>;

GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

Test::CubicBVHTrackingGeometry grid(NBOXES, 1000, 5);

auto volumes = grid.volumes;
auto tg = grid.trackingGeometry;

BOOST_DATA_TEST_CASE(
    bvhnavigation_test,
    bdata::random((bdata::seed = 7, bdata::engine = std::mt19937(),
                   bdata::distribution = std::uniform_real_distribution<>(-5,
                                                                          5))) ^
        bdata::random((bdata::seed = 2, bdata::engine = std::mt19937(),
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 3, bdata::engine = std::mt19937(),
                       bdata::distribution =
                           std::uniform_real_distribution<>(-100, 100))) ^
        bdata::random((bdata::seed = 4, bdata::engine = std::mt19937(),
                       bdata::distribution =
                           std::uniform_real_distribution<>(-100, 100))) ^
        bdata::random((bdata::seed = 5, bdata::engine = std::mt19937(),
                       bdata::distribution =
                           std::uniform_real_distribution<>(-100, 100))) ^
        bdata::xrange(NTESTS),
    eta, phi, x, y, z, index) {
  using namespace Acts::UnitLiterals;
  (void)index;

  // construct ray from parameters
  double theta = 2 * std::atan(std::exp(-eta));
  Acts::Vector3 dir;
  dir << std::cos(phi), std::sin(phi), 1. / std::tan(theta);
  dir.normalize();
  Ray ray({x, y, z}, dir);

  // naive collection: iterate over all the boxes
  std::vector<SurfaceIntersection> hits;
  for (const auto& vol : volumes) {
    const auto& absVol = dynamic_cast<const AbstractVolume&>(*vol);
    auto bndSurfaces = absVol.boundarySurfaces();
    // collect all surfaces that are hit
    for (const auto& bndSrf : bndSurfaces) {
      const auto& srf = bndSrf->surfaceRepresentation();
      auto srmi = srf.intersect(tgContext, ray.origin(), ray.dir(),
                                Acts::BoundaryCheck(true));
      for (const auto& sri : srmi.split()) {
        if (sri && sri.pathLength() >= s_onSurfaceTolerance) {
          // does intersect
          hits.push_back(sri);
        }
      }
    }
  }

  // sort by path length
  std::sort(hits.begin(), hits.end(), SurfaceIntersection::forwardOrder);
  std::vector<const Surface*> expHits;
  expHits.reserve(hits.size());
  for (const auto& hit : hits) {
    expHits.push_back(hit.object());
  }

  // now do the same through a propagator
  using SteppingLogger = Acts::detail::SteppingLogger;
  using Stepper = StraightLineStepper;
  using PropagatorType = Propagator<Stepper, Navigator>;

  Stepper stepper{};
  Navigator navigator({tg});
  PropagatorType propagator(stepper, navigator);

  using ActionList = Acts::ActionList<SteppingLogger>;
  using AbortConditions = Acts::AbortList<>;

  Acts::PropagatorOptions<ActionList, AbortConditions> options(tgContext,
                                                               mfContext);

  options.pathLimit = 20_m;

  Acts::Vector4 pos4 = Acts::Vector4::Zero();
  pos4.segment<3>(Acts::ePos0) = ray.origin();
  // momentum value should be irrelevant.
  Acts::CurvilinearTrackParameters startPar(
      pos4, ray.dir(), 1_e / 50_GeV, std::nullopt, ParticleHypothesis::pion());

  const auto result = propagator.propagate(startPar, options).value();

  // collect surfaces
  std::vector<const Surface*> actHits;
  auto steppingResults =
      result.template get<SteppingLogger::result_type>().steps;
  for (const auto& step : steppingResults) {
    if (!step.surface) {
      continue;
    }

    auto sensitiveID = step.surface->geometryId().sensitive();
    if (sensitiveID != 0) {
      actHits.push_back(step.surface.get());
    }
  }

  BOOST_CHECK_EQUAL(expHits.size(), actHits.size());
  for (size_t i = 0; i < expHits.size(); i++) {
    const Surface* exp = expHits[i];
    const Surface* act = actHits[i];

    BOOST_CHECK_EQUAL(exp, act);
    BOOST_CHECK_EQUAL(exp->geometryId(), act->geometryId());
  }
}
