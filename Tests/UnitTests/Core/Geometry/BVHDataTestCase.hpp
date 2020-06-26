// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace bdata = boost::unit_test::data;
using Box = Acts::Volume::BoundingBox;
using Ray = Acts::Ray<double, 3>;

GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

std::tuple<std::vector<const Volume*>, std::shared_ptr<TrackingGeometry>>
gridBoxFactory(size_t n = NBOXES, double hl = 1000, size_t octd = 5) {
  Box::Size size(Acts::Vector3D(2, 2, 2));

  std::shared_ptr<CuboidVolumeBounds> vbds =
      std::make_shared<CuboidVolumeBounds>(10, 10, 10);

  double min = -hl;
  double max = hl;

  double step = (max - min) / double(n);
  std::vector<std::unique_ptr<const Volume>> volumes;
  std::vector<std::unique_ptr<Box>> boxStore;
  boxStore.reserve((n + 1) * (n + 1) * (n + 1));

  std::cout << "generating: " << (n + 1) * (n + 1) * (n + 1)
            << " bounding boxes" << std::endl;

  std::vector<Box*> boxes;
  boxes.reserve(boxStore.size());

  for (size_t i = 0; i <= n; i++) {
    for (size_t j = 0; j <= n; j++) {
      for (size_t k = 0; k <= n; k++) {
        Vector3D pos(min + i * step, min + j * step, min + k * step);

        auto trf = std::make_shared<Transform3D>(Translation3D(pos));
        auto vol = std::make_unique<AbstractVolume>(trf, vbds);

        volumes.push_back(std::move(vol));
        boxStore.push_back(
            std::make_unique<Box>(volumes.back()->boundingBox()));
        boxes.push_back(boxStore.back().get());
      }
    }
  }

  Box* top = make_octree(boxStore, boxes, octd);

  // create trackingvolume
  // will own the volumes, so make non-owning copy first
  std::vector<const Volume*> volumeCopy;
  volumeCopy.reserve(volumes.size());
  for (auto& vol : volumes) {
    volumeCopy.push_back(vol.get());
  }

  // box like overall shape
  auto tvTrf = std::make_shared<Transform3D>(Transform3D::Identity());
  auto tvBounds =
      std::make_shared<CuboidVolumeBounds>(hl * 1.1, hl * 1.1, hl * 1.1);

  auto tv =
      TrackingVolume::create(tvTrf, tvBounds, std::move(boxStore),
                             std::move(volumes), top, nullptr, "TheVolume");

  auto tg = std::make_shared<TrackingGeometry>(tv);

  return {std::move(volumeCopy), tg};
}

auto [volumes, tg] = gridBoxFactory();

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
  Acts::Vector3D dir;
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
      auto sri = srf.intersect(tgContext, ray.origin(), ray.dir(), true);
      if (sri and sri.intersection.pathLength >= s_onSurfaceTolerance) {
        // does intersect
        hits.push_back(std::move(sri));
      }
    }
  }

  // sort by path length
  std::sort(hits.begin(), hits.end());
  std::vector<const Surface*> expHits;
  expHits.reserve(hits.size());
  for (const auto& hit : hits) {
    expHits.push_back(hit.object);
  }

  // now do the same through a propagator
  using SteppingLogger = Acts::detail::SteppingLogger;
  using Stepper = StraightLineStepper;
  using PropagatorType = Propagator<Stepper, Navigator>;

  Stepper stepper{};
  Navigator navigator(tg);
  PropagatorType propagator(std::move(stepper), navigator);

  using DebugOutput = Acts::DebugOutputActor;
  using ActionList = Acts::ActionList<SteppingLogger, DebugOutput>;
  using AbortConditions = Acts::AbortList<>;

  Acts::PropagatorOptions<ActionList, AbortConditions> options(tgContext,
                                                               mfContext);

  options.debug = false;
  options.pathLimit = 20_m;

  // this should be irrelevant.
  double mom = 50_GeV;

  Acts::CurvilinearParameters startPar(std::nullopt, ray.origin(),
                                       ray.dir() * mom, +1, 0.);

  const auto result = propagator.propagate(startPar, options).value();

  const auto debugString =
      result.template get<DebugOutput::result_type>().debugString;

  if (options.debug) {
    std::cout << debugString << std::endl;
  }

  // collect surfaces
  std::vector<const Surface*> actHits;
  auto steppingResults =
      result.template get<SteppingLogger::result_type>().steps;
  for (const auto& step : steppingResults) {
    if (!step.surface) {
      continue;
    }

    auto sensitiveID = step.surface->geoID().sensitive();
    if (sensitiveID != 0) {
      actHits.push_back(step.surface.get());
    }
  }

  BOOST_CHECK_EQUAL(expHits.size(), actHits.size());
  for (size_t i = 0; i < expHits.size(); i++) {
    const Surface* exp = expHits[i];
    const Surface* act = actHits[i];

    BOOST_CHECK_EQUAL(exp, act);
    BOOST_CHECK_EQUAL(exp->geoID(), act->geoID());
  }
}
