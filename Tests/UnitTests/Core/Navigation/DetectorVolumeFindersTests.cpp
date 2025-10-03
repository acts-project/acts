// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

using namespace Acts;

// A test context
GeometryContext tContext;

Experimental::NavigationState nState;

double r0 = 0.;
double r1 = 10.;
double r2 = 100.;
double r3 = 200.;
double zHalfL = 200.;

Transform3 nominal = Transform3::Identity();

// Create a bunch of volumes
auto cyl0Bounds = std::make_unique<CylinderVolumeBounds>(r0, r1, zHalfL);

auto cyl1Bounds = std::make_unique<CylinderVolumeBounds>(r1, r2, zHalfL);

auto cyl2Bounds = std::make_unique<CylinderVolumeBounds>(r2, r3, zHalfL);

auto portalGenerator = Experimental::defaultPortalGenerator();

auto cyl0 = Experimental::DetectorVolumeFactory::construct(
    portalGenerator, tContext, "Cyl0", nominal, std::move(cyl0Bounds),
    Experimental::tryAllPortals());

auto cyl1 = Experimental::DetectorVolumeFactory::construct(
    portalGenerator, tContext, "Cyl1", nominal, std::move(cyl1Bounds),
    Experimental::tryAllPortals());

auto cyl2 = Experimental::DetectorVolumeFactory::construct(
    portalGenerator, tContext, "Cyl2", nominal, std::move(cyl2Bounds),
    Experimental::tryAllPortals());

std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes012 = {
    cyl0, cyl1, cyl2};

Experimental::GeometryIdGenerator::Config generatorConfig;
Experimental::GeometryIdGenerator generator(
    generatorConfig,
    getDefaultLogger("SequentialIdGenerator", Logging::VERBOSE));

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(NavigationSuite)

// Test finding detectors by trial and error
BOOST_AUTO_TEST_CASE(RootVolumeFinder) {
  auto cache = generator.generateCache();
  for (auto& vol : volumes012) {
    generator.assignGeometryId(cache, *vol);
  }

  auto det012 = Experimental::Detector::makeShared(
      "Det012", volumes012, Experimental::tryRootVolumes());

  nState.currentDetector = det012.get();
  Experimental::RootVolumeFinder rvf;
  // Cylinder 0
  nState.position = Vector3(5., 0., 0.);
  rvf.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, cyl0.get());
  // Cylinder 1
  nState.position = Vector3(50., 0., 0.);
  rvf.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, cyl1.get());
  // Cylinder 2
  nState.position = Vector3(150., 0., 0.);
  rvf.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, cyl2.get());

  nState.currentDetector = nullptr;
  BOOST_CHECK_THROW(rvf.update(tContext, nState), std::runtime_error);
}

// Test finding detectors beu trial and error
BOOST_AUTO_TEST_CASE(IndexedDetectorVolumeFinder) {
  auto cache = generator.generateCache();
  for (auto& vol : volumes012) {
    generator.assignGeometryId(cache, *vol);
  }

  auto det012 = Experimental::Detector::makeShared(
      "Det012", volumes012, Experimental::tryRootVolumes());

  nState.currentDetector = det012.get();

  using SingleIndex = std::size_t;

  using Axis = Axis<AxisType::Variable, AxisBoundaryType::Bound>;
  using Grid = Grid<SingleIndex, Axis>;

  std::vector<double> b = {r0, r1, r2, r3};
  Axis a(b);
  Grid g(std::make_tuple(a));

  g.atPosition(std::array<double, 1u>{5.}) = 0u;
  g.atPosition(std::array<double, 1u>{50.}) = 1u;
  g.atPosition(std::array<double, 1u>{150.}) = 2u;

  Experimental::IndexedDetectorVolumesImpl<decltype(g)> idv(
      std::move(g), {AxisDirection::AxisR});

  // Cylinder 0
  nState.position = Vector3(5., 0., 0.);
  idv.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, cyl0.get());
  // Cylinder 1
  nState.position = Vector3(50., 0., 0.);
  idv.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, cyl1.get());
  // Cylinder 2
  nState.position = Vector3(150., 0., 0.);
  idv.update(tContext, nState);
  BOOST_CHECK_EQUAL(nState.currentVolume, cyl2.get());

  nState.currentDetector = nullptr;
  BOOST_CHECK_THROW(idv.update(tContext, nState), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
