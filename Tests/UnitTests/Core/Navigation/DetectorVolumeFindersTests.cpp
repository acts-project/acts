// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/SurfaceCandidatesDelegates.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <array>
#include <memory>
#include <vector>

// A test context
Acts::GeometryContext tContext;

Acts::Experimental::NavigationState nState;

Acts::ActsScalar r0 = 0.;
Acts::ActsScalar r1 = 10.;
Acts::ActsScalar r2 = 100.;
Acts::ActsScalar r3 = 200.;
Acts::ActsScalar zHalfL = 200.;

Acts::Transform3 nominal = Acts::Transform3::Identity();

// Create a bunch of volumes
auto cyl0Bounds = std::make_unique<Acts::CylinderVolumeBounds>(r0, r1, zHalfL);

auto cyl1Bounds = std::make_unique<Acts::CylinderVolumeBounds>(r1, r2, zHalfL);

auto cyl2Bounds = std::make_unique<Acts::CylinderVolumeBounds>(r2, r3, zHalfL);

auto cyl0 = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "Cyl0", nominal, std::move(cyl0Bounds),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());

auto cyl1 = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "Cyl1", nominal, std::move(cyl1Bounds),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());

auto cyl2 = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "Cyl2", nominal, std::move(cyl2Bounds),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());

std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes012 = {
    cyl0, cyl1, cyl2};
auto det012 = Acts::Experimental::Detector::makeShared(
    "Det012", volumes012,
    Acts::Experimental::makeDetectorVolumeFinder<
        const Acts::Experimental::RootVolumeFinder>());

auto volumeA = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "A", Acts::Transform3::Identity(),
    std::make_unique<Acts::CylinderVolumeBounds>(0, 10, 200),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());
auto volumeB = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "B", Acts::Transform3::Identity(),
    std::make_unique<Acts::CylinderVolumeBounds>(10, 100, 200),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());
auto volumeC = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "C", Acts::Transform3::Identity(),
    std::make_unique<Acts::CylinderVolumeBounds>(10200, 10, 200),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());
auto volumeD = Acts::Experimental::DetectorVolumeFactory::construct(
    Acts::Experimental::makePortalGenerator<
        const Acts::Experimental::DefaultPortalGenerator>(),
    tContext, "D", Acts::Transform3::Identity(),
    std::make_unique<Acts::CylinderVolumeBounds>(0, 10, 200),
    Acts::Experimental::makeSurfaceCandidatesDelegate<
        const Acts::Experimental::AllPortals>());

BOOST_AUTO_TEST_SUITE(Experimental)

// Test finding detectors by trial and error
BOOST_AUTO_TEST_CASE(RootVolumeFinder) {
  nState.currentDetector = det012.get();
  Acts::Experimental::RootVolumeFinder rvf;
  // Cylinder 0
  nState.position = Acts::Vector3(5., 0., 0.);
  rvf.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl0.get());
  // Cylinder 1
  nState.position = Acts::Vector3(50., 0., 0.);
  rvf.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl1.get());
  // Cylinder 2
  nState.position = Acts::Vector3(150., 0., 0.);
  rvf.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl2.get());

  nState.currentDetector = nullptr;
  BOOST_CHECK_THROW(rvf.find(tContext, nState), std::runtime_error);
}

// Test finding detectors beu trial and error
BOOST_AUTO_TEST_CASE(IndexedDetectorVolumeFinder) {
  nState.currentDetector = det012.get();

  using SingleIndex = std::size_t;

  using Axis = Acts::detail::Axis<Acts::detail::AxisType::Variable,
                                  Acts::detail::AxisBoundaryType::Bound>;
  using Grid = Acts::detail::Grid<SingleIndex, Axis>;

  std::vector<Acts::ActsScalar> b = {r0, r1, r2, r3};
  Axis a(b);
  Grid g(std::make_tuple(a));

  g.atPosition(std::array<Acts::ActsScalar, 1u>{5.}) = 0u;
  g.atPosition(std::array<Acts::ActsScalar, 1u>{50.}) = 1u;
  g.atPosition(std::array<Acts::ActsScalar, 1u>{150.}) = 2u;

  auto idv = Acts::Experimental::makeIndexedDetectorVolumeFinder<decltype(g)>(
      std::move(g), {Acts::binR});

  // Cylinder 0
  nState.position = Acts::Vector3(5., 0., 0.);
  idv.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl0.get());
  // Cylinder 1
  nState.position = Acts::Vector3(50., 0., 0.);
  idv.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl1.get());
  // Cylinder 2
  nState.position = Acts::Vector3(150., 0., 0.);
  idv.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl2.get());

  nState.currentDetector = nullptr;
  BOOST_CHECK_THROW(idv.find(tContext, nState), std::runtime_error);
}

// These tests check the behavior of the volume finders, i.e. the
// helper delegates that set/reset the volume raw pointer in the
// NavigaitonState according to some given information.
BOOST_AUTO_TEST_CASE(Unconnected) {
  Acts::Experimental::DetectorVolumeFinder ucFinder;
  BOOST_CHECK(not ucFinder.connected());
}

// The end of world is reached
BOOST_AUTO_TEST_CASE(EndOfWorld) {
  nState.currentVolume = volumeA.get();
  BOOST_CHECK(nState.currentVolume == volumeA.get());

  Acts::Experimental::EndOfWorldVolume eow;
  eow.find(tContext, nState);

  BOOST_CHECK(nState.currentVolume == nullptr);
}

// A single link exists and this is set
BOOST_AUTO_TEST_CASE(SingleVolume) {
  nState.currentVolume = volumeA.get();
  BOOST_CHECK(nState.currentVolume == volumeA.get());

  Acts::Experimental::SingleDetectorVolume svu(volumeB.get());
  svu.find(tContext, nState);

  BOOST_CHECK(nState.currentVolume == volumeB.get());

  BOOST_CHECK_THROW(Acts::Experimental::SingleDetectorVolume(nullptr),
                    std::invalid_argument);
}

// A typlical volume array in 1 dimension (bound, not closed)
BOOST_AUTO_TEST_CASE(VolumeArray) {
  std::vector<Acts::ActsScalar> zArray = {-200, -100, 100, 400, 1000};

  std::vector<const Acts::Experimental::DetectorVolume*> volumes = {
      volumeA.get(), volumeB.get(), volumeC.get(), volumeD.get()};
  Acts::Experimental::BoundVolumesGrid1 bvg(zArray, Acts::binZ, volumes);
  // Reset the navigation state
  nState.currentVolume = nullptr;

  // Check the volume retrieval
  nState.position = Acts::Vector3(0., 0., -150.);
  bvg.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeA.get());

  nState.position = Acts::Vector3(0., 0., 600.);
  bvg.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeD.get());

  // Check a shifted one
  Acts::Transform3 shift300 = Acts::Transform3::Identity();
  shift300.pretranslate(Acts::Vector3(0, 0, 300));

  Acts::Experimental::BoundVolumesGrid1 bvgs(zArray, Acts::binZ, volumes,
                                             shift300.inverse());

  // 150 (-300) -> transforms to -150, hence it yields A
  nState.position = Acts::Vector3(0., 0., 150.);
  bvgs.find(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeA.get());
}

BOOST_AUTO_TEST_SUITE_END()
