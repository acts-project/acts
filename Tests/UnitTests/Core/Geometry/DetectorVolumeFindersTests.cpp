// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Detector.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Geometry/detail/DetectorVolumeFinders.hpp"
#include "Acts/Geometry/detail/PortalGenerators.hpp"
#include "Acts/Geometry/detail/SurfaceCandidatesUpdators.hpp"
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

auto portalGenerator = Acts::Experimental::detail::defaultPortalGenerator();

auto cyl0 = Acts::Experimental::DetectorVolumeFactory::construct(
    portalGenerator, tContext, "Cyl0", nominal, std::move(cyl0Bounds),
    Acts::Experimental::detail::allPortals());

auto cyl1 = Acts::Experimental::DetectorVolumeFactory::construct(
    portalGenerator, tContext, "Cyl1", nominal, std::move(cyl1Bounds),
    Acts::Experimental::detail::allPortals());

auto cyl2 = Acts::Experimental::DetectorVolumeFactory::construct(
    portalGenerator, tContext, "Cyl2", nominal, std::move(cyl2Bounds),
    Acts::Experimental::detail::allPortals());

std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes012 = {
    cyl0, cyl1, cyl2};

auto det012 = Acts::Experimental::Detector::makeShared(
    "Det012", volumes012, Acts::Experimental::detail::tryAllVolumes());

BOOST_AUTO_TEST_SUITE(Experimental)

// Test finding detectors beu trial and error
BOOST_AUTO_TEST_CASE(TrialAndErrorDetectorVolumeFinder) {
  nState.currentDetector = det012.get();
  Acts::Experimental::detail::TrialAndErrorImpl tae;
  // Cylinder 0
  nState.position = Acts::Vector3(5., 0., 0.);
  tae.update(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl0.get());
  // Cylinder 1
  nState.position = Acts::Vector3(50., 0., 0.);
  tae.update(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl1.get());
  // Cylinder 2
  nState.position = Acts::Vector3(150., 0., 0.);
  tae.update(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl2.get());

  nState.currentDetector = nullptr;
  BOOST_CHECK_THROW(tae.update(tContext, nState), std::runtime_error);
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

  Acts::Experimental::detail::IndexedDetectorVolumeImpl<decltype(g)> idv(
      std::move(g), {Acts::binR});

  // Cylinder 0
  nState.position = Acts::Vector3(5., 0., 0.);
  idv.update(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl0.get());
  // Cylinder 1
  nState.position = Acts::Vector3(50., 0., 0.);
  idv.update(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl1.get());
  // Cylinder 2
  nState.position = Acts::Vector3(150., 0., 0.);
  idv.update(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl2.get());

  nState.currentDetector = nullptr;
  BOOST_CHECK_THROW(idv.update(tContext, nState), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
