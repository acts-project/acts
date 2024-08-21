// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/NavigationStreamHelper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

namespace {

std::vector<std::shared_ptr<Acts::Surface>> createPlaneSurfaces() {
  auto rectangle = std::make_shared<Acts::RectangleBounds>(10., 10.);
  // Surface A:
  // This surface should not be reachable from (0.,0.,0.) position along z
  Acts::Transform3 aTransform = Acts::Transform3::Identity();
  aTransform.pretranslate(Acts::Vector3(0., 0., -20.));
  auto surfaceA =
      Acts::Surface::makeShared<Acts::PlaneSurface>(aTransform, rectangle);
  // Surface B:
  // This surface should not be reachable from (0.,0.,0.) position along z with
  // boundary check
  Acts::Transform3 bTransform = Acts::Transform3::Identity();
  bTransform.pretranslate(Acts::Vector3(50., 50., 100.));
  auto surfaceB =
      Acts::Surface::makeShared<Acts::PlaneSurface>(bTransform, rectangle);
  // Surface C:
  Acts::Transform3 cTransform = Acts::Transform3::Identity();
  cTransform.pretranslate(Acts::Vector3(0., 0., 200.));
  auto surfaceC =
      Acts::Surface::makeShared<Acts::PlaneSurface>(cTransform, rectangle);
  // Surface D:
  Acts::Transform3 dTransform = Acts::Transform3::Identity();
  dTransform.pretranslate(Acts::Vector3(0., 0., 400.));
  auto surfaceD =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);

  // Let's fill them shuffled
  return {surfaceC, surfaceA, surfaceD, surfaceB};
}

}  // namespace

using namespace Acts;

auto gContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Navigation)

BOOST_AUTO_TEST_CASE(NavigationStream_InitializePlanes) {
  // Create the punch of surfaces
  auto surfaces = createPlaneSurfaces();

  NavigationStream nStreamTemplate;
  for (const auto& surface : surfaces) {
    NavigationStreamHelper::fillSurfaces(nStreamTemplate, {surface.get()},
                                         Acts::BoundaryTolerance::None());
  }
  BOOST_CHECK_EQUAL(nStreamTemplate.candidates.size(), 4);

  // (1) Run an initial update
  // - from a position where all are reachable and valid
  // - with infinite boundary tolerance
  NavigationStream nStream = nStreamTemplate;
  BOOST_CHECK(NavigationStreamHelper::initializeStream(
      nStream, gContext, Acts::Vector3(0., 0., -30.), Acts::Vector3(0., 0., 1.),
      BoundaryTolerance::Infinite()));

  BOOST_CHECK_EQUAL(nStream.activeCandidates(), 4u);
  BOOST_CHECK_EQUAL(&nStream.currentCandidate().surface(), surfaces[1u].get());

  // (2) Run an initial update
  // - from a position where all but one are reachable
  // - with infinite boundary tolerance
  nStream = nStreamTemplate;
  BOOST_CHECK(NavigationStreamHelper::initializeStream(
      nStream, gContext, Acts::Vector3(0., 0., 0.), Acts::Vector3(0., 0., 1.),
      BoundaryTolerance::Infinite()));
  BOOST_CHECK_EQUAL(nStream.activeCandidates(), 3u);
  BOOST_CHECK_EQUAL(&nStream.currentCandidate().surface(), surfaces[3u].get());

  // (3) Run an initial update
  // - from a position where all would be reachable, but
  // - with no boundary tolerance
  nStream = nStreamTemplate;
  BOOST_CHECK(NavigationStreamHelper::initializeStream(
      nStream, gContext, Acts::Vector3(0., 0., -100.),
      Acts::Vector3(0., 0., 1.), BoundaryTolerance::None()));
  BOOST_CHECK_EQUAL(nStream.activeCandidates(), 3u);

  // (4) Run an initial update
  // - none of the surfaces should be reachable
  nStream = nStreamTemplate;
  BOOST_CHECK(!NavigationStreamHelper::initializeStream(
      nStream, gContext, Acts::Vector3(0., 0., 0.), Acts::Vector3(1., 0., 0.),
      BoundaryTolerance::Infinite()));
  BOOST_CHECK_EQUAL(nStream.activeCandidates(), 0u);
  BOOST_CHECK_THROW(nStream.currentCandidate(), std::out_of_range);


}

BOOST_AUTO_TEST_SUITE_END()