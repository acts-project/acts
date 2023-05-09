// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/detail/SupportHelper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(CylindricalSupport) {
  // As a single cylinder
  auto singleSupport =
      Acts::Experimental::detail::SupportHelper::cylindricalSupport(
          Acts::Transform3::Identity(), {100., 400., M_PI, 0., 0., 0.}, 1u);
  BOOST_CHECK(singleSupport.size() == 1u);
  BOOST_CHECK(singleSupport[0u]->type() ==
              Acts::Surface::SurfaceType::Cylinder);

  // As a split cylinder
  auto splitSupport =
      Acts::Experimental::detail::SupportHelper::cylindricalSupport(
          Acts::Transform3::Identity(), {100., 400., M_PI, 0., 0., 0.}, 32u);
  BOOST_CHECK(splitSupport.size() == 32u);
  for (const auto& ss : splitSupport) {
    BOOST_CHECK(ss->type() == Acts::Surface::SurfaceType::Plane);
  }

  // As a split cylinder - sectoral
  auto splitSectoralSupport =
      Acts::Experimental::detail::SupportHelper::cylindricalSupport(
          Acts::Transform3::Identity(),
          {100., 400., 0.25 * M_PI, 0.75 * M_PI, 0., 0.}, 128u);
  BOOST_CHECK(splitSectoralSupport.size() == 128u);
  for (const auto& ss : splitSectoralSupport) {
    BOOST_CHECK(ss->type() == Acts::Surface::SurfaceType::Plane);
  }
}

BOOST_AUTO_TEST_CASE(DiscSupport) {
  // As a single disc
  auto singleSupport = Acts::Experimental::detail::SupportHelper::discSupport(
      Acts::Transform3::Identity(), {100., 400., M_PI, 0.}, 1u);
  BOOST_CHECK(singleSupport.size() == 1u);
  BOOST_CHECK(singleSupport[0u]->type() == Acts::Surface::SurfaceType::Disc);

  // As a split disc
  auto splitSupport = Acts::Experimental::detail::SupportHelper::discSupport(
      Acts::Transform3::Identity(), {100., 400., M_PI, 0.}, 32u);
  BOOST_CHECK(splitSupport.size() == 32u);
  for (const auto& ss : splitSupport) {
    BOOST_CHECK(ss->type() == Acts::Surface::SurfaceType::Plane);
  }

  // As a split disc - sectoral
  auto splitSectoralSupport =
      Acts::Experimental::detail::SupportHelper::discSupport(
          Acts::Transform3::Identity(), {100., 400., 0.5 * M_PI, 0.}, 16u);
  BOOST_CHECK(splitSectoralSupport.size() == 16u);
  for (const auto& ss : splitSectoralSupport) {
    BOOST_CHECK(ss->type() == Acts::Surface::SurfaceType::Plane);
  }
}

BOOST_AUTO_TEST_CASE(addCylinderSupport) {
  std::vector<std::shared_ptr<Acts::Surface>> lSurfaces;
  std::vector<size_t> assignToAll;

  // The Extent
  Acts::Extent lExtent;
  lExtent.set(Acts::binR, 100., 110.);
  lExtent.set(Acts::binZ, -400., -400.);

  // Get the main support parameters:
  // - doff .. offset (in r.z)
  // - demin, demax .. envelop min, max (in z,r)
  // - dphimin, dphimin .. envelop min, max (in phi)
  std::array<Acts::ActsScalar, 5u> sValues = {10, 5, 5, 0., 0.};
  // Add a single support cylinder
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Cylinder,
      sValues, std::nullopt, 1u);
  BOOST_CHECK(lSurfaces.size() == 1u);
  BOOST_CHECK(assignToAll.size() == 1u);
  BOOST_CHECK(assignToAll[0u] == 0u);
  // The radius of the newly created suport surface should be 10 out of the
  // maximum
  CHECK_CLOSE_ABS(lSurfaces[0u]->bounds().values()[0u], 120, 1e-3);

  // Clear up for the next test
  lSurfaces.clear();
  assignToAll.clear();
  // Add split surfaces as support to already exisint surfaces
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Cylinder,
      sValues, std::nullopt, 16u);
  BOOST_CHECK(lSurfaces.size() == 16u);
  BOOST_CHECK(assignToAll.empty());
}

BOOST_AUTO_TEST_CASE(addDiscSupport) {
  std::vector<std::shared_ptr<Acts::Surface>> lSurfaces;
  std::vector<size_t> assignToAll;

  // The Extent
  Acts::Extent lExtent;
  lExtent.set(Acts::binR, 100., 400.);
  lExtent.set(Acts::binZ, -110., -100.);

  // Get the main support parameters:
  // - doff .. offset (in r.z)
  // - demin, demax .. envelop min, max (in z,r)
  // - dphimin, dphimin .. envelop min, max (in phi)
  std::array<Acts::ActsScalar, 5u> sValues = {-10, 10, 20, 0., 0.};
  // Add a single disc as a support cylinder
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Disc,
      sValues, std::nullopt, 1u);
  BOOST_CHECK(lSurfaces.size() == 1u);
  BOOST_CHECK(assignToAll.size() == 1u);
  BOOST_CHECK(assignToAll[0u] == 0u);
  // The radius of the newly created suport surface should be 10 out of the
  // minimum
  CHECK_CLOSE_ABS(lSurfaces[0u]->transform(tContext).translation().z(), -120,
                  1e-3);

  // Clear up for the next test
  lSurfaces.clear();
  assignToAll.clear();
  // Add split surfaces as support disc
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Disc,
      sValues, std::nullopt, 16u);
  BOOST_CHECK(lSurfaces.size() == 16u);
  BOOST_CHECK(assignToAll.empty());
}

BOOST_AUTO_TEST_CASE(addMisconfiguredSupport) {
  std::vector<std::shared_ptr<Acts::Surface>> lSurfaces;
  std::vector<size_t> assignToAll;

  // Get the main support parameters:
  // - doff .. offset (in r.z)
  // - demin, demax .. envelop min, max (in z,r)
  // - dphimin, dphimin .. envelop min, max (in phi)
  std::array<Acts::ActsScalar, 5u> sValues = {-10, 10, 20, 0., 0.};

  // Unconstrainted extent
  Acts::Extent lExtent;

  // R - Z are not constrained
  // Add a single disc as a support cylinder
  BOOST_CHECK_THROW(
      Acts::Experimental::detail::SupportHelper::addSupport(
          lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Cylinder,
          sValues, std::nullopt, 1u),
      std::runtime_error);

  // The Extent
  lExtent.set(Acts::binR, 100., 400.);
  lExtent.set(Acts::binZ, -110., -100.);

  // Add a single disc as a support cylinder
  BOOST_CHECK_THROW(
      Acts::Experimental::detail::SupportHelper::addSupport(
          lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Cone,
          sValues, std::nullopt, 1u),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
