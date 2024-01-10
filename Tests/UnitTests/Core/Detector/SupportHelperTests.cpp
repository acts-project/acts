// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/SupportHelper.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

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

BOOST_AUTO_TEST_CASE(RectangularSupport) {
  // As a single rectangle
  auto rectangleSupport =
      Acts::Experimental::detail::SupportHelper::rectangularSupport(
          Acts::Transform3::Identity(), {100., 400});
  BOOST_CHECK_EQUAL(rectangleSupport.size(), 1u);
  BOOST_CHECK_EQUAL(rectangleSupport[0u]->type(),
              Acts::Surface::SurfaceType::Plane);
}

BOOST_AUTO_TEST_CASE(addCylinderSupport) {
  std::vector<std::shared_ptr<Acts::Surface>> lSurfaces;
  std::vector<size_t> assignToAll;

  // The Extent
  Acts::Extent lExtent;
  lExtent.set(Acts::binR, 100., 110.);
  lExtent.set(Acts::binZ, -400., 400.);

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
  // The radius of the newly created support surface should be 10 out of the
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

  // Special test to check that support can be smaller than extent
  // Clear up for the next test
  sValues = {10, -5, -5, 0., 0.};
  lSurfaces.clear();
  assignToAll.clear();
  // Add split surfaces as support to already exisint surfaces
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Cylinder,
      sValues, std::nullopt, 1u);
  BOOST_CHECK_EQUAL(lSurfaces.size(), 1u);
  auto surface = lSurfaces[0u];

  auto boundValues = surface->bounds().values();
  CHECK_CLOSE_ABS(boundValues[1u], 395, 1e-3);

  // Invalid test catch
  sValues = {10, -500, -600, 0., 0.};
  BOOST_CHECK_THROW(
      Acts::Experimental::detail::SupportHelper::addSupport(
          lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Cylinder,
          sValues, std::nullopt, 1u),
      std::runtime_error);
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
  // The radius of the newly created support surface should be 10 out of the
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

  // Special test that the support can be smaller than the estimated layer
  // structure
  sValues = {-10, -10, -20, 0., 0.};

  lSurfaces.clear();
  assignToAll.clear();
  // Add split surfaces as support disc
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Disc,
      sValues, std::nullopt, 1u);

  BOOST_CHECK(lSurfaces.size() == 1u);
  auto surface = lSurfaces[0u];

  auto boundValues = surface->bounds().values();
  CHECK_CLOSE_ABS(boundValues[0u], 110, 1e-3);
  CHECK_CLOSE_ABS(boundValues[1u], 380, 1e-3);

  // Invalid test catch
  sValues = {-10, -500, -600, 0., 0.};

  BOOST_CHECK_THROW(
      Acts::Experimental::detail::SupportHelper::addSupport(
          lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Disc,
          sValues, std::nullopt, 1u),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(addRectangularSupport) {
  // Possible binning iterations
  using BinningSet = std::array<Acts::BinningValue, 3u>;

  std::vector<BinningSet> binningSets = {
      BinningSet{Acts::binX, Acts::binY, Acts::binZ},
      BinningSet{Acts::binY, Acts::binZ, Acts::binX},
      BinningSet{Acts::binZ, Acts::binX, Acts::binY}};

  for (auto [sBinning, l0Binning, l1Binning] : binningSets) {
    std::vector<std::shared_ptr<Acts::Surface>> lSurfaces;
    std::vector<size_t> assignToAll;

    // The Extent
    Acts::Extent lExtent;
    lExtent.set(sBinning, -20., 20.);
    lExtent.set(l0Binning, -300., 300.);
    lExtent.set(l1Binning, -200., 200.);

    // The support parameters:
    std::array<Acts::ActsScalar, 5u> sValues = {-10, 10, 10, 20., 20.};

    // Add a single rectangle plane as a support rectantle
    Acts::Experimental::detail::SupportHelper::addSupport(
        lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Plane,
        sValues, std::nullopt, 1u, sBinning);

    BOOST_CHECK(lSurfaces.size() == 1u);
    BOOST_CHECK(assignToAll.size() == 1u);
    BOOST_CHECK(assignToAll[0u] == 0u);

    auto surface = lSurfaces[0u];
    auto boundValues = surface->bounds().values();
    CHECK_CLOSE_ABS(boundValues[0u], -310, 1e-3);
    CHECK_CLOSE_ABS(boundValues[1u], -220, 1e-3);
    CHECK_CLOSE_ABS(boundValues[2u], 310, 1e-3);
    CHECK_CLOSE_ABS(boundValues[3u], 220, 1e-3);

    CHECK_CLOSE_ABS(surface->transform(tContext).translation()[sBinning], -10,
                    1e-3);
  }

  // Dedicated test that the support can also be forced to be smaller
  std::vector<std::shared_ptr<Acts::Surface>> lSurfaces;
  std::vector<size_t> assignToAll;

  // The Extent
  Acts::Extent lExtent;
  lExtent.set(Acts::binZ, -20., 20.);
  lExtent.set(Acts::binX, -300., 300.);
  lExtent.set(Acts::binY, -200., 200.);

  // The support parameters:
  std::array<Acts::ActsScalar, 5u> sValues = {-10, -10, -10, -20., -20.};

  // Add a single rectangle plane as a support cylinder
  Acts::Experimental::detail::SupportHelper::addSupport(
      lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Plane,
      sValues, std::nullopt, 1u, Acts::binZ);

  auto surface = lSurfaces[0u];
  auto boundValues = surface->bounds().values();
  CHECK_CLOSE_ABS(boundValues[0u], -290, 1e-3);
  CHECK_CLOSE_ABS(boundValues[1u], -180, 1e-3);
  CHECK_CLOSE_ABS(boundValues[2u], 290, 1e-3);
  CHECK_CLOSE_ABS(boundValues[3u], 180, 1e-3);

  // Catch invalid test
  sValues = {-10, -500, -600, 0., 0.};

  BOOST_CHECK_THROW(
      Acts::Experimental::detail::SupportHelper::addSupport(
          lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Plane,
          sValues, std::nullopt, 1u, Acts::binZ),
      std::runtime_error);
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

  // Non-constrained extent for rectangular support
  BOOST_CHECK_THROW(
      Acts::Experimental::detail::SupportHelper::addSupport(
          lSurfaces, assignToAll, lExtent, Acts::Surface::SurfaceType::Plane,
          sValues, std::nullopt, 1u, Acts::binZ),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
