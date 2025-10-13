// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/detail/PortalHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/Json/PortalJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>

#include <nlohmann/json.hpp>

namespace Acts::Experimental {
class DetectorVolume {};
}  // namespace Acts::Experimental

using namespace Acts;

GeometryContext tContext;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(PortalUnconnected) {
  std::ofstream out;

  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3(0., 0., 0.), Vector3(0., 1., 0.))
          .planeSurface();

  auto portal = std::make_shared<Experimental::Portal>(std::move(surface));

  BOOST_CHECK_NE(portal, nullptr);

  auto jPortal = PortalJsonConverter::toJson(tContext, *portal, {});

  out.open("portal.json");
  out << jPortal.dump(4);
  out.close();

  // Now read it back in
  auto in =
      std::ifstream("portal.json", std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json jPortalIn;
  in >> jPortalIn;
  in.close();

  auto portalIn = PortalJsonConverter::fromJson(tContext, jPortalIn, {});

  BOOST_CHECK_NE(portalIn, nullptr);
}

BOOST_AUTO_TEST_CASE(PortalSingleConnected) {
  std::ofstream out;

  auto forwardVolume = std::make_shared<Experimental::DetectorVolume>();
  auto backwardVolume = std::make_shared<Experimental::DetectorVolume>();

  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3(0., 0., 0.), Vector3(0., 1., 0.))
          .planeSurface();

  auto portal = std::make_shared<Experimental::Portal>(std::move(surface));
  BOOST_CHECK_NE(portal, nullptr);
  // Attaching the portals
  Experimental::detail::PortalHelper::attachExternalNavigationDelegate(
      *portal, forwardVolume, Direction::Forward());
  Experimental::detail::PortalHelper::attachExternalNavigationDelegate(
      *portal, backwardVolume, Direction::Backward());

  std::vector<const Experimental::DetectorVolume*> detectorVolumes = {
      forwardVolume.get(), backwardVolume.get()};
  // No volumes provided, must bail
  BOOST_CHECK_THROW(PortalJsonConverter::toJson(tContext, *portal, {}),
                    std::runtime_error);
  auto jPortal =
      PortalJsonConverter::toJson(tContext, *portal, detectorVolumes);

  out.open("portal-single-connected.json");
  out << jPortal.dump(4);
  out.close();

  // Now read it back in
  auto in = std::ifstream("portal-single-connected.json",
                          std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json jPortalIn;
  in >> jPortalIn;
  in.close();

  auto portalIn = PortalJsonConverter::fromJson(
      tContext, jPortalIn, {forwardVolume, backwardVolume});
  BOOST_CHECK_NE(portalIn, nullptr);
}

BOOST_AUTO_TEST_CASE(PortalMultiConnected) {
  std::ofstream out;

  auto forwardVolumeA = std::make_shared<Experimental::DetectorVolume>();
  auto forwardVolumeB = std::make_shared<Experimental::DetectorVolume>();
  auto forwardVolumeC = std::make_shared<Experimental::DetectorVolume>();

  auto backwardVolume = std::make_shared<Experimental::DetectorVolume>();

  std::shared_ptr<PlaneSurface> surface =
      CurvilinearSurface(Vector3(0., 0., 0.), Vector3(0., 1., 0.))
          .planeSurface();

  auto portal = std::make_shared<Experimental::Portal>(std::move(surface));
  BOOST_CHECK_NE(portal, nullptr);

  // Attaching the portals
  Experimental::detail::PortalHelper::attachExternalNavigationDelegate(
      *portal, backwardVolume, Direction::Backward());

  Experimental::detail::PortalHelper::attachDetectorVolumesUpdater(
      tContext, *portal, {forwardVolumeA, forwardVolumeB, forwardVolumeC},
      Direction::Forward(), {-100, 10, 20, 200}, AxisDirection::AxisX);

  std::vector<const Experimental::DetectorVolume*> detectorVolumes = {
      forwardVolumeA.get(), forwardVolumeB.get(), forwardVolumeC.get(),
      backwardVolume.get()};

  auto jPortal =
      PortalJsonConverter::toJson(tContext, *portal, detectorVolumes);

  out.open("portal-multi-connected.json");
  out << jPortal.dump(4);
  out.close();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
