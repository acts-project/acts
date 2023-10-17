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
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/// Unpack to shared - simply to test the getSharedPtr mechanism
///
/// @tparam referenced_type is the type of the referenced object
///
/// @param rt is the referenced object
///
/// @returns a shared pointer
template <typename referenced_type>
std::shared_ptr<referenced_type> unpackToShared(referenced_type& rt) {
  return rt.getSharedPtr();
}

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(DetectorConstruction) {
  Acts::ActsScalar r0 = 0.;
  Acts::ActsScalar r1 = 10.;
  Acts::ActsScalar r2 = 100.;
  Acts::ActsScalar r3 = 200.;
  Acts::ActsScalar zHalfL = 200.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  // Create a bunch of volumes
  auto cyl0Bounds =
      std::make_unique<Acts::CylinderVolumeBounds>(r0, r1, zHalfL);

  auto cyl0BoundsCopy =
      std::make_unique<Acts::CylinderVolumeBounds>(r0, r1, zHalfL);

  auto cyl1Bounds =
      std::make_unique<Acts::CylinderVolumeBounds>(r1, r2, zHalfL);

  auto cyl2Bounds =
      std::make_unique<Acts::CylinderVolumeBounds>(r2, r3, zHalfL);

  auto portalGenerator = Acts::Experimental::defaultPortalGenerator();

  auto cyl0 = Acts::Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl0", nominal, std::move(cyl0Bounds),
      Acts::Experimental::tryAllPortals());

  auto cyl0nameDup = Acts::Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl0", nominal, std::move(cyl0BoundsCopy),
      Acts::Experimental::tryAllPortals());

  auto cyl1 = Acts::Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl1", nominal, std::move(cyl1Bounds),
      Acts::Experimental::tryAllPortals());

  auto cyl2 = Acts::Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "Cyl2", nominal, std::move(cyl2Bounds),
      Acts::Experimental::tryAllPortals());

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes012 =
      {cyl0, cyl1, cyl2};
  auto det012 = Acts::Experimental::Detector::makeShared(
      "Det012", volumes012, Acts::Experimental::tryRootVolumes());

  // Check the basic return functions
  BOOST_CHECK(det012->name() == "Det012");
  BOOST_CHECK(det012->volumes().size() == 3u);
  BOOST_CHECK(det012->volumePtrs().size() == 3u);

  // Check the shared pointer mechanism
  BOOST_CHECK(det012 == unpackToShared<Acts::Experimental::Detector>(*det012));
  BOOST_CHECK(det012 ==
              unpackToShared<const Acts::Experimental::Detector>(*det012));

  // Check the inside function with positions
  Acts::Experimental::NavigationState nState;
  nState.position = Acts::Vector3(5., 0., 0.);
  nState.currentDetector = det012.get();
  det012->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == cyl0.get());

  auto find1 = det012->findDetectorVolume(tContext, Acts::Vector3(15., 0., 0.));
  BOOST_CHECK(find1 == cyl1.get());

  auto find2 =
      det012->findDetectorVolume(tContext, Acts::Vector3(150., 0., 0.));
  BOOST_CHECK(find2 == cyl2.get());

  auto findNull =
      det012->findDetectorVolume(tContext, Acts::Vector3(1500., 0., 0.));
  BOOST_CHECK(findNull == nullptr);

  /// Find by name
  auto find0 = det012->findDetectorVolume("Cyl0");
  BOOST_CHECK(find0 == cyl0.get());

  findNull = det012->findDetectorVolume("Null");
  BOOST_CHECK(findNull == nullptr);

  // Misconfigured - unkonnected finder
  Acts::Experimental::DetectorVolumeUpdator unconnected;
  BOOST_CHECK_THROW(
      Acts::Experimental::Detector::makeShared("Det012_unconnected", volumes012,
                                               std::move(unconnected)),
      std::invalid_argument);

  // Misconfigured - duplicate name
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes002 =
      {cyl0, cyl0nameDup, cyl2};
  BOOST_CHECK_THROW(Acts::Experimental::Detector::makeShared(
                        "Det002_name_duplicate", volumes002,
                        Acts::Experimental::tryRootVolumes()),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DetectorConstructionWithHierarchyMap) {
  auto portalGenerator = Acts::Experimental::defaultPortalGenerator();

  std::vector<std::unique_ptr<Acts::Test::DetectorElementStub>> detStore;
  std::vector<Acts::ActsScalar> radii = {100, 102, 104, 106, 108, 110};
  auto cylinderVoumeBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(80, 130, 200);
  std::vector<std::shared_ptr<Acts::Surface>> surfaces = {};
  for (auto [ir, r] : Acts::enumerate(radii)) {
    auto detElement = std::make_unique<Acts::Test::DetectorElementStub>(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::CylinderBounds>(r, 190.), 0.1);
    auto surface = detElement->surface().getSharedPtr();
    surface->assignGeometryId(Acts::GeometryIdentifier{}.setSensitive(ir + 1));
    surfaces.push_back(std::move(surface));
    detStore.push_back(std::move(detElement));
  }

  auto cylVolume = Acts::Experimental::DetectorVolumeFactory::construct(
      portalGenerator, tContext, "CylinderVolume", Acts::Transform3::Identity(),
      std::move(cylinderVoumeBounds), surfaces, {},
      Acts::Experimental::tryNoVolumes(),
      Acts::Experimental::tryAllPortalsAndSurfaces());

  auto det = Acts::Experimental::Detector::makeShared(
      "DetWithSurfaces", {cylVolume}, Acts::Experimental::tryRootVolumes());

  const auto& sensitiveHierarchyMap = det->sensitiveHierarchyMap();
  BOOST_CHECK(sensitiveHierarchyMap.size() == 6u);
}

BOOST_AUTO_TEST_SUITE_END()
