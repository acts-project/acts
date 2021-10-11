// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <vector>

#include "GeometryHelper.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Experimental)

/// This test tests the structure and consistency
/// of the internal object store of the DetectorVolume
BOOST_AUTO_TEST_CASE(DetectorVolumeObjectStore) {
  /// An object struct
  struct A {
    int i = 0;
    A(int _i) : i(_i) {}
    A() = default;
  };

  // Test for shared vector
  std::vector<std::shared_ptr<A>> aSharedVector = {
      std::make_shared<A>(1), std::make_shared<A>(2), std::make_shared<A>(3)};

  DetectorVolume::ObjectStore<std::shared_ptr<A>> aSharedStore(
      std::move(aSharedVector));
  // Check size and pointer consistency
  BOOST_CHECK(aSharedStore.internal.size() == aSharedStore.external.size());
  for (auto [i, refobject] : enumerate(aSharedStore.internal)) {
    BOOST_CHECK(refobject.get() == aSharedStore.external[i]);
  }
}

BOOST_AUTO_TEST_CASE(DetectorVolumeInside) {
  // Create some unique Bounds
  auto cylVolBounds = std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 10.);

  Transform3 translated = Transform3::Identity();
  translated.pretranslate(Vector3{0.5, 0., 25.});

  auto tranlatedVolume = DetectorVolume::makeShared(
      translated, std::move(cylVolBounds), "TranslatedVolume");

  /// Origin (0, 0, 0) is outside
  BOOST_CHECK(not tranlatedVolume->inside(Vector3{0., 0., 0.}));
  // Just outside
  BOOST_CHECK(not tranlatedVolume->inside(Vector3{1.49, 0., 30.}));
  // Inside because of tolerance
  BOOST_CHECK(tranlatedVolume->inside(Vector3{1.49, 0., 30.}, 0.2));
  // Properly inside
  BOOST_CHECK(tranlatedVolume->inside(Vector3{0.5, 3., 20}));
}

BOOST_AUTO_TEST_CASE(DetectorVolumePortalAccess) {
  // First cylindrical volumes
  auto cylVol = std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 10.);
  Transform3 transform = Transform3::Identity();
  auto volume =
      DetectorVolume::makeShared(transform, std::move(cylVol), "Volume");

  // Access portals (cont access), let's copy
  std::vector<const Portal*> portals = volume->portals();
  BOOST_TEST(portals.size(), 4u);

  auto portalPtrs = volume->portalPtrs();
  BOOST_TEST(portalPtrs.size(), 4u);

  // Raw pointers are identical
  for (auto [i, p] : enumerate(portals)) {
    BOOST_CHECK(p == portalPtrs[i].get());
  }

  // Now let us exchange the portal disc at positive z
  Transform3 pDiscTransform = Transform3::Identity();
  pDiscTransform.pretranslate(Vector3{0., 0., 10.});
  auto pDiscBounds = std::make_shared<RadialBounds>(1., 4.);
  auto pDisc = Surface::makeShared<DiscSurface>(pDiscTransform, pDiscBounds);
  auto pPortal = std::make_shared<Portal>(std::move(pDisc));
  volume->updatePortalPtr(pPortal, 1);

  auto uPortals = volume->portals();
  BOOST_TEST(uPortals.size(), 4u);
  // Test against all reference
  BOOST_CHECK(uPortals[0] == portals[0]);
  BOOST_CHECK(uPortals[1] != portals[1]);
  BOOST_CHECK(uPortals[2] == portals[2]);
  BOOST_CHECK(uPortals[3] == portals[3]);

  auto uPortalPtrs = volume->portalPtrs();
  BOOST_TEST(uPortalPtrs.size(), 4u);

  // Raw pointers are identical
  for (auto [i, p] : enumerate(uPortals)) {
    BOOST_CHECK(p == uPortalPtrs[i].get());
  }
}

BOOST_AUTO_TEST_CASE(DetectorVolumeCylContainerInZ) {
  // Three cylindrical volumes
  auto cylVol0 = std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 10.);
  auto cylVol1 = std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 5.);
  auto cylVol2 = std::make_unique<Acts::CylinderVolumeBounds>(1., 4., 10.);

  Transform3 transform0 = Transform3::Identity();
  transform0.pretranslate(Vector3{0., 0., -15.});
  Transform3 transform1 = Transform3::Identity();
  transform1.pretranslate(Vector3{0., 0., 0.});
  Transform3 transform2 = Transform3::Identity();
  transform2.pretranslate(Vector3{0., 0., 15.});

  // Create the three volumes ...
  auto volume0 =
      DetectorVolume::makeShared(transform0, std::move(cylVol0), "Volume0");
  auto volume1 =
      DetectorVolume::makeShared(transform1, std::move(cylVol1), "Volume1");
  auto volume2 =
      DetectorVolume::makeShared(transform2, std::move(cylVol2), "Volume2");

  std::vector<std::shared_ptr<DetectorVolume>> volumes = {volume0, volume1,
                                                          volume2};
}

BOOST_AUTO_TEST_CASE(DetectorVolumeEnvironment) {
  // Create a single barrel volume with surfaces
  auto volumeWithSurfaces = createBarrelVolume(20.);
  GeometryContext gctx = GeometryContext();

  // Start inside the volume
  Vector3 position(22., 0., 0.);
  Vector3 direction = Vector3(0.75, 0.75, 0.1).normalized();
  BoundaryCheck bCheck = true;

  auto nEnvironment =
      volumeWithSurfaces->environment(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == volumeWithSurfaces.get());
  // There should be 4 portals to try
  BOOST_CHECK(nEnvironment.portals.size() == 4u);
  // There should be 1 surface accessible (empirically)
  BOOST_CHECK(nEnvironment.surfaces.size() == 3u);
}

BOOST_AUTO_TEST_CASE(DetectorVolumePortalEntryExit) {
  // Create a single barrel volume with surfaces
  ActsScalar vRmin = 20.;
  ActsScalar vRmax = 60.;
  ActsScalar vZhalf = 500.;
  // Create the volume now
  auto volumeWithSurfaces = createBarrelVolume(vRmin, vRmax, vZhalf);
  GeometryContext gctx = GeometryContext();

  const auto& portals = volumeWithSurfaces->portals();

  // ENTRY TESTS

  // Entry through negative disc
  const Portal* negativeDiscPortal = portals[0];

  // Start on the disc
  Vector3 position(vRmin + 2., 0., -vZhalf);
  Vector3 direction = Vector3(0.25, 0.15, 1.0).normalized();
  BoundaryCheck bCheck = true;

  // Get the new environement
  auto nEnvironment =
      negativeDiscPortal->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == volumeWithSurfaces.get());
  // There should be 4 portals to try
  BOOST_CHECK(nEnvironment.portals.size() == 4u);
  // There should be 1 surface accessible (empirically)
  BOOST_CHECK(nEnvironment.surfaces.size() == 1u);

  // Entry through positive disc
  const Portal* positiveDiscPortal = portals[1];

  // Start on the disc
  position = Vector3(vRmin + 2, 0., vZhalf);
  direction = Vector3(0.25, 0.15, -1.0).normalized();
  nEnvironment = positiveDiscPortal->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == volumeWithSurfaces.get());
  // There should be 4 portals to try
  BOOST_CHECK(nEnvironment.portals.size() == 4u);
  // There should be 1 surface accessible (empirically)
  BOOST_CHECK(nEnvironment.surfaces.size() == 1u);

  // Entry through the inner tube
  const Portal* innerCylinder = portals[3];

  // Start on the inner cylinder
  position = Vector3(vRmin, 0., 5.);
  direction = Vector3(0.8, 0.4, 0.05).normalized();
  nEnvironment = innerCylinder->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == volumeWithSurfaces.get());
  // There should be 4 portals to try
  BOOST_CHECK(nEnvironment.portals.size() == 4u);
  // There should be 4 surfaces accessible (empirically)
  BOOST_CHECK(nEnvironment.surfaces.size() == 4u);

  // Start on the outer cylinder
  const Portal* outerCylinder = portals[2];
  // Start on the inner cylinder
  position = Vector3(vRmax, 0., 5.);
  direction = Vector3(-0.8, -0.4, 0.05).normalized();
  nEnvironment = outerCylinder->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == volumeWithSurfaces.get());
  // There should be 4 portals to try
  BOOST_CHECK(nEnvironment.portals.size() == 4u);
  // There should be 2 surface accessible (empirically)
  BOOST_CHECK(nEnvironment.surfaces.size() == 2u);

  // EXIT TESTS
  position = Vector3(vRmin + 2., 0., -vZhalf);
  direction = Vector3(0.25, 0.15, -1.0).normalized();

  // Get the new environement
  nEnvironment = negativeDiscPortal->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == nullptr);
  BOOST_CHECK(nEnvironment.portals.empty());
  BOOST_CHECK(nEnvironment.surfaces.empty());

  position = Vector3(vRmin + 2., 0., vZhalf);
  direction = Vector3(0.25, 0.15, 1.0).normalized();

  // Get the new environement
  nEnvironment = positiveDiscPortal->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == nullptr);
  BOOST_CHECK(nEnvironment.portals.empty());
  BOOST_CHECK(nEnvironment.surfaces.empty());

  // Start on the inner cylinder
  position = Vector3(vRmin, 0., 5.);
  direction = Vector3(-0.8, -0.4, 0.05).normalized();
  nEnvironment = innerCylinder->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == nullptr);
  BOOST_CHECK(nEnvironment.portals.empty());
  BOOST_CHECK(nEnvironment.surfaces.empty());

  // Start on the inner cylinder
  position = Vector3(vRmax, 0., 5.);
  direction = Vector3(0.8, 0.4, 0.05).normalized();
  nEnvironment = outerCylinder->next(gctx, position, direction, bCheck);
  BOOST_CHECK(nEnvironment.volume == nullptr);
  BOOST_CHECK(nEnvironment.portals.empty());
  BOOST_CHECK(nEnvironment.surfaces.empty());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
