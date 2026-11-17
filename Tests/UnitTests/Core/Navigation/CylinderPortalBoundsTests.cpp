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
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/detail/CylinderPortalBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <array>
#include <memory>
#include <numbers>
#include <utility>

namespace ActsTests {
namespace {
using namespace Acts;
using enum CylinderVolumeBounds::Face;

struct Fixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  std::unique_ptr<TrackingVolume> annulus(
      const Transform3& transform = Transform3::Identity(),
      double halfPhi = std::numbers::pi) const {
    return std::make_unique<TrackingVolume>(
        transform,
        std::make_shared<CylinderVolumeBounds>(100., 400., 300., halfPhi),
        "Annulus");
  }
};
}  // namespace

BOOST_FIXTURE_TEST_SUITE(CylinderPortalBounds, Fixture)

BOOST_AUTO_TEST_CASE(AnnulusAndStreamLifetime) {
  auto volume = annulus();
  SingleCylinderPortalShell shell(gctx, *volume);
  shell.applyToVolume();
  detail::CylinderPortalBounds bounds(gctx, *volume);

  BOOST_REQUIRE(bounds.valid());
  BOOST_CHECK_EQUAL(bounds.cachedVolume, volume.get());
  BOOST_CHECK_EQUAL(bounds.innerCylinder, shell.portal(InnerCylinder).get());
  BOOST_CHECK_EQUAL(bounds.outerCylinder, shell.portal(OuterCylinder).get());
  BOOST_CHECK_EQUAL(bounds.lowerDisc, shell.portal(NegativeDisc).get());
  BOOST_CHECK_EQUAL(bounds.upperDisc, shell.portal(PositiveDisc).get());
  BOOST_CHECK_EQUAL(bounds.innerRadius2, 100. * 100.);
  BOOST_CHECK_EQUAL(bounds.outerRadius2, 400. * 400.);
  BOOST_CHECK_EQUAL(bounds.lowerZ, -300.);
  BOOST_CHECK_EQUAL(bounds.upperZ, 300.);
  BOOST_CHECK(!bounds.nestedContainer);

  NavigationStream stream;
  AppendOnlyNavigationStream append(stream);
  BOOST_CHECK(detail::CylinderPortalBounds::get(stream) == nullptr);
  bounds.attach(append);
  BOOST_CHECK_EQUAL(detail::CylinderPortalBounds::get(stream), &bounds);
  stream.reset();
  BOOST_CHECK(detail::CylinderPortalBounds::get(stream) == nullptr);
}

BOOST_AUTO_TEST_CASE(AxialTranslationsAndRotationsRemainSupported) {
  for (const Vector3& offset :
       std::array<Vector3, 3>{Vector3(0., 0., 50.), Vector3(20., -30., 0.),
                              Vector3(20., -30., 50.)}) {
    Transform3 transform = Transform3::Identity();
    transform.translate(offset);
    transform.rotate(AngleAxis3(0.7, Vector3::UnitZ()));
    auto volume = annulus(transform);
    SingleCylinderPortalShell shell(gctx, *volume);
    shell.applyToVolume();
    detail::CylinderPortalBounds bounds(gctx, *volume);

    BOOST_REQUIRE(bounds.valid());
    BOOST_CHECK_EQUAL(bounds.center, offset.head<2>());
    BOOST_CHECK_EQUAL(bounds.lowerZ, offset.z() - 300.);
    BOOST_CHECK_EQUAL(bounds.upperZ, offset.z() + 300.);
  }
}

BOOST_AUTO_TEST_CASE(TiltedVolumeIsUnsupported) {
  Transform3 transform = Transform3::Identity();
  transform.rotate(AngleAxis3(0.2, Vector3::UnitX()));
  auto volume = annulus(transform);
  SingleCylinderPortalShell shell(gctx, *volume);
  shell.applyToVolume();
  detail::CylinderPortalBounds bounds(gctx, *volume);
  BOOST_CHECK(!bounds.valid());
}

BOOST_AUTO_TEST_CASE(SectorVolumeIsUnsupported) {
  auto volume = annulus(Transform3::Identity(), std::numbers::pi / 2.);
  SingleCylinderPortalShell shell(gctx, *volume);
  shell.applyToVolume();
  BOOST_REQUIRE_EQUAL(volume->portals().size(), 6u);
  detail::CylinderPortalBounds bounds(gctx, *volume);
  BOOST_CHECK(!bounds.valid());

  NavigationStream stream;
  AppendOnlyNavigationStream append(stream);
  bounds.attach(append);
  BOOST_CHECK(detail::CylinderPortalBounds::get(stream) == nullptr);
}

BOOST_AUTO_TEST_CASE(RejectUnsupportedPortalSupports) {
  // Replacing one portal leaves the owner's bounds and portal count intact.
  // Each support must therefore be validated independently of the owner.
  for (bool cylinder : {false, true}) {
    for (unsigned int defect = 0; defect < 3; ++defect) {
      BOOST_TEST_CONTEXT("cylinder=" << cylinder << ", defect=" << defect) {
        auto volume = annulus();
        SingleCylinderPortalShell shell(gctx, *volume);
        Transform3 transform = Transform3::Identity();
        if (!cylinder) {
          transform.translate(Vector3(0., 0., 300.));
        }
        if (defect == 0) {
          transform.translate(Vector3(10., 0., 0.));
        } else if (defect == 1) {
          transform.rotate(AngleAxis3(0.2, Vector3::UnitX()));
        }
        const double halfPhi =
            defect == 2 ? std::numbers::pi / 2. : std::numbers::pi;
        std::shared_ptr<RegularSurface> surface;
        if (cylinder) {
          surface = Surface::makeShared<CylinderSurface>(transform, 100., 300.,
                                                         halfPhi);
        } else {
          surface =
              Surface::makeShared<DiscSurface>(transform, 100., 400., halfPhi);
        }
        auto portal = std::make_shared<Portal>(Direction::AlongNormal(),
                                               std::move(surface), *volume);
        shell.setPortal(std::move(portal),
                        cylinder ? InnerCylinder : PositiveDisc);
        shell.applyToVolume();
        BOOST_REQUIRE_EQUAL(volume->portals().size(), 4u);
        detail::CylinderPortalBounds bounds(gctx, *volume);
        BOOST_CHECK(!bounds.valid());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(NestedCylinderHasFourDiscInterfaces) {
  TrackingVolume parent(Transform3::Identity(),
                        std::make_shared<CylinderVolumeBounds>(0., 400., 300.),
                        "Parent");
  Transform3 childTransform = Transform3::Identity();
  childTransform.translate(Vector3(0., 0., 50.));
  TrackingVolume child(childTransform,
                       std::make_shared<CylinderVolumeBounds>(0., 100., 100.),
                       "Child");
  SingleCylinderPortalShell parentShell(gctx, parent);
  SingleCylinderPortalShell childShell(gctx, child);
  parentShell.applyToVolume();
  // A container sees both its own outer faces and its child's three faces.
  for (auto face : {NegativeDisc, PositiveDisc, OuterCylinder}) {
    parent.addPortal(childShell.portal(face));
  }
  BOOST_REQUIRE_EQUAL(parent.portals().size(), 6u);
  detail::CylinderPortalBounds bounds(gctx, parent);
  BOOST_REQUIRE(bounds.valid());
  BOOST_CHECK(bounds.nestedContainer);
  BOOST_CHECK_EQUAL(bounds.containerDiscCount, 4u);
  BOOST_CHECK_EQUAL(bounds.innerCylinder,
                    childShell.portal(OuterCylinder).get());
  BOOST_CHECK_EQUAL(bounds.innerMinZ, -50.);
  BOOST_CHECK_EQUAL(bounds.innerMaxZ, 150.);
}

BOOST_AUTO_TEST_CASE(PortalCountAloneDoesNotIdentifyAContainer) {
  auto volume = annulus();
  SingleCylinderPortalShell shell(gctx, *volume);
  shell.applyToVolume();
  for (unsigned int i = 0; i < 2; ++i) {
    volume->addPortal(std::make_shared<Portal>(
        Direction::AlongNormal(),
        Surface::makeShared<PlaneSurface>(Transform3::Identity()), *volume));
  }
  BOOST_REQUIRE_EQUAL(volume->portals().size(), 6u);
  detail::CylinderPortalBounds bounds(gctx, *volume);
  BOOST_CHECK(!bounds.valid());
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
