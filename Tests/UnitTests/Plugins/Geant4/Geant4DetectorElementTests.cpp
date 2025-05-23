// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/AlignmentDelegate.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <utility>

/// Mockup class
class G4VPhysicalVolume {};

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Geant4Plugin)

BOOST_AUTO_TEST_CASE(Geant4DetectorElementConstruction) {
  // G4 phys volume
  auto g4physVol = std::make_shared<G4VPhysicalVolume>();

  // A surface
  auto rBounds = std::make_shared<Acts::RectangleBounds>(10., 10.);
  auto rTransform = Acts::Transform3::Identity();
  rTransform.pretranslate(Acts::Vector3(0., 0., 10));
  auto rSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      rTransform, std::move(rBounds));
  // A detector element
  Acts::Geant4DetectorElement g4DetElement(rSurface, *g4physVol, rTransform,
                                           0.1);

  BOOST_CHECK_EQUAL(g4DetElement.thickness(), 0.1);
  BOOST_CHECK_EQUAL(&g4DetElement.surface(), rSurface.get());
  BOOST_CHECK_EQUAL(&g4DetElement.g4PhysicalVolume(), g4physVol.get());
  BOOST_CHECK(
      g4DetElement.transform(tContext).isApprox(rSurface->transform(tContext)));
}

BOOST_AUTO_TEST_CASE(Geant4DetectorElementAlginment) {
  // G4 phys volume
  auto g4physVol = std::make_shared<G4VPhysicalVolume>();

  // A surface
  auto rBounds = std::make_shared<Acts::RectangleBounds>(10., 10.);
  auto rTransform = Acts::Transform3::Identity();
  rTransform.pretranslate(Acts::Vector3(0., 0., 10));

  auto rSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      rTransform, std::move(rBounds));
  // A detector element
  Acts::Geant4DetectorElement g4DetElement(rSurface, *g4physVol, rTransform,
                                           0.1);

  auto alignedTransform = Acts::Transform3::Identity();
  alignedTransform.pretranslate(Acts::Vector3(1., 2., 3.));

  // Test nominal alignment context
  Acts::GeometryContext nominalContext;
  BOOST_CHECK(
      g4DetElement.surface().transform(nominalContext).isApprox(rTransform));

  // Test aligned context via AlignmentDelegate
  // Mockup struct for the delegate
  struct AlignedTransform {
    Acts::Transform3 data = Acts::Transform3::Identity();

    const Acts::Transform3* provide(
        [[maybe_unused]] const Acts::Surface& surface) const {
      return &data;
    }
  };

  AlignedTransform aTransform{alignedTransform};
  // Create and connect the delegate, and finally the GeometryContext
  Acts::AlignmentDelegate alignmentDelegate;
  alignmentDelegate.connect<&AlignedTransform::provide>(&aTransform);
  Acts::GeometryContext alignedContext(alignmentDelegate);

  BOOST_CHECK(g4DetElement.surface()
                  .transform(alignedContext)
                  .isApprox(alignedTransform));
}

BOOST_AUTO_TEST_SUITE_END()
