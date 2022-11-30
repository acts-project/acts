// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Geant4Plugin)

BOOST_AUTO_TEST_CASE(Geant4AlgebraConversion) {
  G4ThreeVector g4Translation(10., 20., 30.);

  auto translated = Acts::Geant4AlgebraConverter{}.transform(g4Translation);
  auto actsTranslation = translated.translation();
  BOOST_CHECK(actsTranslation[0] == 10.);
  BOOST_CHECK(actsTranslation[1] == 20.);
  BOOST_CHECK(actsTranslation[2] == 30.);

  auto translatedScaled =
      Acts::Geant4AlgebraConverter{10.}.transform(g4Translation);
  auto actsTranslationScaled = translatedScaled.translation();
  BOOST_CHECK(actsTranslationScaled[0] == 100.);
  BOOST_CHECK(actsTranslationScaled[1] == 200.);
  BOOST_CHECK(actsTranslationScaled[2] == 300.);
}

BOOST_AUTO_TEST_CASE(Geant4CylinderConversion) {
  G4Tubs cylinder("Cylinder", 399., 401., 800., 0., 2 * M_PI);
  auto bounds = Acts::Geant4ShapeConverter{}.cylinderBounds(cylinder);
  CHECK_CLOSE_ABS(bounds->get(Acts::CylinderBounds::BoundValues::eR), 400.,
                  10e-10);
  CHECK_CLOSE_ABS(bounds->get(Acts::CylinderBounds::BoundValues::eHalfLengthZ),
                  800., 10e-10);
  CHECK_CLOSE_ABS(
      bounds->get(Acts::CylinderBounds::BoundValues::eHalfPhiSector), M_PI,
      10e-10);
  CHECK_CLOSE_ABS(bounds->get(Acts::CylinderBounds::BoundValues::eAveragePhi),
                  0., 10e-10);
}

BOOST_AUTO_TEST_CASE(Geant4RadialConversion) {
  G4Tubs disk("disk", 40., 400., 2., 0., 2 * M_PI);
  auto bounds = Acts::Geant4ShapeConverter{}.radialBounds(disk);
  CHECK_CLOSE_ABS(bounds->get(Acts::RadialBounds::BoundValues::eMinR), 40.,
                  10e-10);
  CHECK_CLOSE_ABS(bounds->get(Acts::RadialBounds::BoundValues::eMaxR), 400.,
                  10e-10);
  CHECK_CLOSE_ABS(bounds->get(Acts::RadialBounds::BoundValues::eHalfPhiSector),
                  M_PI, 10e-10);
  CHECK_CLOSE_ABS(bounds->get(Acts::RadialBounds::BoundValues::eAveragePhi), 0.,
                  10e-10);
}

BOOST_AUTO_TEST_CASE(Geant4BoxConversion) {
  // Test the standard orientations
  G4Box sensorXY("SensorXY", 23., 34., 1.);
  auto rectangleXY = Acts::Geant4ShapeConverter{}.rectangleBounds(sensorXY);
  auto boundsXY = std::get<0u>(rectangleXY);
  auto axesXY = std::get<1u>(rectangleXY);
  CHECK_CLOSE_ABS(boundsXY->halfLengthX(), 23., 10e-10);
  CHECK_CLOSE_ABS(boundsXY->halfLengthY(), 34., 10e-10);
  auto refXY = std::array<int, 2u>{0, 1};
  BOOST_CHECK(axesXY == refXY);

  G4Box sensorYZ("SensorYZ", 2., 45., 56.);
  auto rectangleYZ = Acts::Geant4ShapeConverter{}.rectangleBounds(sensorYZ);
  auto boundsYZ = std::get<0u>(rectangleYZ);
  auto axesYZ = std::get<1u>(rectangleYZ);
  CHECK_CLOSE_ABS(boundsYZ->halfLengthX(), 45., 10e-10);
  CHECK_CLOSE_ABS(boundsYZ->halfLengthY(), 56., 10e-10);
  auto refYZ = std::array<int, 2u>{1, 2};
  BOOST_CHECK(axesYZ == refYZ);

  G4Box sensorZX("SensorZX", 78., 2., 67.);
  auto rectangleZX = Acts::Geant4ShapeConverter{}.rectangleBounds(sensorZX);
  auto boundsZX = std::get<0u>(rectangleZX);
  auto axesZX = std::get<1u>(rectangleZX);
  CHECK_CLOSE_ABS(boundsZX->halfLengthX(), 67., 10e-10);
  CHECK_CLOSE_ABS(boundsZX->halfLengthY(), 78., 10e-10);
  auto refZX = std::array<int, 2u>{2, 0};
  BOOST_CHECK(axesZX == refZX);

  // Test the flipped axis
  G4Box sensorXz("SensorXz", 78., 2., 67.);
  auto rectangleXz =
      Acts::Geant4ShapeConverter{1, true}.rectangleBounds(sensorXz);
  auto boundsXz = std::get<0u>(rectangleXz);
  auto axesXz = std::get<1u>(rectangleXz);
  CHECK_CLOSE_ABS(boundsXz->halfLengthX(), 78., 10e-10);
  CHECK_CLOSE_ABS(boundsXz->halfLengthY(), 67., 10e-10);
  auto refXz = std::array<int, 2u>{0, -2};
  BOOST_CHECK(axesXz == refXz);
}

BOOST_AUTO_TEST_CASE(Geant4TrapzoidConversion) {
  // Standard TRD: XY are already well defined
  G4Trd trdXY("trdXY", 100, 150, 200, 200, 2);
  auto trapezoidXY = Acts::Geant4ShapeConverter{}.trapezoidBounds(trdXY);
  auto boundsXY = std::get<0u>(trapezoidXY);
  auto axesXY = std::get<1u>(trapezoidXY);
  CHECK_CLOSE_ABS(
      boundsXY->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 100,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsXY->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 150,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsXY->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 200,
      10e-10);
  auto refXY = std::array<int, 2u>{0, 1};
  BOOST_CHECK(axesXY == refXY);

  // Flipped, yX are the coordinates
  G4Trd trdyX("trdyX", 200, 200, 100, 150, 2);
  auto trapezoidyX = Acts::Geant4ShapeConverter{}.trapezoidBounds(trdyX);
  auto boundsyX = std::get<0u>(trapezoidyX);
  auto axesyX = std::get<1u>(trapezoidyX);
  CHECK_CLOSE_ABS(
      boundsyX->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 100,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsyX->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 150,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsyX->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 200,
      10e-10);
  auto refyX = std::array<int, 2u>{-1, 0};
  BOOST_CHECK(axesyX == refyX);

  // YZ span the trapezoid
  G4Trd trdYZ("trdYZ", 2, 2, 120, 140, 200);
  auto trapezoidYZ = Acts::Geant4ShapeConverter{}.trapezoidBounds(trdYZ);
  auto boundsYZ = std::get<0u>(trapezoidYZ);
  auto axesYZ = std::get<1u>(trapezoidYZ);
  CHECK_CLOSE_ABS(
      boundsYZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 120.,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsYZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 140.,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsYZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 200.,
      10e-10);
  auto refYZ = std::array<int, 2u>{1, 2};
  BOOST_CHECK(axesYZ == refYZ);

  // Xz span the trapezoid
  G4Trd trdXz("trdXz", 50, 75, 1, 1, 200);
  auto trapezoidXz = Acts::Geant4ShapeConverter{}.trapezoidBounds(trdXz);
  auto boundsXz = std::get<0u>(trapezoidXz);
  auto axesXz = std::get<1u>(trapezoidXz);
  CHECK_CLOSE_ABS(
      boundsXz->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 50.,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsXz->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 75.,
      10e-10);
  CHECK_CLOSE_ABS(
      boundsXz->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 200.,
      10e-10);
  auto refXz = std::array<int, 2u>{0, -2};
  BOOST_CHECK(axesXz == refXz);
}

BOOST_AUTO_TEST_CASE(Geant4PlanarConversion) {
  G4Box boxXY("boxXY", 23., 34., 1.);
  auto pBoundsBox = std::get<0u>(Acts::Geant4ShapeConverter{}.planarBounds(boxXY));
  auto rBounds = dynamic_cast<const Acts::RectangleBounds*>(pBoundsBox.get());
  BOOST_CHECK(rBounds != nullptr);

  G4Trd trdXY("trdXY", 100, 150, 200, 200, 2);
  auto pBoundsTrd = std::get<0u>(Acts::Geant4ShapeConverter{}.planarBounds(trdXY));
  auto tBounds = dynamic_cast<const Acts::TrapezoidBounds*>(pBoundsTrd.get());
  BOOST_CHECK(tBounds != nullptr);

}

BOOST_AUTO_TEST_SUITE_END()
