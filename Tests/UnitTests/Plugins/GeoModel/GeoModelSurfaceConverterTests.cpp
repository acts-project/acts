// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/GeoModel/GeoModelSurfaceConverter.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoTrap.h>
#include <GeoModelKernel/GeoTrd.h>

Acts::GeometryContext tContext;
Acts::RotationMatrix3 idRotation = Acts::RotationMatrix3::Identity();
Acts::Transform3 idTransform = Acts::Transform3::Identity();

BOOST_AUTO_TEST_SUITE(GeoModelPlugin)

// GeoBox conversion test case
BOOST_AUTO_TEST_CASE(GeoBoxConversion) {
  auto material = new GeoMaterial("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object) - XY
  auto boxXY = new GeoBox(100, 200, 2);
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto fphysXY = new GeoFullPhysVol(logXY);

  auto [elementXY, surfaceXY] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysXY);

  BOOST_CHECK(surfaceXY->type() == Acts::Surface::SurfaceType::Plane);

  // Check the bounds
  const Acts::RectangleBounds* rBoundsXY =
      dynamic_cast<const Acts::RectangleBounds*>(&(surfaceXY->bounds()));
  BOOST_CHECK(rBoundsXY != nullptr);
  CHECK_CLOSE_ABS(rBoundsXY->halfLengthX(), 100, 1e-6);
  CHECK_CLOSE_ABS(rBoundsXY->halfLengthY(), 200, 1e-6);

  // Check the transform -> should be identity transform
  const Acts::Transform3& transformXY = surfaceXY->transform(tContext);
  BOOST_CHECK(transformXY.isApprox(idTransform));

  // (BOX object) - YZ
  auto boxYZ = new GeoBox(2, 200, 300);
  auto logYZ = new GeoLogVol("LogVolumeYZ", boxYZ, material);
  auto fphysYZ = new GeoFullPhysVol(logYZ);

  auto [elementYZ, surfaceYZ] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysYZ);

  BOOST_CHECK(surfaceYZ->type() == Acts::Surface::SurfaceType::Plane);
  const Acts::RectangleBounds* rBoundsYZ =
      dynamic_cast<const Acts::RectangleBounds*>(&(surfaceYZ->bounds()));
  BOOST_CHECK(rBoundsYZ != nullptr);
  CHECK_CLOSE_ABS(rBoundsYZ->halfLengthX(), 200, 1e-6);
  CHECK_CLOSE_ABS(rBoundsYZ->halfLengthY(), 300, 1e-6);

  // Check the transform -> should be cyclic permutation of the identity
  const Acts::Transform3& transformYZ = surfaceYZ->transform(tContext);

  Acts::RotationMatrix3 rotationYZ = transformYZ.rotation();
  BOOST_CHECK(rotationYZ.col(0).isApprox(idRotation.col(1)));
  BOOST_CHECK(rotationYZ.col(1).isApprox(idRotation.col(2)));
  BOOST_CHECK(rotationYZ.col(2).isApprox(idRotation.col(0)));

  // (BOX object) - XZ
  auto boxXZ = new GeoBox(400, 2, 300);
  auto logXZ = new GeoLogVol("LogVolumeXZ", boxXZ, material);
  auto fphysXZ = new GeoFullPhysVol(logXZ);

  auto [elementXZ, surfaceXZ] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysXZ);

  BOOST_CHECK(surfaceXZ->type() == Acts::Surface::SurfaceType::Plane);

  // Check the bounds
  const Acts::RectangleBounds* rBoundsXZ =
      dynamic_cast<const Acts::RectangleBounds*>(&(surfaceXZ->bounds()));

  BOOST_CHECK(rBoundsXZ != nullptr);
  CHECK_CLOSE_ABS(rBoundsXZ->halfLengthX(), 300, 1e-6);
  CHECK_CLOSE_ABS(rBoundsXZ->halfLengthY(), 400, 1e-6);

  // Check the transform -> should be cyclic permutation of the identity
  const Acts::Transform3& transformXZ = surfaceXZ->transform(tContext);

  Acts::RotationMatrix3 rotationXZ = transformXZ.rotation();
  BOOST_CHECK(rotationXZ.col(0).isApprox(idRotation.col(2)));
  BOOST_CHECK(rotationXZ.col(1).isApprox(idRotation.col(0)));
  BOOST_CHECK(rotationXZ.col(2).isApprox(idRotation.col(1)));
}

// GeoTrap conversion test case
BOOST_AUTO_TEST_CASE(GeoTrapConversion) {
  auto material = new GeoMaterial("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (BOX object) - XY
  auto boxXY = new GeoBox(100, 200, 2);
  auto logXY = new GeoLogVol("LogVolumeXY", boxXY, material);
  auto fphysXY = new GeoFullPhysVol(logXY);

  auto [elementXY, surfaceXY] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysXY);
}

// GeoTrap conversion test case
BOOST_AUTO_TEST_CASE(GeoTrdConversion) {
  auto material = new GeoMaterial("Material", 1.0);
  // Let's create a GeoFullPhysVol object

  // (Trapezoid object) - YZ
  //  GeoTrd (double XHalfLength1, double XHalfLength2, double YHalfLength1,
  //  double YHalfLength2, double ZHalfLength);
  auto trapYZ = new GeoTrd(2, 2, 50, 80, 60);
  auto logYZ = new GeoLogVol("LogVolumeYZ", trapYZ, material);
  auto fphysYZ = new GeoFullPhysVol(logYZ);

  auto [elementYZ, surfaceYZ] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysYZ);

  // Check the bounds
  const Acts::TrapezoidBounds* tBoundsYZ =
      dynamic_cast<const Acts::TrapezoidBounds*>(&(surfaceYZ->bounds()));

  BOOST_CHECK(tBoundsYZ != nullptr);
  CHECK_CLOSE_ABS(
      tBoundsYZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 50,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsYZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 80,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsYZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 60,
      1e-6);

  // Check the transform -> should be cyclic permutation of the identity
  const Acts::Transform3& transformYZ = surfaceYZ->transform(tContext);

  Acts::RotationMatrix3 rotationYZ = transformYZ.rotation();
  BOOST_CHECK(rotationYZ.col(0).isApprox(idRotation.col(1)));
  BOOST_CHECK(rotationYZ.col(1).isApprox(idRotation.col(2)));
  BOOST_CHECK(rotationYZ.col(2).isApprox(idRotation.col(0)));

  // (Trapezoid object) - YZ swapped
  //  GeoTrd (double XHalfLength1, double XHalfLength2, double YHalfLength1,
  //  double YHalfLength2, double ZHalfLength);
  auto trapYZs = new GeoTrd(2, 2, 80, 50, 60);
  auto logYZs = new GeoLogVol("LogVolumeYZs", trapYZs, material);
  auto fphysYZs = new GeoFullPhysVol(logYZs);

  auto [elementYZs, surfaceYZs] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysYZs);

  // Check the bounds
  const Acts::TrapezoidBounds* tBoundsYZs =
      dynamic_cast<const Acts::TrapezoidBounds*>(&(surfaceYZs->bounds()));

  BOOST_CHECK(tBoundsYZs != nullptr);
  CHECK_CLOSE_ABS(
      tBoundsYZs->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 50,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsYZs->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 80,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsYZs->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 60,
      1e-6);

  // Check the transform -> should be cyclic permutation of the identity
  const Acts::Transform3& transformYZs = surfaceYZs->transform(tContext);

  Acts::RotationMatrix3 rotationYZs = transformYZs.rotation();
  BOOST_CHECK(rotationYZs.col(0).isApprox(idRotation.col(1)));
  BOOST_CHECK(rotationYZs.col(1).isApprox(-idRotation.col(2)));
  BOOST_CHECK(rotationYZs.col(2).isApprox(-idRotation.col(0)));

  // (Trapezoid object) - XZ
  auto trapXZ = new GeoTrd(50, 80, 2, 2, 60);
  auto logXZ = new GeoLogVol("LogVolumeXZ", trapXZ, material);
  auto fphysXZ = new GeoFullPhysVol(logXZ);

  auto [elementXZ, surfaceXZ] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysXZ);

  // Check the bounds
  const Acts::TrapezoidBounds* tBoundsXZ =
      dynamic_cast<const Acts::TrapezoidBounds*>(&(surfaceXZ->bounds()));

  BOOST_CHECK(tBoundsXZ != nullptr);
  CHECK_CLOSE_ABS(
      tBoundsXZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 50,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsXZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 80,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsXZ->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 60,
      1e-6);

  // Check the transform -> cyclic permuttation not possible
  const Acts::Transform3& transformXZ = surfaceXZ->transform(tContext);

  Acts::RotationMatrix3 rotationXZ = transformXZ.rotation();
  BOOST_CHECK(rotationXZ.col(0).isApprox(idRotation.col(0)));
  BOOST_CHECK(rotationXZ.col(1).isApprox(idRotation.col(2)));
  BOOST_CHECK(rotationXZ.col(2).isApprox(-1 * idRotation.col(1)));

  // (Trapezoid object) - XZs (swapped)
  auto trapXZs = new GeoTrd(80, 50, 2, 2, 60);
  auto logXZs = new GeoLogVol("LogVolumeXZs", trapXZs, material);
  auto fphysXZs = new GeoFullPhysVol(logXZs);

  auto [elementXZs, surfaceXZs] =
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysXZs);

  // Check the bounds
  const Acts::TrapezoidBounds* tBoundsXZs =
      dynamic_cast<const Acts::TrapezoidBounds*>(&(surfaceXZs->bounds()));

  BOOST_CHECK(tBoundsXZs != nullptr);
  CHECK_CLOSE_ABS(
      tBoundsXZs->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY), 50,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsXZs->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY), 80,
      1e-6);
  CHECK_CLOSE_ABS(
      tBoundsXZs->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY), 60,
      1e-6);

  // Check the transform -> cyclic permuttation not possible
  const Acts::Transform3& transformXZs = surfaceXZs->transform(tContext);

  Acts::RotationMatrix3 rotationXZs = transformXZs.rotation();
  BOOST_CHECK(rotationXZs.col(0).isApprox(idRotation.col(0)));
  BOOST_CHECK(rotationXZs.col(1).isApprox(-idRotation.col(2)));
  BOOST_CHECK(rotationXZs.col(2).isApprox(idRotation.col(1)));

  // Double - trazoid -> throw exception
  auto trapDouble = new GeoTrd(50, 80, 50, 80, 60);
  auto logDouble = new GeoLogVol("LogVolumeDouble", trapDouble, material);
  auto fphysDouble = new GeoFullPhysVol(logDouble);

  BOOST_CHECK_THROW(
      Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(*fphysDouble),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
