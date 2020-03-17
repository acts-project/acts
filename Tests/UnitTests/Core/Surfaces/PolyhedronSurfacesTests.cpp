// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

// Helper
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/ObjTestWriter.hpp"

// The class to test
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/Polyhedron.hpp"

// Cone surface
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"

// Cylinder surface
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

// Disc Surface
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

// Plane Surface
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

// Straw Surface
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"

#include "Acts/Utilities/ObjHelper.hpp"
#include "Acts/Utilities/Units.hpp"

#include <fstream>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {

using namespace UnitLiterals;

using IdentifiedPolyderon = std::tuple<std::string, bool, Polyhedron>;

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

std::vector<std::tuple<std::string, bool, unsigned int>> testModes = {
    {"", false, 72}, {"Triangulate", true, 72}, {"Extremas", false, 1}};

auto transform = std::make_shared<Transform3D>(Transform3D::Identity());

BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit tests for Cone Surfaces
BOOST_AUTO_TEST_CASE(ConeSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  double hzpmin = 10_mm;
  double hzpos = 35_mm;
  double hzneg = -20_mm;
  double alpha = 0.234;
  double phiSector = 0.358;
  ObjTestWriter::writeSectorPlanesObj("ConeSectorPlanes", phiSector, 0., hzpos,
                                      hzpos);

  for (const auto& mode : testModes) {
    unsigned int segments = std::get<unsigned int>(mode);
    std::string modename = std::get<std::string>(mode);
    bool modetrg = std::get<bool>(mode);

    /// The full cone on one side
    auto cone = std::make_shared<ConeBounds>(alpha, 0_mm, hzpos);
    auto oneCone = Surface::makeShared<ConeSurface>(transform, cone);
    auto oneConePh = oneCone->polyhedronRepresentation(tgContext, segments);
    size_t expectedFaces = segments < 4 ? 4 : segments;
    BOOST_CHECK_EQUAL(oneConePh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(oneConePh.vertices.size(), expectedFaces + 1);
    // Check the extent in space
    double r = hzpos * std::tan(alpha);
    auto extent = oneConePh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hzpos, 1e-6);
    testTypes.push_back({"ConeOneFull" + modename, modetrg, oneConePh});

    /// The full cone on one side
    auto conePiece = std::make_shared<ConeBounds>(alpha, hzpmin, hzpos);
    auto oneConePiece = Surface::makeShared<ConeSurface>(transform, conePiece);
    auto oneConePiecePh =
        oneConePiece->polyhedronRepresentation(tgContext, segments);
    expectedFaces = segments < 4 ? 4 : segments;
    // Check the extent in space
    double rmin = hzpmin * std::tan(alpha);
    extent = oneConePiecePh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, rmin, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, hzpmin, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hzpos, 1e-6);
    testTypes.push_back(
        {"ConeOnePieceFull" + modename, modetrg, oneConePiecePh});

    // The full cone on both sides
    auto coneBoth = std::make_shared<ConeBounds>(alpha, hzneg, hzpos);
    auto twoCones = Surface::makeShared<ConeSurface>(transform, coneBoth);
    auto twoConesPh = twoCones->polyhedronRepresentation(tgContext, segments);
    expectedFaces = segments < 4 ? 8 : 2 * segments;
    BOOST_CHECK_EQUAL(twoConesPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(twoConesPh.vertices.size(), expectedFaces + 1);
    extent = twoConesPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, hzneg, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hzpos, 1e-6);
    testTypes.push_back({"ConesTwoFull" + modename, modetrg, twoConesPh});

    // A centered sectoral cone on both sides
    auto sectoralBoth =
        std::make_shared<ConeBounds>(alpha, hzneg, hzpos, phiSector, 0.);
    auto sectoralCones =
        Surface::makeShared<ConeSurface>(transform, sectoralBoth);
    auto sectoralConesPh =
        sectoralCones->polyhedronRepresentation(tgContext, segments);
    extent = sectoralConesPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0_mm, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, hzneg, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hzpos, 1e-6);
    testTypes.push_back({"ConesSectoral" + modename, modetrg, sectoralConesPh});
  }
  ObjTestWriter::writeObj(testTypes);
}

/// Unit tests for Cylinder Surfaces
BOOST_AUTO_TEST_CASE(CylinderSurfacePolyhedrons) {
  double r = 25_mm;
  double hZ = 35_mm;

  double phiSector = 0.458;
  double averagePhi = -1.345;
  ObjTestWriter::writeSectorPlanesObj("CylinderCentralSectorPlanes", phiSector,
                                      0., 1.5 * r, 1.5 * hZ);
  ObjTestWriter::writeSectorPlanesObj("CylinderShiftedSectorPlanes", phiSector,
                                      averagePhi, 1.5 * r, 1.5 * hZ);

  std::vector<IdentifiedPolyderon> testTypes;

  for (const auto& mode : testModes) {
    unsigned int segments = std::get<unsigned int>(mode);
    std::string modename = std::get<std::string>(mode);
    bool modetrg = std::get<bool>(mode);

    size_t expectedFaces = segments < 4 ? 4 : segments;
    size_t expectedVertices = segments < 4 ? 8 : 2 * segments;

    /// The full cone on one side
    auto cylinder = std::make_shared<CylinderBounds>(r, hZ);
    auto fullCylinder =
        Surface::makeShared<CylinderSurface>(transform, cylinder);
    auto fullCylinderPh =
        fullCylinder->polyhedronRepresentation(tgContext, segments);

    BOOST_CHECK_EQUAL(fullCylinderPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(fullCylinderPh.vertices.size(), expectedVertices);
    // Check the extent in space
    auto extent = fullCylinderPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, -hZ, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hZ, 1e-6);
    testTypes.push_back({"CylinderFull" + modename, modetrg, fullCylinderPh});

    /// The full cone on one side
    auto sectorCentered = std::make_shared<CylinderBounds>(r, phiSector, hZ);
    auto centerSectoredCylinder =
        Surface::makeShared<CylinderSurface>(transform, sectorCentered);
    auto centerSectoredCylinderPh =
        centerSectoredCylinder->polyhedronRepresentation(tgContext, segments);

    // Check the extent in space
    extent = centerSectoredCylinderPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, r * std::cos(phiSector), 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -r * std::sin(phiSector), 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, r * std::sin(phiSector), 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, -hZ, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hZ, 1e-6);
    testTypes.push_back({"CylinderSectorCentered" + modename, modetrg,
                         centerSectoredCylinderPh});

    /// The full cone on one side
    auto sectorShifted =
        std::make_shared<CylinderBounds>(r, averagePhi, phiSector, hZ);
    auto shiftedSectoredCylinder =
        Surface::makeShared<CylinderSurface>(transform, sectorShifted);
    auto shiftedSectoredCylinderPh =
        shiftedSectoredCylinder->polyhedronRepresentation(tgContext, segments);

    // Check the extent in space
    extent = shiftedSectoredCylinderPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binR].first, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, r, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, -hZ, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, hZ, 1e-6);
    testTypes.push_back({"CylinderSectorShifted" + modename, modetrg,
                         shiftedSectoredCylinderPh});
  }

  ObjTestWriter::writeObj(testTypes);
}

/// Unit tests for Disc Surfaces
BOOST_AUTO_TEST_CASE(DiscSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  double innerR = 10_mm;
  double outerR = 25_mm;

  double phiSector = 0.345;
  double averagePhi = -1.0;

  double cphi = std::cos(phiSector);
  double sphi = std::sin(phiSector);

  std::pair<Vector3D, Vector3D> lineA = {
      Vector3D(0., 0., 0.), Vector3D(outerR * cphi, outerR * sphi, 0.)};
  std::pair<Vector3D, Vector3D> lineB = {
      Vector3D(0., 0., 0.), Vector3D(outerR * cphi, -outerR * sphi, 0.)};
  ObjTestWriter::writeSectorLinesObj("DiscSectorLines", lineA, lineB);

  double minPhi = averagePhi - phiSector;
  double maxPhi = averagePhi + phiSector;
  lineA = {Vector3D(0., 0., 0.),
           Vector3D(outerR * std::cos(minPhi), outerR * std::sin(minPhi), 0.)};
  lineB = {Vector3D(0., 0., 0.),
           Vector3D(outerR * std::cos(maxPhi), outerR * std::sin(maxPhi), 0.)};
  ObjTestWriter::writeSectorLinesObj("DiscSectorLinesShifted", lineA, lineB);

  for (const auto& mode : testModes) {
    unsigned int segments = std::get<unsigned int>(mode);
    std::string modename = std::get<std::string>(mode);
    bool modetrg = std::get<bool>(mode);

    // Full disc
    auto disc = std::make_shared<RadialBounds>(0_mm, outerR);
    auto fullDisc = Surface::makeShared<DiscSurface>(transform, disc);
    auto fullDiscPh = fullDisc->polyhedronRepresentation(tgContext, segments);

    unsigned int expectedVertices = segments > 4 ? segments : 4;
    unsigned int expectedFaces = 1;

    BOOST_CHECK_EQUAL(fullDiscPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(fullDiscPh.vertices.size(), expectedVertices);

    auto extent = fullDiscPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);

    testTypes.push_back({"DiscFull" + modename, modetrg, fullDiscPh});

    // Ring disc
    auto radial = std::make_shared<RadialBounds>(innerR, outerR);
    auto radialDisc = Surface::makeShared<DiscSurface>(transform, radial);
    auto radialPh = radialDisc->polyhedronRepresentation(tgContext, segments);
    extent = radialPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, innerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back({"DiscRing" + modename, modetrg, radialPh});

    // Sectoral disc - around 0.
    auto sector = std::make_shared<RadialBounds>(0., outerR, phiSector);
    auto sectorDisc = Surface::makeShared<DiscSurface>(transform, sector);
    auto sectorPh = sectorDisc->polyhedronRepresentation(tgContext, segments);
    extent = sectorPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back({"DiscSectorCentered" + modename, modetrg, sectorPh});

    // Sectoral ring - around 0.
    auto sectorRing = std::make_shared<RadialBounds>(innerR, outerR, phiSector);
    auto sectorRingDisc =
        Surface::makeShared<DiscSurface>(transform, sectorRing);
    auto sectorRingDiscPh =
        sectorRingDisc->polyhedronRepresentation(tgContext, segments);
    extent = sectorRingDiscPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, innerR * std::cos(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, innerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back(
        {"DiscRingSectorCentered" + modename, modetrg, sectorRingDiscPh});

    // Sectoral disc - shifted
    auto sectorRingShifted =
        std::make_shared<RadialBounds>(innerR, outerR, averagePhi, phiSector);
    auto sectorRingDiscShifted =
        Surface::makeShared<DiscSurface>(transform, sectorRingShifted);
    auto sectorRingDiscShiftedPh =
        sectorRingDiscShifted->polyhedronRepresentation(tgContext, segments);
    extent = sectorRingDiscShiftedPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binR].first, innerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back(
        {"DiscRingSectorShifted" + modename, modetrg, sectorRingDiscShiftedPh});

    // Trapezoid for a disc
    double halfXmin = 10_mm;
    double halfXmax = 20_mm;
    auto trapezoidDisc = std::make_shared<DiscTrapezoidBounds>(
        halfXmin, halfXmax, innerR, outerR, 0.);
    auto trapezoidDiscSf =
        Surface::makeShared<DiscSurface>(transform, trapezoidDisc);
    auto trapezoidDiscSfPh =
        trapezoidDiscSf->polyhedronRepresentation(tgContext, segments);
    extent = trapezoidDiscSfPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binR].first, innerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back(
        {"DiscTrapezoidCentered" + modename, modetrg, trapezoidDiscSfPh});

    auto trapezoidDiscShifted = std::make_shared<DiscTrapezoidBounds>(
        halfXmin, halfXmax, innerR, outerR, averagePhi);
    auto trapezoidDiscShiftedSf =
        Surface::makeShared<DiscSurface>(transform, trapezoidDiscShifted);
    auto trapezoidDiscShiftedSfPh =
        trapezoidDiscShiftedSf->polyhedronRepresentation(tgContext, segments);
    extent = trapezoidDiscShiftedSfPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binR].first, innerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, outerR, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back(
        {"DiscTrapezoidShifted" + modename, modetrg, trapezoidDiscShiftedSfPh});

    double minRadius = 7.;
    double maxRadius = 12.;
    double minPhiA = 0.75;
    double maxPhiA = 1.4;

    Vector2D offset(-2., 2.);

    auto annulus = std::make_shared<AnnulusBounds>(minRadius, maxRadius,
                                                   minPhiA, maxPhiA, offset);
    auto annulusDisc = Surface::makeShared<DiscSurface>(transform, annulus);
    auto annulusDiscPh =
        annulusDisc->polyhedronRepresentation(tgContext, segments);

    testTypes.push_back(
        {"DiscAnnulus" + modename, modetrg, trapezoidDiscShiftedSfPh});
  }

  ObjTestWriter::writeObj(testTypes);
}

/// Unit tests for Plane Surfaces
BOOST_AUTO_TEST_CASE(PlaneSurfacePolyhedrons) {
  std::vector<IdentifiedPolyderon> testTypes;

  double rhX = 10_mm;
  double rhY = 25_mm;
  double shiftY = 50_mm;
  auto rectangular = std::make_shared<RectangleBounds>(rhX, rhY);

  // Special test for shifted plane to check rMin/rMax
  Vector3D shift(0., shiftY, 0.);
  auto shiftedTransform =
      std::make_shared<Transform3D>(Transform3D::Identity());
  shiftedTransform->pretranslate(shift);
  auto shiftedPlane =
      Surface::makeShared<PlaneSurface>(shiftedTransform, rectangular);
  auto shiftedPh = shiftedPlane->polyhedronRepresentation(tgContext, 1);
  auto shiftedExtent = shiftedPh.extent();
  // Let's check the extent
  CHECK_CLOSE_ABS(shiftedExtent.ranges[binX].first, -rhX, 1e-6);
  CHECK_CLOSE_ABS(shiftedExtent.ranges[binX].second, rhX, 1e-6);
  CHECK_CLOSE_ABS(shiftedExtent.ranges[binY].first, -rhY + shiftY, 1e-6);
  CHECK_CLOSE_ABS(shiftedExtent.ranges[binY].second, rhY + shiftY, 1e-6);

  for (const auto& mode : testModes) {
    unsigned int segments = std::get<unsigned int>(mode);
    std::string modename = std::get<std::string>(mode);
    bool modetrg = std::get<bool>(mode);

    /// Rectangular Plane
    auto rectangularPlane =
        Surface::makeShared<PlaneSurface>(transform, rectangular);
    auto rectangularPh =
        rectangularPlane->polyhedronRepresentation(tgContext, segments);
    auto extent = rectangularPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -rhX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, rhX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -rhY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, rhY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second,
                    std::sqrt(rhX * rhX + rhY * rhY), 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    BOOST_CHECK(rectangularPh.vertices.size() == 4);
    BOOST_CHECK(rectangularPh.faces.size() == 1);
    std::vector<size_t> expectedRect = {0, 1, 2, 3};
    BOOST_CHECK(rectangularPh.faces[0] == expectedRect);
    testTypes.push_back({"PlaneRectangle" + modename, modetrg, rectangularPh});

    /// Trapezoidal Plane
    double thX1 = 10_mm;
    double thX2 = 25_mm;
    double thY = 35_mm;

    auto trapezoid = std::make_shared<TrapezoidBounds>(thX1, thX2, thY);
    auto trapezoidalPlane =
        Surface::makeShared<PlaneSurface>(transform, trapezoid);
    auto trapezoidalPh =
        trapezoidalPlane->polyhedronRepresentation(tgContext, segments);
    extent = trapezoidalPh.extent();

    double thX = std::max(thX1, thX2);
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -thX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, thX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -thY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, thY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second,
                    std::sqrt(thX * thX + thY * thY), 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    BOOST_CHECK(trapezoidalPh.vertices.size() == 4);
    BOOST_CHECK(trapezoidalPh.faces.size() == 1);
    std::vector<size_t> expectedTra = {0, 1, 2, 3};
    BOOST_CHECK(trapezoidalPh.faces[0] == expectedTra);
    testTypes.push_back({"PlaneTrapezoid" + modename, modetrg, trapezoidalPh});

    /// Ring-like ellispoidal Plane
    double rMaxX = 30_mm;
    double rMaxY = 40_mm;
    auto ellipse = std::make_shared<EllipseBounds>(0., 0., rMaxX, rMaxY);
    auto ellipsoidPlane = Surface::makeShared<PlaneSurface>(transform, ellipse);
    auto ellispoidPh =
        ellipsoidPlane->polyhedronRepresentation(tgContext, segments);
    extent = ellispoidPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -rMaxX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, rMaxX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -rMaxY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -rMaxY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, rMaxY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);

    testTypes.push_back({"PlaneFullEllipse" + modename, modetrg, ellispoidPh});

    double rMinX = 10_mm;
    double rMinY = 20_mm;
    auto ellipseRing =
        std::make_shared<EllipseBounds>(rMinX, rMaxX, rMinY, rMaxY);
    auto ellipsoidRingPlane =
        Surface::makeShared<PlaneSurface>(transform, ellipseRing);
    auto ellispoidRingPh =
        ellipsoidRingPlane->polyhedronRepresentation(tgContext, segments);

    extent = ellispoidPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -rMaxX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, rMaxX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -rMaxY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, rMinX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second, rMaxY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);

    // BOOST_CHECK(ellispoidPh.vertices.size()==72);
    testTypes.push_back(
        {"PlaneRingEllipse" + modename, modetrg, ellispoidRingPh});

    /// ConvextPolygonBounds test
    std::vector<Vector2D> vtxs = {
        Vector2D(-40_mm, -10_mm), Vector2D(-10_mm, -30_mm),
        Vector2D(30_mm, -20_mm),  Vector2D(10_mm, 20_mm),
        Vector2D(-20_mm, 50_mm),  Vector2D(-30_mm, 30_mm)};

    auto sextagon = std::make_shared<ConvexPolygonBounds<6>>(vtxs);
    auto sextagonPlane = Surface::makeShared<PlaneSurface>(transform, sextagon);
    auto sextagonPlanePh =
        sextagonPlane->polyhedronRepresentation(tgContext, segments);
    testTypes.push_back({"PlaneSextagon" + modename, modetrg, sextagonPlanePh});

    /// Diamond shaped plane
    double hMinX = 10_mm;
    double hMedX = 20_mm;
    double hMaxX = 15_mm;
    double hMinY = 40_mm;
    double hMaxY = 50_mm;
    auto diamond =
        std::make_shared<DiamondBounds>(hMinX, hMedX, hMaxX, hMinY, hMaxY);
    auto diamondPlane = Surface::makeShared<PlaneSurface>(transform, diamond);
    auto diamondPh =
        diamondPlane->polyhedronRepresentation(tgContext, segments);
    BOOST_CHECK(diamondPh.vertices.size() == 6);
    BOOST_CHECK(diamondPh.faces.size() == 1);
    extent = diamondPh.extent();
    CHECK_CLOSE_ABS(extent.ranges[binX].first, -hMedX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binX].second, hMedX, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].first, -hMinY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binY].second, hMaxY, 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binR].second,
                    std::sqrt(hMaxX * hMaxX + hMaxY * hMaxY), 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].first, 0., 1e-6);
    CHECK_CLOSE_ABS(extent.ranges[binZ].second, 0., 1e-6);
    testTypes.push_back({"PlaneDiamond" + modename, modetrg, diamondPh});
  }
  ObjTestWriter::writeObj(testTypes);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts