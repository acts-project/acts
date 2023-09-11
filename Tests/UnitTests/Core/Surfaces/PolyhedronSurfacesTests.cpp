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
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <fstream>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {

using namespace UnitLiterals;

using IdentifiedPolyhedron = std::tuple<std::string, bool, Polyhedron>;

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

std::vector<std::tuple<std::string, bool, unsigned int>> testModes = {
    {"", false, 72}, {"Triangulate", true, 72}, {"Extremas", false, 1}};

auto transform = std::make_shared<Transform3>(Transform3::Identity());

BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit tests for Cone Surfaces
BOOST_AUTO_TEST_CASE(ConeSurfacePolyhedrons) {
  double hzpmin = 10_mm;
  double hzpos = 35_mm;
  double hzneg = -20_mm;
  double alpha = 0.234;
  double phiSector = 0.358;

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
    CHECK_CLOSE_ABS((extent.range(binX).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0_mm, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0_mm, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hzpos, 1e-6);

    /// The full cone on one side
    auto conePiece = std::make_shared<ConeBounds>(alpha, hzpmin, hzpos);
    auto oneConePiece = Surface::makeShared<ConeSurface>(transform, conePiece);
    auto oneConePiecePh =
        oneConePiece->polyhedronRepresentation(tgContext, segments);
    expectedFaces = segments < 4 ? 4 : segments;
    // Check the extent in space
    double rmin = hzpmin * std::tan(alpha);
    extent = oneConePiecePh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() rmin, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() hzpmin, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hzpos, 1e-6);

    // The full cone on both sides
    auto coneBoth = std::make_shared<ConeBounds>(alpha, hzneg, hzpos);
    auto twoCones = Surface::makeShared<ConeSurface>(transform, coneBoth);
    auto twoConesPh = twoCones->polyhedronRepresentation(tgContext, segments);
    expectedFaces = segments < 4 ? 8 : 2 * segments;
    BOOST_CHECK_EQUAL(twoConesPh.faces.size(), expectedFaces);
    BOOST_CHECK_EQUAL(twoConesPh.vertices.size(), expectedFaces + 1);
    extent = twoConesPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0_mm, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() hzneg, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hzpos, 1e-6);

    // A centered sectoral cone on both sides
    auto sectoralBoth =
        std::make_shared<ConeBounds>(alpha, hzneg, hzpos, phiSector, 0.);
    auto sectoralCones =
        Surface::makeShared<ConeSurface>(transform, sectoralBoth);
    auto sectoralConesPh =
        sectoralCones->polyhedronRepresentation(tgContext, segments);
    extent = sectoralConesPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0_mm, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() hzneg, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hzpos, 1e-6);
  }
}

/// Unit tests for Cylinder Surfaces
BOOST_AUTO_TEST_CASE(CylinderSurfacePolyhedrons) {
  double r = 25_mm;
  double hZ = 35_mm;

  double phiSector = 0.458;
  double averagePhi = -1.345;

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
    CHECK_CLOSE_ABS((extent.range(binX).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() -hZ, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hZ, 1e-6);

    /// The full cone on one side
    auto sectorCentered = std::make_shared<CylinderBounds>(r, phiSector, hZ);
    auto centerSectoredCylinder =
        Surface::makeShared<CylinderSurface>(transform, sectorCentered);
    auto centerSectoredCylinderPh =
        centerSectoredCylinder->polyhedronRepresentation(tgContext, segments);

    // Check the extent in space
    extent = centerSectoredCylinderPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() r * std::cos(phiSector), 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -r * std::sin(phiSector), 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() r * std::sin(phiSector), 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() -hZ, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hZ, 1e-6);

    /// The full cone on one side
    auto sectorShifted =
        std::make_shared<CylinderBounds>(r, averagePhi, phiSector, hZ);
    auto shiftedSectoredCylinder =
        Surface::makeShared<CylinderSurface>(transform, sectorShifted);
    auto shiftedSectoredCylinderPh =
        shiftedSectoredCylinder->polyhedronRepresentation(tgContext, segments);

    // Check the extent in space
    extent = shiftedSectoredCylinderPh.extent();
    CHECK_CLOSE_ABS((extent.range(binR).min() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() r, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() -hZ, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() hZ, 1e-6);
  }
}

/// Unit tests for Disc Surfaces
BOOST_AUTO_TEST_CASE(DiscSurfacePolyhedrons) {
  double innerR = 10_mm;
  double outerR = 25_mm;

  double phiSector = 0.345;
  double averagePhi = -1.0;

  double cphi = std::cos(phiSector);
  double sphi = std::sin(phiSector);

  std::pair<Vector3, Vector3> lineA = {
      Vector3(0., 0., 0.), Vector3(outerR * cphi, outerR * sphi, 0.)};
  std::pair<Vector3, Vector3> lineB = {
      Vector3(0., 0., 0.), Vector3(outerR * cphi, -outerR * sphi, 0.)};

  double minPhi = averagePhi - phiSector;
  double maxPhi = averagePhi + phiSector;
  lineA = {Vector3(0., 0., 0.),
           Vector3(outerR * std::cos(minPhi), outerR * std::sin(minPhi), 0.)};
  lineB = {Vector3(0., 0., 0.),
           Vector3(outerR * std::cos(maxPhi), outerR * std::sin(maxPhi), 0.)};

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
    CHECK_CLOSE_ABS((extent.range(binX).min() -outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    // Ring disc
    auto radial = std::make_shared<RadialBounds>(innerR, outerR);
    auto radialDisc = Surface::makeShared<DiscSurface>(transform, radial);
    auto radialPh = radialDisc->polyhedronRepresentation(tgContext, segments);
    extent = radialPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() -outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() innerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    // Sectoral disc - around 0.
    auto sector = std::make_shared<RadialBounds>(0., outerR, phiSector);
    auto sectorDisc = Surface::makeShared<DiscSurface>(transform, sector);
    auto sectorPh = sectorDisc->polyhedronRepresentation(tgContext, segments);
    extent = sectorPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    // Sectoral ring - around 0.
    auto sectorRing = std::make_shared<RadialBounds>(innerR, outerR, phiSector);
    auto sectorRingDisc =
        Surface::makeShared<DiscSurface>(transform, sectorRing);
    auto sectorRingDiscPh =
        sectorRingDisc->polyhedronRepresentation(tgContext, segments);
    extent = sectorRingDiscPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() innerR * std::cos(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() outerR * std::sin(phiSector),
                    1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() innerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    // Sectoral disc - shifted
    auto sectorRingShifted =
        std::make_shared<RadialBounds>(innerR, outerR, averagePhi, phiSector);
    auto sectorRingDiscShifted =
        Surface::makeShared<DiscSurface>(transform, sectorRingShifted);
    auto sectorRingDiscShiftedPh =
        sectorRingDiscShifted->polyhedronRepresentation(tgContext, segments);
    extent = sectorRingDiscShiftedPh.extent();
    CHECK_CLOSE_ABS((extent.range(binR).min() innerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

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
    CHECK_CLOSE_ABS((extent.range(binR).min() innerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    auto trapezoidDiscShifted = std::make_shared<DiscTrapezoidBounds>(
        halfXmin, halfXmax, innerR, outerR, averagePhi);
    auto trapezoidDiscShiftedSf =
        Surface::makeShared<DiscSurface>(transform, trapezoidDiscShifted);
    auto trapezoidDiscShiftedSfPh =
        trapezoidDiscShiftedSf->polyhedronRepresentation(tgContext, segments);
    extent = trapezoidDiscShiftedSfPh.extent();
    CHECK_CLOSE_ABS((extent.range(binR).min() innerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() outerR, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    double minRadius = 7.;
    double maxRadius = 12.;
    double minPhiA = 0.75;
    double maxPhiA = 1.4;

    Vector2 offset(-2., 2.);

    auto annulus = std::make_shared<AnnulusBounds>(minRadius, maxRadius,
                                                   minPhiA, maxPhiA, offset);
    auto annulusDisc = Surface::makeShared<DiscSurface>(transform, annulus);
    auto annulusDiscPh =
        annulusDisc->polyhedronRepresentation(tgContext, segments);
  }
}

/// Unit tests for Plane Surfaces
BOOST_AUTO_TEST_CASE(PlaneSurfacePolyhedrons) {
  double rhX = 10_mm;
  double rhY = 25_mm;
  double shiftY = 50_mm;
  auto rectangular = std::make_shared<RectangleBounds>(rhX, rhY);

  // Special test for shifted plane to check rMin/rMax
  Vector3 shift(0., shiftY, 0.);
  auto shiftedTransform = std::make_shared<Transform3>(Transform3::Identity());
  shiftedTransform->pretranslate(shift);
  auto shiftedPlane =
      Surface::makeShared<PlaneSurface>(shiftedTransform, rectangular);
  auto shiftedPh = shiftedPlane->polyhedronRepresentation(tgContext, 1);
  auto shiftedExtent = shiftedPh.extent();
  // Let's check the extent
  CHECK_CLOSE_ABS(shifted(extent.range(binX).min() -rhX, 1e-6);
  CHECK_CLOSE_ABS(shifted(extent.range(binX).max() rhX, 1e-6);
  CHECK_CLOSE_ABS(shifted(extent.range(binY).min() -rhY + shiftY, 1e-6);
  CHECK_CLOSE_ABS(shifted(extent.range(binY).max() rhY + shiftY, 1e-6);

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
    CHECK_CLOSE_ABS((extent.range(binX).min() -rhX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() rhX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -rhY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() rhY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max()
                    std::hypot(rhX, rhY), 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);
    BOOST_CHECK(rectangularPh.vertices.size() == 4);
    BOOST_CHECK(rectangularPh.faces.size() == 1);
    std::vector<size_t> expectedRect = {0, 1, 2, 3};
    BOOST_CHECK(rectangularPh.faces[0] == expectedRect);

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
    CHECK_CLOSE_ABS((extent.range(binX).min() -thX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() thX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -thY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() thY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max()
                    std::hypot(thX, thY), 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);
    BOOST_CHECK(trapezoidalPh.vertices.size() == 4);
    BOOST_CHECK(trapezoidalPh.faces.size() == 1);
    std::vector<size_t> expectedTra = {0, 1, 2, 3};
    BOOST_CHECK(trapezoidalPh.faces[0] == expectedTra);

    /// Ring-like ellispoidal Plane
    double rMaxX = 30_mm;
    double rMaxY = 40_mm;
    auto ellipse = std::make_shared<EllipseBounds>(0., 0., rMaxX, rMaxY);
    auto ellipsoidPlane = Surface::makeShared<PlaneSurface>(transform, ellipse);
    auto ellispoidPh =
        ellipsoidPlane->polyhedronRepresentation(tgContext, segments);
    extent = ellispoidPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() -rMaxX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() rMaxX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -rMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() rMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() rMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    double rMinX = 10_mm;
    double rMinY = 20_mm;
    auto ellipseRing =
        std::make_shared<EllipseBounds>(rMinX, rMaxX, rMinY, rMaxY);
    auto ellipsoidRingPlane =
        Surface::makeShared<PlaneSurface>(transform, ellipseRing);
    auto ellispoidRingPh =
        ellipsoidRingPlane->polyhedronRepresentation(tgContext, segments);

    extent = ellispoidPh.extent();
    CHECK_CLOSE_ABS((extent.range(binX).min() -rMaxX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() rMaxX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -rMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() rMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() rMinX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max() rMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);

    /// ConvextPolygonBounds test
    std::vector<Vector2> vtxs = {
        Vector2(-40_mm, -10_mm), Vector2(-10_mm, -30_mm),
        Vector2(30_mm, -20_mm),  Vector2(10_mm, 20_mm),
        Vector2(-20_mm, 50_mm),  Vector2(-30_mm, 30_mm)};

    auto sextagon = std::make_shared<ConvexPolygonBounds<6>>(vtxs);
    auto sextagonPlane = Surface::makeShared<PlaneSurface>(transform, sextagon);
    auto sextagonPlanePh =
        sextagonPlane->polyhedronRepresentation(tgContext, segments);

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
    CHECK_CLOSE_ABS((extent.range(binX).min() -hMedX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binX).max() hMedX, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).min() -hMinY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binY).max() hMaxY, 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binR).max()
                    std::hypot(hMaxX, hMaxY), 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).min() 0., 1e-6);
    CHECK_CLOSE_ABS((extent.range(binZ).max() 0., 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
