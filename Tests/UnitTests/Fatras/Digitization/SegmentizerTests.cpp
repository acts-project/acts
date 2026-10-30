// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "ActsFatras/Digitization/Segmentizer.hpp"

#include <cmath>
#include <fstream>
#include <memory>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "DigitizationCsvOutput.hpp"
#include "PlanarSurfaceTestBeds.hpp"

namespace bdata = boost::unit_test::data;

using namespace Acts;
using namespace ActsFatras;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DigitizationSuite)

BOOST_AUTO_TEST_CASE(SegmentizerCartesian) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();

  auto rectangleBounds = std::make_shared<RectangleBounds>(1., 1.);
  auto planeSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                        rectangleBounds);

  // The segmentation
  const auto segmentation = IMultiAxis::create(
      *IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 20,
                                AxisDirection::AxisX),
      *IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 20,
                                AxisDirection::AxisY));

  Segmentizer cl;

  // Test: Normal hit into the surface
  Vector2 nPosition(0.37, 0.76);
  auto nSegments =
      cl.segments(geoCtx, *planeSurface, *segmentation, {nPosition, nPosition});
  BOOST_CHECK_EQUAL(nSegments.size(), 1);
  BOOST_CHECK_EQUAL(nSegments[0].bin[0], 13);
  BOOST_CHECK_EQUAL(nSegments[0].bin[1], 17);

  // Test: Inclined hit into the surface - negative x direction
  Vector2 ixPositionS(0.37, 0.76);
  Vector2 ixPositionE(0.02, 0.73);
  auto ixSegments = cl.segments(geoCtx, *planeSurface, *segmentation,
                                {ixPositionS, ixPositionE});
  BOOST_CHECK_EQUAL(ixSegments.size(), 4);

  // Test: Inclined hit into the surface - positive y direction
  Vector2 iyPositionS(0.37, 0.76);
  Vector2 iyPositionE(0.39, 0.91);
  auto iySegments = cl.segments(geoCtx, *planeSurface, *segmentation,
                                {iyPositionS, iyPositionE});
  BOOST_CHECK_EQUAL(iySegments.size(), 3);

  // Test: Inclined hit into the surface - x/y direction
  Vector2 ixyPositionS(-0.27, 0.76);
  Vector2 ixyPositionE(-0.02, -0.73);
  auto ixySegments = cl.segments(geoCtx, *planeSurface, *segmentation,
                                 {ixyPositionS, ixyPositionE});
  BOOST_CHECK_EQUAL(ixySegments.size(), 18);
}

BOOST_AUTO_TEST_CASE(SegmentizerPolarRadial) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();

  auto radialBounds = std::make_shared<const RadialBounds>(5., 10., 0.25, 0.);
  auto radialDisc =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), radialBounds);

  // The segmentation
  const auto segmentation = IMultiAxis::create(
      *IAxis::createEquidistant(AxisBoundaryType::Bound, 5., 10., 2,
                                AxisDirection::AxisR),
      *IAxis::createEquidistant(AxisBoundaryType::Bound, -0.25, 0.25, 250,
                                AxisDirection::AxisPhi));

  Segmentizer cl;

  // Test: Normal hit into the surface
  Vector2 nPosition(6.76, 0.5);
  auto nSegments =
      cl.segments(geoCtx, *radialDisc, *segmentation, {nPosition, nPosition});
  BOOST_CHECK_EQUAL(nSegments.size(), 1);
  BOOST_CHECK_EQUAL(nSegments[0].bin[0], 0);
  BOOST_CHECK_EQUAL(nSegments[0].bin[1], 161);

  // Test: now opver more phi strips
  Vector2 sPositionS(6.76, 0.5);
  Vector2 sPositionE(7.03, -0.3);
  auto sSegment =
      cl.segments(geoCtx, *radialDisc, *segmentation, {sPositionS, sPositionE});
  BOOST_CHECK_EQUAL(sSegment.size(), 59);

  // Test: jump over R boundary, but stay in phi bin
  sPositionS = Vector2(6.76, 0.);
  sPositionE = Vector2(7.83, 0.);
  sSegment =
      cl.segments(geoCtx, *radialDisc, *segmentation, {sPositionS, sPositionE});
  BOOST_CHECK_EQUAL(sSegment.size(), 2);
}

BOOST_AUTO_TEST_CASE(SegmentizerCylinderPhiSeam) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();

  const double radius = 5.;
  const double halfZ = 10.;
  const double halfRPhi = std::numbers::pi * radius;
  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(radius, halfZ));

  // Unrolled (rPhi, z) readout frame, closed in rPhi
  const std::size_t nRPhi = 100;
  const auto segmentation = IMultiAxis::create(
      *IAxis::createEquidistant(AxisBoundaryType::Closed, -halfRPhi, halfRPhi,
                                nRPhi, AxisDirection::AxisRPhi),
      *IAxis::createEquidistant(AxisBoundaryType::Bound, -halfZ, halfZ, 20,
                                AxisDirection::AxisZ));
  const auto openSegmentation = IMultiAxis::create(
      *IAxis::createEquidistant(AxisBoundaryType::Bound, -halfRPhi, halfRPhi,
                                nRPhi, AxisDirection::AxisRPhi),
      *IAxis::createEquidistant(AxisBoundaryType::Bound, -halfZ, halfZ, 20,
                                AxisDirection::AxisZ));

  Segmentizer cl;

  auto totalPath = [](const auto& segments) {
    double path = 0.;
    for (const auto& s : segments) {
      path += s.activation;
    }
    return path;
  };

  // Test: away from the seam, the closed axis behaves like an open one
  Vector2 aStart(0.1, 0.3);
  Vector2 aEnd(1.4, 0.9);
  auto aSegments =
      cl.segments(geoCtx, *cylinder, *segmentation, {aStart, aEnd});
  auto aOpenSegments =
      cl.segments(geoCtx, *cylinder, *openSegmentation, {aStart, aEnd});
  BOOST_REQUIRE_EQUAL(aSegments.size(), aOpenSegments.size());
  for (std::size_t i = 0; i < aSegments.size(); ++i) {
    BOOST_CHECK_EQUAL(aSegments[i].bin[0], aOpenSegments[i].bin[0]);
    BOOST_CHECK_EQUAL(aSegments[i].bin[1], aOpenSegments[i].bin[1]);
    BOOST_CHECK_EQUAL(aSegments[i].activation, aOpenSegments[i].activation);
  }

  // Test: crossing the seam at +pi R in rPhi gives two channels, not the ring
  Vector2 sStart(halfRPhi - 0.05, 0.3);
  Vector2 sEnd(halfRPhi + 0.05, 0.3);
  auto sSegments =
      cl.segments(geoCtx, *cylinder, *segmentation, {sStart, sEnd});
  BOOST_REQUIRE_EQUAL(sSegments.size(), 2);
  BOOST_CHECK_EQUAL(sSegments[0].bin[0], nRPhi - 1);
  BOOST_CHECK_EQUAL(sSegments[1].bin[0], 0);
  BOOST_CHECK_CLOSE(totalPath(sSegments), 0.1, 1e-6);

  // Test: the same in the other direction
  auto rSegments =
      cl.segments(geoCtx, *cylinder, *segmentation, {sEnd, sStart});
  BOOST_REQUIRE_EQUAL(rSegments.size(), 2);
  BOOST_CHECK_EQUAL(rSegments[0].bin[0], 0);
  BOOST_CHECK_EQUAL(rSegments[1].bin[0], nRPhi - 1);

  // Test: crossing the seam at -pi R in rPhi
  Vector2 nStart(-halfRPhi + 0.05, 0.3);
  Vector2 nEnd(-halfRPhi - 0.05, 0.3);
  auto nSegments =
      cl.segments(geoCtx, *cylinder, *segmentation, {nStart, nEnd});
  BOOST_REQUIRE_EQUAL(nSegments.size(), 2);
  BOOST_CHECK_EQUAL(nSegments[0].bin[0], 0);
  BOOST_CHECK_EQUAL(nSegments[1].bin[0], nRPhi - 1);

  // Test: an inclined segment across the seam is channelised like the same
  // segment rotated by half a turn, with the rPhi bins shifted by nRPhi / 2
  Vector2 iStart(halfRPhi - 0.83, -0.47);
  Vector2 iEnd(halfRPhi + 0.61, 1.38);
  Vector2 shift(halfRPhi, 0.);
  auto iSegments =
      cl.segments(geoCtx, *cylinder, *segmentation, {iStart, iEnd});
  auto iShiftedSegments = cl.segments(geoCtx, *cylinder, *segmentation,
                                      {iStart - shift, iEnd - shift});
  BOOST_REQUIRE_EQUAL(iSegments.size(), iShiftedSegments.size());
  BOOST_CHECK_LT(iSegments.size(), 15u);
  for (std::size_t i = 0; i < iSegments.size(); ++i) {
    BOOST_CHECK_EQUAL(iSegments[i].bin[0],
                      (iShiftedSegments[i].bin[0] + nRPhi / 2) % nRPhi);
    BOOST_CHECK_EQUAL(iSegments[i].bin[1], iShiftedSegments[i].bin[1]);
    BOOST_CHECK_CLOSE(iSegments[i].activation, iShiftedSegments[i].activation,
                      1e-6);
  }
  BOOST_CHECK_CLOSE(totalPath(iSegments), (iEnd - iStart).norm(), 1e-6);
}

/// Unit test for testing the Segmentizer
BOOST_DATA_TEST_CASE(
    RandomSegmentizerTest,
    bdata::random((
        bdata::engine = std::mt19937(), bdata::seed = 1,
        bdata::distribution = std::uniform_real_distribution<double>(0., 1.))) ^
        bdata::random((bdata::engine = std::mt19937(), bdata::seed = 2,
                       bdata::distribution =
                           std::uniform_real_distribution<double>(0., 1.))) ^
        bdata::random((bdata::engine = std::mt19937(), bdata::seed = 3,
                       bdata::distribution =
                           std::uniform_real_distribution<double>(0., 1.))) ^
        bdata::random((bdata::engine = std::mt19937(), bdata::seed = 4,
                       bdata::distribution =
                           std::uniform_real_distribution<double>(0., 1.))) ^
        bdata::xrange(25),
    startR0, startR1, endR0, endR1, index) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  Segmentizer cl;

  // Test beds with random numbers generated inside
  PlanarSurfaceTestBeds pstd;
  auto testBeds = pstd(1.);

  DigitizationCsvOutput csvHelper;

  for (const auto& tb : testBeds) {
    const auto& name = std::get<0>(tb);
    const auto* surface = (std::get<1>(tb)).get();
    const auto& segmentation = std::get<2>(tb);
    const auto& randomizer = std::get<3>(tb);

    if (index == 0) {
      std::ofstream shape;
      std::ofstream grid;
      const auto centerXY = surface->center(geoCtx).segment<2>(0);
      // 0 - write the shape
      shape.open("Segmentizer" + name + "Borders.csv");
      if (surface->type() == Surface::Plane) {
        const auto* pBounds =
            static_cast<const PlanarBounds*>(&(surface->bounds()));
        csvHelper.writePolygon(shape, pBounds->vertices(1), -centerXY);
      } else if (surface->type() == Surface::Disc) {
        const auto* dBounds =
            static_cast<const DiscBounds*>(&(surface->bounds()));
        csvHelper.writePolygon(shape, dBounds->vertices(72), -centerXY);
      }
      // 1 - write the grid
      grid.open("Segmentizer" + name + "Grid.csv");
      const IAxis& axis0 = segmentation->getAxis(0);
      const IAxis& axis1 = segmentation->getAxis(1);
      if (axis0.getDirection() == AxisDirection::AxisX &&
          axis1.getDirection() == AxisDirection::AxisY) {
        double bxmin = axis0.getMin();
        double bxmax = axis0.getMax();
        double bymin = axis1.getMin();
        double bymax = axis1.getMax();
        const std::vector<double> xboundaries = axis0.getBinEdges();
        const std::vector<double> yboundaries = axis1.getBinEdges();
        for (const double xval : xboundaries) {
          csvHelper.writeLine(grid, {xval, bymin}, {xval, bymax});
        }
        for (const double yval : yboundaries) {
          csvHelper.writeLine(grid, {bxmin, yval}, {bxmax, yval});
        }
      } else if (axis0.getDirection() == AxisDirection::AxisR &&
                 axis1.getDirection() == AxisDirection::AxisPhi) {
        double brmin = axis0.getMin();
        double brmax = axis0.getMax();
        double bphimin = axis1.getMin();
        double bphimax = axis1.getMax();
        const std::vector<double> rboundaries = axis0.getBinEdges();
        const std::vector<double> phiboundaries = axis1.getBinEdges();
        for (const double r : rboundaries) {
          csvHelper.writeArc(grid, r, bphimin, bphimax);
        }
        for (const double phi : phiboundaries) {
          double cphi = std::cos(phi);
          double sphi = std::sin(phi);
          csvHelper.writeLine(grid, {brmin * cphi, brmin * sphi},
                              {brmax * cphi, brmax * sphi});
        }
      }
    }

    auto start = randomizer(startR0, startR1);
    auto end = randomizer(endR0, endR1);

    std::ofstream segments;
    segments.open("Segmentizer" + name + "Segments_n" + std::to_string(index) +
                  ".csv");

    std::ofstream cluster;
    cluster.open("Segmentizer" + name + "Cluster_n" + std::to_string(index) +
                 ".csv");

    /// Run the Segmentizer
    auto cSegments = cl.segments(geoCtx, *surface, *segmentation, {start, end});

    for (const auto& cs : cSegments) {
      csvHelper.writeLine(segments, cs.path2D[0], cs.path2D[1]);
    }

    segments.close();
    cluster.close();
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
