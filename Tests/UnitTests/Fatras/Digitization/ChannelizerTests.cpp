// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsFatras/Digitization/Channelizer.hpp"
#include "BoundRandomValues.hpp"

#include <array>
#include <fstream>
#include <functional>
#include <vector>

namespace bdata = boost::unit_test::data;

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(ChannelizerCartesian) {
  Acts::GeometryContext geoCtx;

  auto rectangleBounds = std::make_shared<Acts::RectangleBounds>(1., 1.);
  auto planeSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), rectangleBounds);

  // The segementation
  Acts::BinUtility pixelated(20, -1., 1., Acts::open, Acts::binX);
  pixelated += Acts::BinUtility(20, -1., 1., Acts::open, Acts::binY);

  Channelizer cl;

  // Test: Normal hit into the surface
  Acts::Vector2D nPosition(0.37, 0.76);
  auto nSegments =
      cl.segments(geoCtx, *planeSurface, pixelated, nPosition, nPosition);
  BOOST_CHECK(nSegments.size() == 1);
  BOOST_CHECK(nSegments[0].bin[0] == 13);
  BOOST_CHECK(nSegments[0].bin[1] == 17);

  // Test: Inclined hit into the surface - negative x direction
  Acts::Vector2D ixPositionS(0.37, 0.76);
  Acts::Vector2D ixPositionE(0.02, 0.73);
  auto ixSegments =
      cl.segments(geoCtx, *planeSurface, pixelated, ixPositionS, ixPositionE);
  BOOST_CHECK(ixSegments.size() == 4);

  // Test: Inclined hit into the surface - positive y direction
  Acts::Vector2D iyPositionS(0.37, 0.76);
  Acts::Vector2D iyPositionE(0.39, 0.91);
  auto iySegments =
      cl.segments(geoCtx, *planeSurface, pixelated, iyPositionS, iyPositionE);
  BOOST_CHECK(iySegments.size() == 3);

  // Test: Inclined hit into the surface - x/y direction
  Acts::Vector2D ixyPositionS(-0.27, 0.76);
  Acts::Vector2D ixyPositionE(-0.02, -0.73);
  auto ixySegments =
      cl.segments(geoCtx, *planeSurface, pixelated, ixyPositionS, ixyPositionE);
  BOOST_CHECK(ixySegments.size() == 18);
}

BOOST_AUTO_TEST_CASE(ChannelizerPolarRadial) {
  Acts::GeometryContext geoCtx;

  auto radialBounds =
      std::make_shared<const Acts::RadialBounds>(5., 10., 0.25, 0.);
  auto radialDisc = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), radialBounds);

  // The segementation
  Acts::BinUtility strips(2, 5., 10., Acts::open, Acts::binR);
  strips += Acts::BinUtility(250, -0.25, 0.25, Acts::open, Acts::binPhi);

  Channelizer cl;

  // Test: Normal hit into the surface
  Acts::Vector2D nPosition(6.76, 0.5);
  auto nSegments =
      cl.segments(geoCtx, *radialDisc, strips, nPosition, nPosition);
  BOOST_CHECK(nSegments.size() == 1);
  BOOST_CHECK(nSegments[0].bin[0] == 0);
  BOOST_CHECK(nSegments[0].bin[1] == 161);

  // Test: now opver more phi strips
  Acts::Vector2D sPositionS(6.76, 0.5);
  Acts::Vector2D sPositionE(7.03, -0.3);
  auto sSegment =
      cl.segments(geoCtx, *radialDisc, strips, sPositionS, sPositionE);
  BOOST_CHECK(sSegment.size() == 59);

  // Test: jump over R boundary, but stay in phi bin
  sPositionS = Acts::Vector2D(6.76, 0.);
  sPositionE = Acts::Vector2D(7.83, 0.);
  sSegment = cl.segments(geoCtx, *radialDisc, strips, sPositionS, sPositionE);
  BOOST_CHECK(sSegment.size() == 2);
}

std::vector<std::array<std::ofstream, 3>> out;

/// Unit test for testing the Channelizer
BOOST_DATA_TEST_CASE(RandomChannelizerTest,
                     bdata::random(0., 1.) ^ bdata::random(0., 1.) ^
                         bdata::random(0., 1.) ^ bdata::random(0., 1.) ^
                         bdata::xrange(25),
                     startR0, startR1, endR0, endR1, index) {
  Acts::GeometryContext geoCtx;
  Channelizer cl;

  // Pixel test in Rectangle
  auto rectangle = std::make_shared<Acts::RectangleBounds>(3., 6.5);
  auto rSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), rectangle);
  Acts::BinUtility pixelated(15, -3., 3., Acts::open, Acts::binX);
  pixelated += Acts::BinUtility(26, -6.5, 6.5, Acts::open, Acts::binY);
  RectangleRandom rRandom(3., 6.5);

  // Cartesian strip test in Trapezoid
  auto trapezoid = std::make_shared<Acts::TrapezoidBounds>(2., 3.5, 4.);
  auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), trapezoid);
  Acts::BinUtility stripsX(35, -3.5, 3.5, Acts::open, Acts::binX);
  stripsX += Acts::BinUtility(1, -4., 4., Acts::open, Acts::binY);
  TrapezoidRandom tRandom(2., 3.5, 4.);

  // Phi strip test in DiscTrapezoid
  double rmin = 2.;
  double rmax = 7.5;
  double xmin = 2.;
  double xmax = 3.5;
  double ymax = std::sqrt(rmax * rmax - xmax * xmax);
  double alpha = std::max(atan2(xmin, rmin), atan2(xmax, ymax));

  auto discTrapezoid =
      std::make_shared<Acts::DiscTrapezoidBounds>(xmin, xmax, rmin, rmax);
  auto dtSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), discTrapezoid);
  Acts::BinUtility stripsPhi(1, 2., 7.5, Acts::open, Acts::binR);
  stripsPhi += Acts::BinUtility(25, M_PI_2 - alpha, M_PI_2 + alpha, Acts::open,
                                Acts::binPhi);
  TrapezoidRandom dtRandom(xmin, xmax, rmin, ymax);

  // Raidal disc test
  auto discRadial =
      std::make_shared<Acts::RadialBounds>(rmin, rmax, M_PI_4, M_PI_2);
  auto dSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), discRadial);
  Acts::BinUtility rphiseg(10, rmin, rmax, Acts::open, Acts::binR);
  rphiseg += Acts::BinUtility(20, (M_PI_2 - M_PI_4), (M_PI_2 + M_PI_4),
                              Acts::open, Acts::binPhi);
  DiscRandom dRandom(rmin, rmax, M_PI_2 - M_PI_4, M_PI_2 + M_PI_4);

  // Annulus disc test
  Acts::Vector2D aorigin(0.05, -0.1);
  double phimin = -0.25;
  double phimax = 0.38;
  auto annulus = std::make_shared<Acts::AnnulusBounds>(rmin, rmax, -0.25, 0.38,
                                                       aorigin, M_PI_4);
  auto aSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), annulus);
  Acts::BinUtility stripsPhiA(1, 0., 10., Acts::open, Acts::binR);
  stripsPhiA += Acts::BinUtility(12, -0.25, 0.36, Acts::open, Acts::binPhi);
  AnnulusRandom aRandom(rmin, rmax, phimin, phimax, aorigin.x(), aorigin.y());

  using Randomizer = std::function<std::array<double, 2>(double, double)>;
  using TestBed = std::tuple<std::string, const Acts::Surface*,
                             Acts::BinUtility, Randomizer>;

  std::vector<TestBed> testBeds = {
      {"Rectangle", rSurface.get(), pixelated, rRandom},
      {"Trapezoid", tSurface.get(), stripsX, tRandom},
      {"DiscTrapezoid", dtSurface.get(), stripsPhi, dtRandom},
      {"DiscRadial", dSurface.get(), rphiseg, dRandom},
      {"Annulus", aSurface.get(), stripsPhiA, aRandom}};

  for (const auto& tb : testBeds) {
    const auto& name = std::get<0>(tb);
    const auto* surface = std::get<1>(tb);
    const auto& segmentation = std::get<2>(tb);
    const auto& randomizer = std::get<3>(tb);

    if (index == 0) {
      std::ofstream shape;
      std::ofstream grid;

      // Helper method to write a line
      auto writeLine = [](std::ostream& outf, const Acts::Vector2D& p0,
                          const Acts::Vector2D& p1) -> void {
        outf << "l," << p0.x() << "," << p0.y() << "," << p1.x() << ","
             << p1.y() << "\n";
      };

      // Helper method to write an arc
      auto writeArc = [&](std::ostream& outf, double r, double phiMin,
                          double phiMax) -> void {
        outf << "a," << r << "," << r << "," << phiMin << "," << phiMax << "\n";
      };

      // Helper method to write a polygon
      auto writePolygon =
          [&](std::ostream& outf,
              const std::vector<Acts::Vector2D>& vertices) -> void {
        auto nvertices = vertices.size();
        for (unsigned long iv = 1; iv < nvertices; ++iv) {
          writeLine(outf, vertices[iv - 1], vertices[iv]);
        }
        writeLine(outf, vertices[nvertices - 1], vertices[0]);
      };

      // 0 - write the shape
      shape.open("Channelizer" + name + "Borders.csv");
      if (surface->type() == Acts::Surface::Plane) {
        const auto* pBounds =
            static_cast<const Acts::PlanarBounds*>(&(surface->bounds()));
        writePolygon(shape, pBounds->vertices(1));
      } else if (surface->type() == Acts::Surface::Disc) {
        const auto* dBounds =
            static_cast<const Acts::DiscBounds*>(&(surface->bounds()));
        writePolygon(shape, dBounds->vertices(72));
      }
      // 1 - write the grid
      grid.open("Channelizer" + name + "Grid.csv");
      if (segmentation.binningData()[0].binvalue == Acts::binX and
          segmentation.binningData()[1].binvalue == Acts::binY) {
        double bxmin = segmentation.binningData()[0].min;
        double bxmax = segmentation.binningData()[0].max;
        double bymin = segmentation.binningData()[1].min;
        double bymax = segmentation.binningData()[1].max;
        const auto& xboundaries = segmentation.binningData()[0].boundaries();
        const auto& yboundaries = segmentation.binningData()[1].boundaries();
        for (const auto xval : xboundaries) {
          writeLine(grid, {xval, bymin}, {xval, bymax});
        }
        for (const auto yval : yboundaries) {
          writeLine(grid, {bxmin, yval}, {bxmax, yval});
        }
      } else if (segmentation.binningData()[0].binvalue == Acts::binR and
                 segmentation.binningData()[1].binvalue == Acts::binPhi) {
        double brmin = 0.;  // segmentation.binningData()[0].min;
        double brmax = segmentation.binningData()[0].max;
        double bphimin = segmentation.binningData()[1].min;
        double bphimax = segmentation.binningData()[1].max;
        const auto& rboundaries = segmentation.binningData()[0].boundaries();
        const auto& phiboundaries = segmentation.binningData()[1].boundaries();
        for (const auto r : rboundaries) {
          writeArc(grid, r, bphimin, bphimax);
        }
        for (const auto phi : phiboundaries) {
          double cphi = std::cos(phi);
          double sphi = std::sin(phi);
          writeLine(grid, {brmin * cphi, brmin * sphi},
                    {brmax * cphi, brmax * sphi});
        }
      }
    }

    auto start = randomizer(startR0, startR1);
    auto end = randomizer(endR0, endR1);

    std::ofstream segments;
    segments.open("Channelizer" + name + "Segments_n" + std::to_string(index) +
                  ".csv");

    std::ofstream cluster;
    cluster.open("Channelizer" + name + "Cluster_n" + std::to_string(index) +
                 ".csv");

    /// Run the channelizer
    auto cSegement = cl.segments(geoCtx, *surface, segmentation,
                                 {start[0], start[1]}, {end[0], end[1]});

    for (const auto& cs : cSegement) {
      segments << "l," << cs.path2D[0].x() << "," << cs.path2D[0].y() << ","
               << cs.path2D[1].x() << "," << cs.path2D[1].y() << "\n";
    }

    segments.close();
    cluster.close();
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras