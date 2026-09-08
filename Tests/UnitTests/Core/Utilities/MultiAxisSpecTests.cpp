// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <numbers>
#include <optional>
#include <stdexcept>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(MultiAxisSpecBasics) {
  using enum Acts::AxisBoundaryType;

  MultiAxisSpec maf({AxisSpec::Equidistant(10, 0., 1., Bound),
                     AxisSpec::Variable({0., 1., 3.}, Open)});
  BOOST_CHECK_EQUAL(maf.size(), 2);
  BOOST_CHECK(!maf.isDeferred());
  BOOST_CHECK_EQUAL(maf.axisSpec(0).nBins(), 10);
  BOOST_CHECK_EQUAL(maf.axisSpec(1).nBins(), 2);
  BOOST_CHECK_THROW(maf.axisSpec(2), std::out_of_range);

  auto multiAxis = maf.buildMultiAxis();
  BOOST_CHECK_EQUAL(multiAxis->getNAxes(), 2);
  BOOST_CHECK(multiAxis->getAxis(0).isEquidistant());
  BOOST_CHECK(multiAxis->getAxis(1).isVariable());
  BOOST_CHECK_EQUAL(multiAxis->getAxis(0).getNBins(), 10);
  BOOST_CHECK_EQUAL(multiAxis->getAxis(1).getNBins(), 2);

  // Empty construction is invalid
  BOOST_CHECK_THROW(MultiAxisSpec({}), std::invalid_argument);

  // Mixed with a deferred axis: full resolution is required for all axes
  MultiAxisSpec mixed({AxisSpec::Equidistant(10, 0., 1., Bound),
                       AxisSpec::DeferredEquidistant(5)});
  BOOST_CHECK(mixed.isDeferred());
  BOOST_CHECK_THROW(mixed.buildMultiAxis(), std::domain_error);
}

BOOST_AUTO_TEST_CASE(MultiAxisSpecDeferredResolution) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;

  MultiAxisSpec maf(
      {AxisSpec::DeferredEquidistant(4, AxisRPhi),
       AxisSpec::DeferredVariable({0., 0.5, 1.}, std::nullopt, AxisZ)});
  BOOST_CHECK(maf.isDeferred());
  BOOST_CHECK_THROW(maf.buildMultiAxis(), std::domain_error);

  MultiAxisSpec::Options options = {
      {.min = -3., .max = 3., .boundaryType = Closed, .direction = AxisRPhi},
      {.min = -10., .max = 10., .boundaryType = Bound, .direction = AxisZ}};
  auto axes = maf.buildMultiAxis(options);
  BOOST_CHECK_EQUAL(axes->getNAxes(), 2);
  BOOST_CHECK_EQUAL(axes->getAxis(0).getBoundaryType(), Closed);
  BOOST_CHECK_EQUAL(axes->getAxis(0).getNBins(), 4);
  BOOST_CHECK(axes->getAxis(0).getDirection() == AxisRPhi);
  BOOST_CHECK_EQUAL(axes->getAxis(1).getBoundaryType(), Bound);
  CHECK_CLOSE_ABS(axes->getAxis(1).getBinEdges()[1], 0., 1e-15);
  BOOST_CHECK(axes->getAxis(1).getDirection() == AxisZ);

  // Directions are validated per axis
  BOOST_CHECK_THROW(
      maf.buildMultiAxis(
          {{.min = -3., .max = 3., .boundaryType = Closed, .direction = AxisZ},
           {.min = -10.,
            .max = 10.,
            .boundaryType = Bound,
            .direction = AxisRPhi}}),
      std::invalid_argument);

  // Sizes have to match
  BOOST_CHECK_THROW(
      maf.buildMultiAxis({{.min = 0., .max = 1., .boundaryType = Bound}}),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MultiAxisSpecXDApi) {
  using enum Acts::AxisBoundaryType;

  Acts::MultiAxisSpec1D maf1D({AxisSpec::DeferredEquidistant(8)});
  BOOST_CHECK_EQUAL(maf1D.size(), 1);
  std::unique_ptr<Acts::IMultiAxis1D> ma1D = maf1D.buildMultiAxis(
      {AxisSpec::Options{.min = 0., .max = 4., .boundaryType = Bound}});
  BOOST_CHECK_EQUAL(ma1D->getNAxes(), 1);
  BOOST_CHECK_EQUAL(ma1D->getAxis(0).getNBins(), 8);

  Acts::MultiAxisSpec2D maf2D({AxisSpec::Equidistant(2, 0., 1., Bound),
                               AxisSpec::Equidistant(3, 0., 1., Bound)});
  std::unique_ptr<Acts::IMultiAxis2D> ma2D = maf2D.buildMultiAxis();
  BOOST_CHECK_EQUAL(ma2D->getNAxes(), 2);
  BOOST_CHECK_EQUAL(ma2D->getNTotalBins(), 6);

  // The XD types remain usable through the runtime-dimension base
  const MultiAxisSpec& base = maf2D;
  BOOST_CHECK_EQUAL(base.buildMultiAxis()->getNAxes(), 2);
}

BOOST_AUTO_TEST_CASE(MultiAxisSpecEqualityAndStreams) {
  using enum Acts::AxisBoundaryType;

  MultiAxisSpec a(
      {AxisSpec::DeferredEquidistant(4), AxisSpec::DeferredEquidistant(5)});
  MultiAxisSpec b(
      {AxisSpec::DeferredEquidistant(4), AxisSpec::DeferredEquidistant(5)});
  MultiAxisSpec c({AxisSpec::DeferredEquidistant(4)});
  BOOST_CHECK(a == b);
  BOOST_CHECK(a != c);

  BOOST_CHECK_EQUAL(
      c.toString(),
      "MultiAxisSpec: 1 axes [AxisSpec: 4 bins, equidistant within "
      "deferred range, deferred boundary type]");
}

BOOST_AUTO_TEST_CASE(ResolveMultiAxisAgainstSurface) {
  using enum AxisDirection;
  using enum AxisBoundaryType;

  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(30., 100.));

  // Positional matching without directions
  MultiAxisSpec2D positional(
      {AxisSpec::DeferredEquidistant(10), AxisSpec::DeferredEquidistant(20)});
  auto axes = resolveMultiAxis(positional, *cylinder);
  BOOST_CHECK_EQUAL(axes->getNAxes(), 2);
  BOOST_CHECK(axes->getAxis(0).getDirection() == AxisRPhi);
  BOOST_CHECK_EQUAL(axes->getAxis(0).getBoundaryType(), Closed);
  BOOST_CHECK_EQUAL(axes->getAxis(0).getNBins(), 10);
  BOOST_CHECK(axes->getAxis(1).getDirection() == AxisZ);
  BOOST_CHECK_EQUAL(axes->getAxis(1).getBoundaryType(), Bound);
  CHECK_CLOSE_ABS(axes->getAxis(1).getMin(), -100., 1e-12);

  // Directed matching in canonical order
  MultiAxisSpec2D directed({AxisSpec::DeferredEquidistant(10, AxisRPhi),
                            AxisSpec::DeferredEquidistant(20, AxisZ)});
  auto directedAxes = resolveMultiAxis(directed, *cylinder);
  BOOST_CHECK_EQUAL(directedAxes->getAxis(0).getNBins(), 10);
  BOOST_CHECK_EQUAL(directedAxes->getAxis(1).getNBins(), 20);

  // Swapped order is re-ordered to match the surface
  MultiAxisSpec2D swapped({AxisSpec::DeferredEquidistant(20, AxisZ),
                           AxisSpec::DeferredEquidistant(10, AxisRPhi)});
  auto reordered = resolveMultiAxis(swapped, *cylinder);
  BOOST_CHECK(reordered->getAxis(0).getDirection() == AxisRPhi);
  BOOST_CHECK_EQUAL(reordered->getAxis(0).getNBins(), 10);
  BOOST_CHECK(reordered->getAxis(1).getDirection() == AxisZ);
  BOOST_CHECK_EQUAL(reordered->getAxis(1).getNBins(), 20);

  // Binning along a single direction is a 2D binning with one bin in the
  // other direction
  MultiAxisSpec2D oneD({AxisSpec::DeferredEquidistant(1, AxisRPhi),
                        AxisSpec::DeferredEquidistant(8, AxisZ)});
  auto oneDAxes = resolveMultiAxis(oneD, *cylinder);
  BOOST_CHECK_EQUAL(oneDAxes->getNAxes(), 2);
  BOOST_CHECK_EQUAL(oneDAxes->getAxis(0).getNBins(), 1);
  BOOST_CHECK_EQUAL(oneDAxes->getAxis(1).getNBins(), 8);

  // Deferred variable binning is scaled onto the surface range
  MultiAxisSpec2D variable(
      {AxisSpec::DeferredEquidistant(1, AxisRPhi),
       AxisSpec::DeferredVariable({0., 0.25, 1.}, std::nullopt, AxisZ)});
  auto variableAxes = resolveMultiAxis(variable, *cylinder);
  auto edges = variableAxes->getAxis(1).getBinEdges();
  CHECK_CLOSE_ABS(edges[0], -100., 1e-12);
  CHECK_CLOSE_ABS(edges[1], -50., 1e-12);
  CHECK_CLOSE_ABS(edges[2], 100., 1e-12);

  // Mismatching directions are rejected
  MultiAxisSpec2D wrongDirs({AxisSpec::DeferredEquidistant(10, AxisR),
                             AxisSpec::DeferredEquidistant(20, AxisPhi)});
  BOOST_CHECK_THROW(resolveMultiAxis(wrongDirs, *cylinder),
                    std::invalid_argument);

  // Duplicate directions are rejected
  MultiAxisSpec2D duplicate({AxisSpec::DeferredEquidistant(10, AxisZ),
                             AxisSpec::DeferredEquidistant(20, AxisZ)});
  BOOST_CHECK_THROW(resolveMultiAxis(duplicate, *cylinder),
                    std::invalid_argument);

  // Mixed directed and undirected axes are rejected
  MultiAxisSpec2D mixed({AxisSpec::DeferredEquidistant(10, AxisRPhi),
                         AxisSpec::DeferredEquidistant(20)});
  BOOST_CHECK_THROW(resolveMultiAxis(mixed, *cylinder), std::invalid_argument);

  // Unsupported surface
  auto straw =
      Surface::makeShared<StrawSurface>(Transform3::Identity(), 5., 100.);
  MultiAxisSpec2D strawBinning(
      {AxisSpec::DeferredEquidistant(4), AxisSpec::DeferredEquidistant(4)});
  BOOST_CHECK_THROW(resolveMultiAxis(strawBinning, *straw),
                    std::invalid_argument);

  // A spec that fixes what the surface dictates has to agree with it
  MultiAxisSpec2D full({AxisSpec::Equidistant(10, 0., 1., Bound),
                        AxisSpec::DeferredEquidistant(20)});
  BOOST_CHECK_THROW(resolveMultiAxis(full, *cylinder), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ResolveMultiAxisSurfaceRanges) {
  using enum AxisBoundaryType;
  using enum AxisDirection;

  // Undirected axes bind positionally to the canonical surface directions
  MultiAxisSpec2D binning(
      {AxisSpec::DeferredEquidistant(1), AxisSpec::DeferredEquidistant(1)});

  // Disc: bound in r, closed in phi over the full azimuth
  auto disc = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), std::make_shared<RadialBounds>(50., 75.));
  std::array<AxisDirection, 2> expectedDisc = {AxisR, AxisPhi};
  BOOST_CHECK(static_cast<const Surface&>(*disc).localAxes() == expectedDisc);
  auto discAxes = resolveMultiAxis(binning, *disc);
  CHECK_CLOSE_ABS(discAxes->getAxis(0).getMin(), 50., 1e-12);
  CHECK_CLOSE_ABS(discAxes->getAxis(0).getMax(), 75., 1e-12);
  BOOST_CHECK_EQUAL(discAxes->getAxis(0).getBoundaryType(), Bound);
  BOOST_CHECK_EQUAL(discAxes->getAxis(1).getBoundaryType(), Closed);

  // Rectangle
  auto rectangle = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(20., 30.));
  std::array<AxisDirection, 2> expectedPlane = {AxisX, AxisY};
  BOOST_CHECK(static_cast<const Surface&>(*rectangle).localAxes() ==
              expectedPlane);
  auto rectangleAxes = resolveMultiAxis(binning, *rectangle);
  CHECK_CLOSE_ABS(rectangleAxes->getAxis(0).getMin(), -20., 1e-12);
  CHECK_CLOSE_ABS(rectangleAxes->getAxis(0).getMax(), 20., 1e-12);
  BOOST_CHECK_EQUAL(rectangleAxes->getAxis(0).getBoundaryType(), Bound);

  // Trapezoid: x spans the wider of the two half lengths
  auto trapezoid = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<TrapezoidBounds>(5., 15., 30.));
  auto trapezoidAxes = resolveMultiAxis(binning, *trapezoid);
  CHECK_CLOSE_ABS(trapezoidAxes->getAxis(0).getMin(), -15., 1e-12);
  CHECK_CLOSE_ABS(trapezoidAxes->getAxis(0).getMax(), 15., 1e-12);
  CHECK_CLOSE_ABS(trapezoidAxes->getAxis(1).getMax(), 30., 1e-12);

  // Sectoral cylinder: the azimuthal axis is bound, not closed
  auto sector = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(),
      std::make_shared<CylinderBounds>(25., 50., std::numbers::pi / 4.));
  auto sectorAxes = resolveMultiAxis(binning, *sector);
  CHECK_CLOSE_ABS(sectorAxes->getAxis(0).getMin(), -25. * std::numbers::pi / 4.,
                  1e-12);
  BOOST_CHECK_EQUAL(sectorAxes->getAxis(0).getBoundaryType(), Bound);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
