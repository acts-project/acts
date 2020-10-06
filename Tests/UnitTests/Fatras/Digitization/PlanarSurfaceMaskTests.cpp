// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "ActsFatras/Digitization/DigitizationError.hpp"
#include "ActsFatras/Digitization/detail/PlanarSurfaceMask.hpp"
#include <Acts/Surfaces/AnnulusBounds.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/DiscTrapezoidBounds.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Surfaces/TrapezoidBounds.hpp>
#include <Acts/Tests/CommonHelpers/FloatComparisons.hpp>
#include <Acts/Utilities/Definitions.hpp>

#include <fstream>

namespace bdata = boost::unit_test::data;

namespace ActsFatras {

std::vector<std::array<std::ofstream, 3>> out;

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(PlaneMaskRectangleBounds) {
  auto rectangleBounds = std::make_shared<Acts::RectangleBounds>(2., 3.5);
  auto planeSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), rectangleBounds);

  using Segment = std::pair<Acts::Vector2D, Acts::Vector2D>;

  ActsFatras::detail::PlanarSurfaceMask psm;

  /// Case one : one outside
  Segment segment = {{2.5, -4.5}, {-1., -1.}};
  auto clipped = psm.apply(*planeSurface, segment).value();

  CHECK_CLOSE_ABS(clipped.second.x(), segment.second.x(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.second.y(), segment.second.y(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.first.x(), 1.5, Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.first.y(), -3.5, Acts::s_epsilon);

  /// Case two : two outside
  segment = {{1., 4.}, {3., 2.}};
  clipped = psm.apply(*planeSurface, segment).value();

  CHECK_CLOSE_ABS(clipped.second.x(), 2., Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.second.y(), 3., Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.first.x(), 1.5, Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.first.y(), 3.5, Acts::s_epsilon);

  /// Case two : both inside (most likely case, untouched)
  segment = {{-1., 0.5}, {0., 2.}};
  clipped = psm.apply(*planeSurface, segment).value();

  CHECK_CLOSE_ABS(clipped.first.x(), segment.first.x(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.first.y(), segment.first.y(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.second.x(), segment.second.x(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.second.y(), segment.second.y(), Acts::s_epsilon);
}

BOOST_AUTO_TEST_CASE(DiscMaskRadialBounds) {
  auto discRadial =
      std::make_shared<Acts::RadialBounds>(2., 7.5, M_PI_4, M_PI_2);
  auto discSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), discRadial);

  using Segment = std::pair<Acts::Vector2D, Acts::Vector2D>;

  ActsFatras::detail::PlanarSurfaceMask psm;

  /// Case one : one outside R min
  Segment segment = {{0.5, 1.8}, {0.9, 6.}};
  auto clipped = psm.apply(*discSurface, segment).value();

  CHECK_CLOSE_ABS(clipped.second.x(), segment.second.x(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.second.y(), segment.second.y(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(Acts::VectorHelpers::perp(clipped.first), 2.,
                  5 * Acts::s_epsilon);

  /// Case two : one outside R max
  segment = {{0.5, 2.8}, {0.9, 8.5}};
  clipped = psm.apply(*discSurface, segment).value();

  CHECK_CLOSE_ABS(clipped.first.x(), segment.first.x(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(clipped.first.y(), segment.first.y(), Acts::s_epsilon);
  CHECK_CLOSE_ABS(Acts::VectorHelpers::perp(clipped.second), 7.5,
                  5 * Acts::s_epsilon);

  /// Case three : both outside R min / max
  segment = {{0.5, 1.8}, {0.9, 8.5}};
  clipped = psm.apply(*discSurface, segment).value();
  CHECK_CLOSE_ABS(Acts::VectorHelpers::perp(clipped.first), 2.,
                  5 * Acts::s_epsilon);
  CHECK_CLOSE_ABS(Acts::VectorHelpers::perp(clipped.second), 7.5,
                  5 * Acts::s_epsilon);
  /// Case four: outside phi min
  segment = {{2.8, 2.5}, {0., 3.5}};
  clipped = psm.apply(*discSurface, segment).value();
  CHECK_CLOSE_ABS(Acts::VectorHelpers::phi(clipped.first), M_PI_4,
                  Acts::s_epsilon);

  /// Case four: outside phi max
  segment = {{0., 3.5}, {-8., 5.}};
  clipped = psm.apply(*discSurface, segment).value();
  CHECK_CLOSE_ABS(Acts::VectorHelpers::phi(clipped.second), M_PI_2 + M_PI_4,
                  Acts::s_epsilon);
}

/// Unit test for testing the orientedSurfaces() function
BOOST_DATA_TEST_CASE(RandomMaskingTest,
                     bdata::random(-10., 10.) ^ bdata::random(-10., 10.) ^
                         bdata::random(-10., 10.) ^ bdata::random(-10., 10.) ^
                         bdata::xrange(1000),
                     startX, startY, endX, endY, index) {
  auto rectangle = std::make_shared<Acts::RectangleBounds>(3., 6.5);
  auto rSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), rectangle);

  auto trapezoid = std::make_shared<Acts::TrapezoidBounds>(2., 3.5, 4.);
  auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3D::Identity(), trapezoid);

  auto discTrapezoid =
      std::make_shared<Acts::DiscTrapezoidBounds>(2., 3.5, 2., 7.5);
  auto dtSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), discTrapezoid);

  auto discRing = std::make_shared<Acts::RadialBounds>(2., 8.5);
  auto driSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), discRing);

  auto discRadial =
      std::make_shared<Acts::RadialBounds>(2., 7.5, M_PI_4, M_PI_2);
  auto dSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), discRadial);

  auto annulus = std::make_shared<Acts::AnnulusBounds>(3., 8.5, 0.75, 1.8);
  auto aSurface = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3D::Identity(), annulus);

  std::vector<std::pair<std::string, const Acts::Surface*>> identifiedSurfaces =
      {{"Rectangle", rSurface.get()},
       {"Trapezoid", tSurface.get()},
       {"DiscTrapezoid", dtSurface.get()},
       {"DiscRadial", dSurface.get()},
       {"Annulus", aSurface.get()}};

  using Segment = std::pair<Acts::Vector2D, Acts::Vector2D>;
  ActsFatras::detail::PlanarSurfaceMask psm;

  unsigned int io = 0;
  for (const auto ids : identifiedSurfaces) {
    if (index == 0) {
      out.push_back(std::array<std::ofstream, 3>());
      out[io][0].open(ids.first + "Inside.csv");
      out[io][1].open(ids.first + "Clipped.csv");
      out[io][2].open(ids.first + "Outside.csv");
      for (size_t ifs = 0; ifs < 3; ++ifs) {
        out[io][ifs] << "x0,y0,x1,y1\n";
      }
    }

    /// Case one : one outside
    Segment segment = {{startX, startY}, {endX, endY}};
    auto clipped = psm.apply(*ids.second, segment);

    if (clipped.ok()) {
      Segment csegment = clipped.value();
      out[io][0] << csegment.first.x() << "," << csegment.first.y() << ",";
      out[io][0] << csegment.second.x() << "," << csegment.second.y() << "\n";
      if (segment != csegment) {
        if (not segment.first.isApprox(csegment.first) and
            segment.second.isApprox(csegment.second)) {
          // one side clipped
          out[io][1] << segment.first.x() << "," << segment.first.y() << ",";
          out[io][1] << csegment.first.x() << "," << csegment.first.y() << "\n";
        } else if (segment.first.isApprox(csegment.first)) {
          // second side clipped
          out[io][1] << csegment.second.x() << "," << csegment.second.y()
                     << ",";
          out[io][1] << segment.second.x() << "," << segment.second.y() << "\n";
        }
      }
    } else {
      out[io][2] << segment.first.x() << "," << segment.first.y() << ",";
      out[io][2] << segment.second.x() << "," << segment.second.y() << "\n";
    }

    if (index == 999) {
      for (size_t ifs = 0; ifs < 3; ++ifs) {
        out[io][ifs].close();
      }
    }
    io++;
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras