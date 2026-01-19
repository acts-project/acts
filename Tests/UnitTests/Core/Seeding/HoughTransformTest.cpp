// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <format>
#include <vector>

using namespace Acts;

namespace ActsTests {

auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

struct DriftCircle {
  double y{0.};
  double z{0.};
  double rDrift{0.};
  double rDriftError{0.};

  DriftCircle(const double _y, const double _z, const double _r,
              const double _rUncert)
      : y{_y}, z{_z}, rDrift{_r}, rDriftError{_rUncert} {}
};

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(hough_transform_seeder) {
  // we are using the slope on yz plane with the y coordinate (hardcoded from
  // the csv MuonSimHit data)
  std::vector<std::pair<double, double>> simHits = {
      {-0.0401472 / 0.994974, -422.612}};

  // Define the drift Circles
  constexpr double uncert{0.3};
  std::array<DriftCircle, 6> driftCircles{
      DriftCircle{-427.981, -225.541, 14.5202, uncert},
      DriftCircle{-412.964, -199.53, 1.66237, uncert},
      DriftCircle{-427.981, -173.519, 12.3176, uncert},
      DriftCircle{-427.981, 173.519, 1.5412, uncert},
      DriftCircle{-442.999, 199.53, 12.3937, uncert},
      DriftCircle{-427.981, 225.541, 3.77967, uncert}

  };

  // configure the binning of the hough plane
  HoughTransformUtils::HoughPlaneConfig planeCfg;
  planeCfg.nBinsX = 1000;
  planeCfg.nBinsY = 1000;

  // instantiate the peak finder
  HoughTransformUtils::PeakFinders::IslandsAroundMaxConfig peakFinderCfg;
  peakFinderCfg.fractionCutoff = 0.7;
  peakFinderCfg.threshold = 3.;
  peakFinderCfg.minSpacingBetweenPeaks = {0., 30.};

  // and map the hough plane to parameter ranges.
  // The first coordinate is tan(theta), the second is z0 in mm
  HoughTransformUtils::HoughAxisRanges axisRanges{-3., 3., -2000., 2000.};

  // create the functions parametrising the hough space lines for drift circles.
  // Note that there are two solutions for each drift circle and angle

  // left solution
  auto houghParamFromDCleft = [](double tanTheta, const DriftCircle& DC) {
    return DC.y - tanTheta * DC.z - DC.rDrift / std::cos(std::atan(tanTheta));
  };
  // right solution
  auto houghParamFromDCright = [](double tanTheta, const DriftCircle& DC) {
    return DC.y - tanTheta * DC.z + DC.rDrift / std::cos(std::atan(tanTheta));
  };

  // create the function parametrising the drift radius uncertainty
  auto houghWidthFromDC = [](double, const DriftCircle& DC) {
    return std::min(DC.rDriftError * 3.,
                    1.0);  // scale reported errors up to at least 1mm or 3
                           // times the reported error as drift circle calib not
                           // fully reliable at this stage
  };

  // instantiate the hough plane
  HoughTransformUtils::HoughPlane<GeometryIdentifier::Value> houghPlane(
      planeCfg);

  // also instantiate the peak finder
  HoughTransformUtils::PeakFinders::IslandsAroundMax<GeometryIdentifier::Value>
      peakFinder(peakFinderCfg);

  // loop over the true hits
  for (auto& sh : simHits) {
    houghPlane.reset();

    for (std::size_t k = 0; k < driftCircles.size(); ++k) {
      auto dc = driftCircles[k];

      houghPlane.fill<DriftCircle>(dc, axisRanges, houghParamFromDCleft,
                                   houghWidthFromDC, k);
      houghPlane.fill<DriftCircle>(dc, axisRanges, houghParamFromDCright,
                                   houghWidthFromDC, k);
    }

    // now get the peaks
    auto maxima = peakFinder.findPeaks(houghPlane, axisRanges);
    for (auto& max : maxima) {
      // check the Hough Transforms results
      BOOST_CHECK_CLOSE(max.x, sh.first, 4.);
      BOOST_CHECK_CLOSE(max.y, sh.second, 0.2);
    }
  }
}

BOOST_AUTO_TEST_CASE(hough_transform_sliding_window) {
  // Create a simple 10x10 Hough space and fill it with two triangular
  // (pyramid-shaped) distributions. Each distribution has its maximum
  // value at the specified bin coordinates and decreases linearly with
  // Manhattan distance. When distributions overlap we take the maximum
  // of the two contributions so that the global maximum remains the
  // requested peak value.
  const std::size_t nX = 10;
  const std::size_t nY = 10;
  HoughTransformUtils::HoughPlaneConfig config{nX, nY};
  HoughTransformUtils::HoughPlane<int> plane(config);

  auto addTriangular = [&](const std::vector<std::array<std::size_t, 2>>& peaks,
                           HoughTransformUtils::YieldType peak) {
    for (std::size_t x = 0; x < nX; ++x) {
      for (std::size_t y = 0; y < nY; ++y) {
        HoughTransformUtils::YieldType val = 0.0f;
        for (auto [cx, cy] : peaks) {
          int dist = (x >= cx ? x - cx : cx - x) +
                     (y >= cy ? y - cy : cy - y);  // Manhattan distance
          val = std::max(
              val, peak - static_cast<HoughTransformUtils::YieldType>(dist));
        }
        plane.fillBin(x, y, x, y, val);
      }
    }
  };

  // Add 3 triangular peaks with maxima 4 at specified locations
  std::vector<std::array<std::size_t, 2>> testPeaks({{4, 4}, {4, 5}, {2, 6}});
  addTriangular(testPeaks, 4.0f);

  // Verify that the maxima at the requested locations are exactly 4
  BOOST_CHECK_EQUAL(plane.nHits(4, 4), 4.0f);
  BOOST_CHECK_EQUAL(plane.nHits(4, 5), 4.0f);

  // config for unisolated max finding
  HoughTransformUtils::PeakFinders::SlidingWindowConfig cfg1{4,     0, 0,
                                                             false, 0, 0};
  auto peaks1 = slidingWindowPeaks(plane, cfg1);
  BOOST_CHECK_EQUAL(peaks1[0][0], 2);
  BOOST_CHECK_EQUAL(peaks1[0][1], 6);
  BOOST_CHECK_EQUAL(peaks1[1][0], 4);
  BOOST_CHECK_EQUAL(peaks1[1][1], 4);
  BOOST_CHECK_EQUAL(peaks1[2][0], 4);
  BOOST_CHECK_EQUAL(peaks1[2][1], 5);

  // config for removing duplicates, no recentering
  HoughTransformUtils::PeakFinders::SlidingWindowConfig cfg2{4,     2, 2,
                                                             false, 0, 0};
  auto peaks2 = slidingWindowPeaks(plane, cfg2);
  BOOST_CHECK_EQUAL(peaks2.size(), 1);
  BOOST_CHECK_EQUAL(peaks2[0][0], 4);
  BOOST_CHECK_EQUAL(peaks2[0][1], 4);

  // config for removing duplicates, with recentering
  HoughTransformUtils::PeakFinders::SlidingWindowConfig cfg3{4,    2, 2,
                                                             true, 2, 2};
  auto peaks3 = slidingWindowPeaks(plane, cfg3);
  BOOST_CHECK_EQUAL(peaks3.size(), 1);
  BOOST_CHECK_EQUAL(peaks3[0][0], 3);
  BOOST_CHECK_EQUAL(peaks3[0][1], 5);

  // test behaviour at the edge of the plane (safety check)
  plane.reset();
  addTriangular({{1, 1}, {5, 5}, {8, 9}}, 4.0f);
  peaks3 = slidingWindowPeaks(plane, cfg3);
  BOOST_CHECK_EQUAL(peaks3.size(), 1);

  {
    auto img1 =
        HoughTransformUtils::PeakFinders::hitsCountImage(plane, {1, 1}, 4, 4);
    // the image should reflect this content (peak at 1, 1)
    //      0 1 2 3 4 5 - indices
    //   +--------+
    //   |        |
    // 0 |  2 3 2 |
    // 1 |  3 4 3 |
    // 2 |  2 3 2 |
    //   +________+
    // content outside of valid plane indices is padded with zeros
    std::vector<unsigned char> expected(
        {0, 0, 0, 0, 0, 2, 3, 2, 0, 3, 4, 3, 0, 2, 3, 2});
    BOOST_CHECK_EQUAL(expected.size(), img1.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(img1.begin(), img1.end(), expected.begin(),
                                  expected.end());
  }
  {
    // customized summary function: gets layers bitmask
    auto maskGetter = [](const HoughTransformUtils::HoughPlane<int>& p, int x,
                         int y) -> unsigned {
      unsigned mask = 0;
      for (unsigned l : p.layers(x, y)) {
        mask |= 1 << std::min(l, 16u);
        std::cout << x << " " << y << " " << l << "\n";
      }
      return mask;
    };

    plane.reset();
    plane.fillBin(1, 1, 0, 1, 1);
    plane.fillBin(1, 1, 1, 2, 1);
    plane.fillBin(1, 1, 2, 3, 1);
    plane.fillBin(2, 2, 3, 3, 1);

    auto img = HoughTransformUtils::PeakFinders::hitsCountImage<int, unsigned>(
        plane, {1, 1}, 4, 4, maskGetter);

    std::vector<unsigned char> expected({0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                         1 << 0x1 | 1 << 0x2 | 1 << 0x3, 0, 0,
                                         0, 0, 1 << 0x3});

    BOOST_CHECK_EQUAL(expected.size(), img.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(img.begin(), img.end(), expected.begin(),
                                  expected.end());
  }
}

BOOST_AUTO_TEST_CASE(hough_plane_layers_hits) {
  // Create 1x1 Hough plane and fill it with hits. For each layer add one more
  // hit than to the previous one.

  const std::size_t nX = 1;
  const std::size_t nY = 1;
  HoughTransformUtils::HoughPlaneConfig config{nX, nY};
  HoughTransformUtils::HoughPlane<uint8_t> plane(config);

  uint8_t nHits = 0;
  static constexpr uint8_t nLayers = 10;
  for (uint8_t layer = 1; layer <= nLayers; ++layer) {
    // Add hits equal to the layer number
    for (uint8_t hit = 0; hit < layer; ++hit) {
      plane.fillBin(0, 0, nHits++, layer);
    }
  }

  BOOST_CHECK_EQUAL(nLayers, plane.nLayers(0, 0));
  BOOST_CHECK_EQUAL(nHits, plane.nHits(0, 0));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
