// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "ActsExamples/EventData/DriftCircle.hpp"
#include "ActsExamples/EventData/MuonSimHit.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <array>
#include <memory>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(hough_transform_seeder) {
  // we are using the slope on yz plane with the y coordinate (hardcoded from
  // the csv MuonSimHit data)
  std::vector<std::pair<double, double>> simHits = {
      {-0.0401472 / 0.994974, -422.612}};
  // The drift circles info (hardcoded from the csv file)
  std::vector<std::vector<float>> tubePos{{-4.54747e-13, -427.981, -225.541},
                                          {0, -412.964, -199.53},
                                          {1.36424e-12, -427.981, -173.519},
                                          {-9.09495e-13, -427.981, 173.519},
                                          {0, -442.999, 199.53},
                                          {-4.54747e-13, -427.981, 225.541}};
  std::vector<float> driftRadius{14.5202, 1.66237, 12.3176,
                                 1.5412,  12.3937, 3.77967};
  std::vector<std::vector<int8_t>> stationInfo{
      {5, -1, 4, 1, 1, 39}, {5, -1, 4, 1, 2, 38}, {5, -1, 4, 1, 3, 39},
      {5, -1, 4, 2, 1, 39}, {5, -1, 4, 2, 2, 39}, {5, -1, 4, 2, 3, 39}};

  BOOST_CHECK_EQUAL(tubePos.size(), driftRadius.size());
  BOOST_CHECK_EQUAL(stationInfo.size(), tubePos.size());
  // Define the drift Circles
  std::vector<ActsExamples::DriftCircle> driftCircles;
  for (std::size_t i = 0; i < tubePos.size(); i++) {
    ActsFatras::Hit::Vector3 tube_pos{tubePos[i][0] * Acts::UnitConstants::mm,
                                      tubePos[i][1] * Acts::UnitConstants::mm,
                                      tubePos[i][2] * Acts::UnitConstants::mm};
    driftCircles.push_back(ActsExamples::DriftCircle(
        std::move(tube_pos), driftRadius[i], 0.0f, stationInfo[i][0],
        stationInfo[i][1], stationInfo[i][2], stationInfo[i][3],
        stationInfo[i][4], stationInfo[i][5]));
  }

  // configure the binning of the hough plane
  Acts::HoughTransformUtils::HoughPlaneConfig planeCfg;
  planeCfg.nBinsX = 1000;
  planeCfg.nBinsY = 1000;

  // instantiate the peak finder
  Acts::HoughTransformUtils::PeakFinders::IslandsAroundMaxConfig peakFinderCfg;
  peakFinderCfg.fractionCutoff = 0.7;
  peakFinderCfg.threshold = 3.;
  peakFinderCfg.minSpacingBetweenPeaks = {0., 30.};

  // and map the hough plane to parameter ranges.
  // The first coordinate is tan(theta), the second is z0 in mm
  Acts::HoughTransformUtils::HoughAxisRanges axisRanges{-3., 3., -2000., 2000.};

  // create the functions parametrising the hough space lines for drift circles.
  // Note that there are two solutions for each drift circle and angle

  // left solution
  auto houghParam_fromDC_left = [](double tanTheta,
                                   const ActsExamples::DriftCircle& DC) {
    return DC.y() - tanTheta * DC.z() -
           DC.rDrift() / std::cos(std::atan(tanTheta));
  };
  // right solution
  auto houghParam_fromDC_right = [](double tanTheta,
                                    const ActsExamples::DriftCircle& DC) {
    return DC.y() - tanTheta * DC.z() +
           DC.rDrift() / std::cos(std::atan(tanTheta));
  };

  // create the function parametrising the drift radius uncertainty
  auto houghWidth_fromDC = [](double, const ActsExamples::DriftCircle& DC) {
    return std::min(DC.rDriftError() * 3.,
                    1.0);  // scale reported errors up to at least 1mm or 3
                           // times the reported error as drift circle calib not
                           // fully reliable at this stage
  };

  // store the true parameters
  // std::vector<PatternSeed> truePatterns;

  // instantiate the hough plane
  Acts::HoughTransformUtils::HoughPlane<Acts::GeometryIdentifier::Value>
      houghPlane(planeCfg);
  // also insantiate the peak finder
  Acts::HoughTransformUtils::PeakFinders::IslandsAroundMax<
      Acts::GeometryIdentifier::Value>
      peakFinder(peakFinderCfg);

  // loop over the true hits
  for (auto& sh : simHits) {
    houghPlane.reset();

    for (ActsExamples::DriftCircle& dc : driftCircles) {
      ActsExamples::muonMdtIdentifierFields idf;
      idf.multilayer = dc.multilayer();
      idf.stationEta = dc.stationEta();
      idf.stationPhi = dc.stationPhi();
      idf.stationName = dc.stationName();
      idf.tubeLayer = dc.tubeLayer();
      idf.tube = dc.tube();
      auto identifier = compressId(idf);
      auto effectiveLayer = 3 * (dc.multilayer() - 1) + (dc.tubeLayer() - 1);

      houghPlane.fill<ActsExamples::DriftCircle>(
          dc, axisRanges, houghParam_fromDC_left, houghWidth_fromDC, identifier,
          effectiveLayer, 1.0);
      houghPlane.fill<ActsExamples::DriftCircle>(
          dc, axisRanges, houghParam_fromDC_right, houghWidth_fromDC,
          identifier, effectiveLayer, 1.0);
    }

    // now get the peaks
    auto maxima = peakFinder.findPeaks(houghPlane, axisRanges);

    for (auto& max : maxima) {
      // check the Hough Transforms results
      BOOST_CHECK_CLOSE(max.x, sh.first, 5.);
      BOOST_CHECK_CLOSE(max.y, sh.second, 1.);
    }
  }
}

}  // namespace Acts::Test
