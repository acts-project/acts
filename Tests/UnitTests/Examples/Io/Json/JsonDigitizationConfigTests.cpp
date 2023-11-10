// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(JsonDigitizationConfig)

BOOST_AUTO_TEST_CASE(DigitizationConfigRoundTrip) {
  std::ofstream out;

  // As all SurfaceBounds have the same streaming API only a one is
  // tested here, all others are tests are identical

  ActsExamples::DigiComponentsConfig dcf;

  ActsExamples::GeometricConfig gdc;

  Acts::BinUtility segmentation;
  segmentation += Acts::BinUtility(336, -8.4, 8.4, Acts::open, Acts::binX);
  segmentation += Acts::BinUtility(1280, -36, 36, Acts::open, Acts::binY);

  gdc.segmentation = segmentation;
  gdc.threshold = 0.01;
  gdc.thickness = 0.15;
  gdc.indices = {Acts::eBoundLoc0, Acts::eBoundLoc1};
  gdc.chargeSmearer = ActsExamples::Digitization::Gauss(1.0);

  ActsExamples::DigiComponentsConfig dcRef;
  dcRef.geometricDigiConfig = gdc;

  nlohmann::json dcJsonOut(dcRef);
  out.open("DigiComponentsConfig.json");
  out << dcJsonOut.dump(2);
  out.close();

  auto in = std::ifstream("DigiComponentsConfig.json",
                          std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json dcJsonIn;
  in >> dcJsonIn;
  in.close();

  ActsExamples::DigiComponentsConfig dcTest(dcJsonIn);
  BOOST_CHECK(dcTest.geometricDigiConfig.indices ==
              dcRef.geometricDigiConfig.indices);
  BOOST_CHECK_EQUAL(dcTest.geometricDigiConfig.segmentation.dimensions(),
                    dcRef.geometricDigiConfig.segmentation.dimensions());
}

BOOST_AUTO_TEST_SUITE_END()
