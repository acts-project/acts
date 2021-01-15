// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>

BOOST_AUTO_TEST_SUITE(UtilitiesJsonConverter)

BOOST_AUTO_TEST_CASE(BinUtilityRoundTripTests) {

  Acts::BinUtility reference(Acts::binR, 0, 4);

  std::ofstream out;


  // Test in in one dimension
  nlohmann::json joneDimOut;
  to_json(joneDimOut, reference);
  out.open("oneDim.json");
  out << joneDimOut.dump(2);
  out.close();

  auto in = std::ifstream("oneDim.json", std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json joneDimIn;
  in >> joneDimIn;
  in.close();

  Acts::BinUtility test;
  from_json(joneDimIn, test);

  // Increase to two dimensions
  reference += Acts::BinUtility(10., -M_PI, M_PI, Acts::closed, Acts::binPhi);
  nlohmann::json jtwoDimOut;
  to_json(jtwoDimOut, reference);
  out.open("twoDim.json");
  out << jtwoDimOut.dump(2);
  out.close();

  in = std::ifstream("twoDim.json", std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json jtwoDimIn;
  in >> jtwoDimIn;
  in.close();
  
  test = Acts::BinUtility();
  from_json(jtwoDimIn, test);


}

BOOST_AUTO_TEST_SUITE_END()
