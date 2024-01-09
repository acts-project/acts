// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Plugins/Json/ExtentJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <nlohmann/json.hpp>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(ExtentJsonConverter)

BOOST_AUTO_TEST_CASE(ExtentRoundtripTests) {
  Extent e;
  e.set(binR, 0, 200);
  e.set(binZ, -50, 50);

  nlohmann::json j;
  j["extent"] = e;

  Extent eIn = j["extent"];

  CHECK_CLOSE_ABS(eIn.min(binR), e.min(binR), 10e-5);
  CHECK_CLOSE_ABS(eIn.max(binR), e.max(binR), 10e-5);
  CHECK_CLOSE_ABS(eIn.min(binZ), e.min(binZ), 10e-5);
  CHECK_CLOSE_ABS(eIn.max(binZ), e.max(binZ), 10e-5);
}

BOOST_AUTO_TEST_SUITE_END()
