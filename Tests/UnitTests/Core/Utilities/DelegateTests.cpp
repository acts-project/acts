// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <numeric>
#include <optional>
#include <random>
#include <tuple>

using namespace Acts;

namespace bd = boost::unit_test::data;

BOOST_AUTO_TEST_SUITE(DelegateTests)

BOOST_AUTO_TEST_CASE(ConnectConstexprLambda) {
  Delegate<int(int, int)> sum;

  constexpr int (*f)(int, int) = [](int a, int b) { return a + b; };

  sum.connect<f>();

  BOOST_CHECK_EQUAL(sum(2, 5), 7);
  BOOST_CHECK_NE(sum(2, 3), 7);
}

float multiply(float a, float b) {
  return a * b;
}

BOOST_AUTO_TEST_CASE(ConnectFunctionPointer) {
  Delegate<float(float, float)> mult;

  mult.connect<multiply>();

  CHECK_CLOSE_REL(mult(2, 5.9), 2 * 5.9, 1e-6);
  BOOST_CHECK_NE(mult(2, 3.2), 58.9);
}

struct Subtractor {
  int v;
  int execute(int a) const { return a - v; }
};

BOOST_AUTO_TEST_CASE(ConnectStruct) {
  Delegate<int(int)> sub;

  Subtractor s{18};
  sub.connect<&Subtractor::execute>(&s);

  BOOST_CHECK_EQUAL(sub(7), 7 - 18);
}

BOOST_AUTO_TEST_SUITE_END()
