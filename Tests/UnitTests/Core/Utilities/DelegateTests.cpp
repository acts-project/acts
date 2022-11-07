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

int sumImpl(int a, int b) {
  return a + b;
}

BOOST_AUTO_TEST_CASE(ConnectConstexprLambda) {
  Delegate<int(int, int)> sum;
  BOOST_CHECK(!sum);
  BOOST_CHECK(!sum.connected());

  sum.connect<&sumImpl>();

  BOOST_CHECK_EQUAL(sum(2, 5), 7);
  BOOST_CHECK_NE(sum(2, 3), 7);

  sum.connect([](const void*, int a, int b) -> int { return a + b; });

  BOOST_CHECK(sum);
  BOOST_CHECK(sum.connected());

  BOOST_CHECK_EQUAL(sum(2, 5), 7);
  BOOST_CHECK_NE(sum(2, 3), 7);
}

float multiply(float a, float b) {
  return a * b;
}

BOOST_AUTO_TEST_CASE(ConnectFunctionPointer) {
  Delegate<float(float, float)> mult;

  BOOST_CHECK(!mult);
  BOOST_CHECK(!mult.connected());

  mult.connect<multiply>();

  BOOST_CHECK(mult);
  BOOST_CHECK(mult.connected());

  CHECK_CLOSE_REL(mult(2, 5.9), 2 * 5.9, 1e-6);
  BOOST_CHECK_NE(mult(2, 3.2), 58.9);
}

struct Subtractor {
  int v;
  int execute(int a) const { return a - v; }
};

BOOST_AUTO_TEST_CASE(ConnectStruct) {
  Delegate<int(int)> sub;

  BOOST_CHECK(!sub);
  BOOST_CHECK(!sub.connected());

  Subtractor s{18};
  sub.connect<&Subtractor::execute>(&s);

  BOOST_CHECK(sub);
  BOOST_CHECK(sub.connected());

  BOOST_CHECK_EQUAL(sub(7), 7 - 18);
}

int addition(const void*, int a, int b) {
  return a + b;
}

BOOST_AUTO_TEST_CASE(ConnectRuntime) {
  {
    Delegate<int(int, int)> add;
    BOOST_CHECK(!add);
    BOOST_CHECK(!add.connected());

    add.connect(&addition);
    BOOST_CHECK(add);
    BOOST_CHECK(add.connected());

    BOOST_CHECK_EQUAL(add(4, 4), 8);
  }

  {
    Delegate<int(int, int)> add{&addition};

    BOOST_CHECK(add);
    BOOST_CHECK(add.connected());

    BOOST_CHECK_EQUAL(add(4, 4), 8);
  }

  {
    Delegate<int(int, int)> add;
    BOOST_CHECK(!add);
    BOOST_CHECK(!add.connected());

    add = &addition;
    BOOST_CHECK(add);
    BOOST_CHECK(add.connected());

    BOOST_CHECK_EQUAL(add(4, 4), 8);
  }
}

void modify(int& v, int a) {
  v = a;
}

void noModify(int v, int a) {
  (void)v;
  v = a;
}

BOOST_AUTO_TEST_CASE(DelegateReference) {
  Delegate<void(int&, int)> d;
  d.connect<&modify>();

  int v = 0;
  d(v, 42);
  BOOST_CHECK_EQUAL(v, 42);

  // This should not compile since the signature is not exactly matched
  // d.connect<&noModify>();
}

struct SignatureTest {
  void modify(int& v, int a) const { v = a; }

  void noModify(int v, int a) const {
    (void)v;
    v = a;
  }
};

BOOST_AUTO_TEST_CASE(DelegateReferenceMember) {
  SignatureTest s;
  Delegate<void(int&, int)> d;
  d.connect<&SignatureTest::modify>(&s);

  int v = 0;
  d(v, 42);
  BOOST_CHECK_EQUAL(v, 42);

  // This should not compile since the signature is not exactly matched
  // d.connect<&SignatureTest::noModify>(&s);
}

BOOST_AUTO_TEST_CASE(StatefullLambdas) {
  std::vector<int> v;

  auto lambda = [&](int n) -> int {
    v.push_back(n);
    return v.size();
  };

  Delegate<int(int)> d(lambda);

  BOOST_CHECK(d);
  BOOST_CHECK(d.connected());
  BOOST_CHECK(d(2) == 1);

  d.disconnect();
  d = lambda;

  BOOST_CHECK(d);
  BOOST_CHECK(d.connected());
  BOOST_CHECK(d(2) == 2);

  d.disconnect();
  d.connect(lambda);

  BOOST_CHECK(d);
  BOOST_CHECK(d.connected());
  BOOST_CHECK(d(2) == 3);

  // This should not compile because of deleted && overloads
  // d.connect([&](int a){ v.push_back(a); return v.size(); });
}

BOOST_AUTO_TEST_SUITE_END()
