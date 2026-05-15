// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"

using namespace Acts;

BOOST_AUTO_TEST_SUITE(GeometryContextTests)

BOOST_AUTO_TEST_CASE(DangerouslyDefaultConstruct) {
  auto ctx = GeometryContext::dangerouslyDefaultConstruct();
  BOOST_CHECK(!ctx.hasValue());
}

BOOST_AUTO_TEST_CASE(FactoryCreatesEmptyContext) {
  auto ctx1 = GeometryContext::dangerouslyDefaultConstruct();
  auto ctx2 = GeometryContext::dangerouslyDefaultConstruct();

  // Both should be empty
  BOOST_CHECK(!ctx1.hasValue());
  BOOST_CHECK(!ctx2.hasValue());
}

BOOST_AUTO_TEST_CASE(ExplicitConstructionWithValue) {
  struct TestContext {
    int value = 42;
  };

  GeometryContext ctx{TestContext{}};
  BOOST_CHECK(ctx.hasValue());
  BOOST_CHECK_EQUAL(ctx.get<TestContext>().value, 42);
}

BOOST_AUTO_TEST_CASE(TemplateConstructorMoveSemantics) {
  struct Movable {
    int value;
    explicit Movable(int v) : value(v) {}
    Movable(const Movable&) = default;
    Movable(Movable&&) = default;
  };

  GeometryContext ctx{Movable{99}};
  BOOST_CHECK_EQUAL(ctx.get<Movable>().value, 99);
}

BOOST_AUTO_TEST_CASE(TemplateConstructorCopySemantics) {
  struct TestContext {
    int value = 123;
  };
  TestContext tc;

  GeometryContext ctx{tc};
  BOOST_CHECK_EQUAL(ctx.get<TestContext>().value, 123);
}

BOOST_AUTO_TEST_CASE(AssignmentOperatorWorks) {
  auto ctx = GeometryContext::dangerouslyDefaultConstruct();
  struct TestContext {
    double value = 3.14;
  };

  ctx = TestContext{};
  BOOST_CHECK(ctx.hasValue());
  BOOST_CHECK_EQUAL(ctx.get<TestContext>().value, 3.14);
}

BOOST_AUTO_TEST_SUITE_END()
