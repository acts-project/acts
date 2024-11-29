// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/context.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Utilities/CheckedCast.hpp"

using namespace Acts;

BOOST_AUTO_TEST_SUITE(CheckedCastTest)

struct Base {
  virtual ~Base() = default;
};
struct A : public Base {
  int value = 0;
};

struct B : public Base {};

BOOST_AUTO_TEST_CASE(ValidDowncast) {
  A a;
  a.value = 5;

  Base* aBase = &a;

  A* aDowncast = checked_cast<A*>(aBase);

  BOOST_CHECK_EQUAL(aDowncast->value, 5);

  A& aRef = checked_cast<A&>(*aBase);
  BOOST_CHECK_EQUAL(aRef.value, 5);
}

BOOST_AUTO_TEST_CASE(InvalidDowncast) {
  B b;
  Base* bBase = &b;
  BOOST_CHECK_EQUAL(checked_cast<A*>(bBase), nullptr);
  BOOST_CHECK_THROW(checked_cast<A&>(*bBase), std::bad_cast);
}

BOOST_AUTO_TEST_SUITE_END()
