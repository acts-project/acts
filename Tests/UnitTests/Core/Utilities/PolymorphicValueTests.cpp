// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/PolymorphicValue.hpp"

// BOOST_AUTO_TEST_SUITE(PolymorphicValueTest)

using namespace Acts;

struct Base {
  virtual int foo() = 0;

  virtual ~Base() = default;
};

struct A : public Base {
  virtual int foo() { return 5; }

  ~A() { std::cout << "DESTROY A" << std::endl; }
};

struct B : public Base {
  virtual int foo() { return 6; }
};

BOOST_AUTO_TEST_CASE(Create) {
  {
    PolymorphicValue<Base> pm{PolymorphicValue<Base>::make<A>()};
    BOOST_CHECK_EQUAL(pm->foo(), 5);
  }

  {
    PolymorphicValue<Base> pm{PolymorphicValue<Base>::make<B>()};
    BOOST_CHECK_EQUAL(pm->foo(), 6);
  }

  {
    PolymorphicValue<Base> pm{PolymorphicValue<Base>::make<A>()};
    BOOST_CHECK_EQUAL(pm->foo(), 5);

    pm = PolymorphicValue<Base>::make<B>();
    BOOST_CHECK_EQUAL(pm->foo(), 6);
  }
}

// BOOST_AUTO_TEST_SUITE_END()