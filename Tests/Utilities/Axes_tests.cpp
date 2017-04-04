// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE grid axis tests
#include <boost/test/included/unit_test.hpp>

#include "ACTS/Utilities/detail/Axis.hpp"

namespace Acts {

using namespace detail;

namespace Test {

  BOOST_AUTO_TEST_CASE(equidistant_axis)
  {
    EquidistantAxis a(0.0, 10.0, 10u);

    BOOST_TEST(a.getNBins() == 10u);
    BOOST_TEST(a.getMax() == 10.);
    BOOST_TEST(a.getMin() == 0.);
    BOOST_TEST(a.getBinWidth() == 1.);
    BOOST_TEST(a.getBin(-0.3) == 0u);
    BOOST_TEST(a.getBin(-0.) == 1u);
    BOOST_TEST(a.getBin(0.) == 1u);
    BOOST_TEST(a.getBin(0.7) == 1u);
    BOOST_TEST(a.getBin(1) == 2u);
    BOOST_TEST(a.getBin(1.2) == 2u);
    BOOST_TEST(a.getBin(2.) == 3u);
    BOOST_TEST(a.getBin(2.7) == 3u);
    BOOST_TEST(a.getBin(3.) == 4u);
    BOOST_TEST(a.getBin(3.6) == 4u);
    BOOST_TEST(a.getBin(4.) == 5u);
    BOOST_TEST(a.getBin(4.98) == 5u);
    BOOST_TEST(a.getBin(5.) == 6u);
    BOOST_TEST(a.getBin(5.12) == 6u);
    BOOST_TEST(a.getBin(6.) == 7u);
    BOOST_TEST(a.getBin(6.00001) == 7u);
    BOOST_TEST(a.getBin(7.) == 8u);
    BOOST_TEST(a.getBin(7.5) == 8u);
    BOOST_TEST(a.getBin(8.) == 9u);
    BOOST_TEST(a.getBin(8.1) == 9u);
    BOOST_TEST(a.getBin(9.) == 10u);
    BOOST_TEST(a.getBin(9.999) == 10u);
    BOOST_TEST(a.getBin(10.) == 11u);
    BOOST_TEST(a.getBin(100.3) == 11u);
  }

  BOOST_AUTO_TEST_CASE(variable_axis)
  {
    VariableAxis a({0, 0.5, 3, 4.5, 6});

    BOOST_TEST(a.getNBins() == 4u);
    BOOST_TEST(a.getMax() == 6.);
    BOOST_TEST(a.getMin() == 0.);
    BOOST_TEST(a.getBin(-0.3) == 0u);
    BOOST_TEST(a.getBin(-0.) == 1u);
    BOOST_TEST(a.getBin(0.) == 1u);
    BOOST_TEST(a.getBin(0.3) == 1u);
    BOOST_TEST(a.getBin(0.5) == 2u);
    BOOST_TEST(a.getBin(1.2) == 2u);
    BOOST_TEST(a.getBin(2.7) == 2u);
    BOOST_TEST(a.getBin(3.) == 3u);
    BOOST_TEST(a.getBin(4.49999) == 3u);
    BOOST_TEST(a.getBin(4.5) == 4u);
    BOOST_TEST(a.getBin(5.12) == 4u);
    BOOST_TEST(a.getBin(6.) == 5u);
    BOOST_TEST(a.getBin(6.00001) == 5u);
    BOOST_TEST(a.getBin(7.5) == 5u);
  }
}  // namespace Test

}  // namespace Acts
