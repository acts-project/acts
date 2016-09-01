// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE Unit Conversion Tests
#include <boost/test/included/unit_test.hpp>
#include "ACTS/Utilities/Units.hpp"

using namespace Acts::units;
namespace utf = boost::unit_test;
namespace tt  = boost::test_tools;

BOOST_AUTO_TEST_SUITE(unit_conversion, *utf::tolerance(1e-15));

BOOST_AUTO_TEST_CASE(length_conversions)
{
  BOOST_TEST(_m == 1e-3 * _km);
  BOOST_TEST(_m == 1e3 * _mm);
  BOOST_TEST(_m == 1e6 * _um);
  BOOST_TEST(_m == 1e9 * _nm);
  BOOST_TEST(_m == 1e12 * _pm);
  BOOST_TEST(_m == 1e15 * _fm);
}

BOOST_AUTO_TEST_CASE(mass_conversions)
{
  BOOST_TEST(_kg == 1e3 * _g);
  BOOST_TEST(_kg == 1e6 * _mg);
  BOOST_TEST(1.660539040e-27 * _kg == _u);
}

BOOST_AUTO_TEST_CASE(time_conversions)
{
  BOOST_TEST(_s == 1e3 * _ms);
  BOOST_TEST(3600 * _s == _h);
}

BOOST_AUTO_TEST_CASE(energy_conversions)
{
  BOOST_TEST(_MeV == 1e-3 * _GeV);
  BOOST_TEST(_MeV == 1e-6 * _TeV);
  BOOST_TEST(_MeV == 1e3 * _keV);
  BOOST_TEST(_MeV == 1e6 * _eV);
  BOOST_TEST(_MeV == 1.60217733e-13 * _J);
}

BOOST_AUTO_TEST_SUITE_END();
