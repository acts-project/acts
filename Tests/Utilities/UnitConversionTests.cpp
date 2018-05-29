// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE Unit Conversion Tests
#include <boost/test/included/unit_test.hpp>
#include "Acts/Utilities/Units.hpp"

using namespace Acts::units;
namespace utf = boost::unit_test;
namespace tt  = boost::test_tools;

BOOST_AUTO_TEST_SUITE(unit_conversion, *utf::tolerance(1e-15))

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
}

// clang-format off
BOOST_AUTO_TEST_CASE(si_nat_conversions, *utf::tolerance(1e-8))
{
  BOOST_TEST(SI2Nat<ENERGY>(1 * _J) == 1. / 1.60217733e-19 * _eV);
  BOOST_TEST(SI2Nat<ENERGY>(1 * _g * _m * _m / (_s * _s)) == 1. / 1.60217733e-16 * _eV);
  BOOST_TEST(SI2Nat<ENERGY>(1 * _J) == 1. / 1.60217733e-10 * _GeV);
  BOOST_TEST(Nat2SI<ENERGY>(100 * _keV) == 1.60217733e-14 * _J);
  BOOST_TEST(Nat2SI<ENERGY>(SI2Nat<ENERGY>(100 * _J)) == 100 * _J);
  BOOST_TEST(SI2Nat<ENERGY>(Nat2SI<ENERGY>(57.3 * _MeV)) == 57300 * _keV);

  BOOST_TEST(SI2Nat<LENGTH>(1 * _m) == 5.0677289e15 / _GeV);
  BOOST_TEST(Nat2SI<LENGTH>(1. / (197.3270523 * _MeV)) == 1 * _fm);
  BOOST_TEST(Nat2SI<LENGTH>(SI2Nat<LENGTH>(10 * _m)) == 1e-2 * _km);
  BOOST_TEST(SI2Nat<LENGTH>(Nat2SI<LENGTH>(1. / _keV)) == 1 / (1e-3 * _MeV));

  BOOST_TEST(SI2Nat<MOMENTUM>(2 * _g * _m / _s) == 3.742313068429198e15 * _GeV);
  BOOST_TEST(Nat2SI<MOMENTUM>(13 * _TeV) == 6.947574808569e-15 * _kg * _m / _s);
  BOOST_TEST(Nat2SI<MOMENTUM>(SI2Nat<MOMENTUM>(3.6 * _mg * _km / _h)) == 1e-6 * _kg * _m / _s);
  BOOST_TEST(SI2Nat<MOMENTUM>(Nat2SI<MOMENTUM>(6.5 * _GeV)) == 6.5e6 * _keV);

  BOOST_TEST(SI2Nat<MASS>(1.67262263806e-27 * _kg) == 938.2720813 * _MeV);
  BOOST_TEST(Nat2SI<MASS>(125.09 * _GeV) == 2.229932766466e-19 * _mg);
  BOOST_TEST(Nat2SI<MASS>(SI2Nat<MASS>(3 * _mg)) == 3e-3 * _g);
  BOOST_TEST(SI2Nat<MASS>(Nat2SI<MASS>(1.77 * _GeV)) == 1770 * _MeV);

  BOOST_TEST(SI2Nat<ENERGY>(_hbar * _c) == 0.197327052356 * _fm * _GeV);

  BOOST_TEST(Nat2SI<MOMENTUM>(_GeV) / (_el_charge * _C * 2 * _T) == 10./6. * _m, tt::tolerance(1e-3));
}
// clang-format on

BOOST_AUTO_TEST_SUITE_END()
