// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Unit Conversion Tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::units;
namespace utf = boost::unit_test;
namespace tt  = boost::test_tools;

BOOST_AUTO_TEST_SUITE(unit_conversion)

BOOST_AUTO_TEST_CASE(length_conversions)
{
  CHECK_CLOSE_REL(_m, 1e-3 * _km, 1e-15);
  CHECK_CLOSE_REL(_m, 1e3 * _mm, 1e-15);
  CHECK_CLOSE_REL(_m, 1e6 * _um, 1e-15);
  CHECK_CLOSE_REL(_m, 1e9 * _nm, 1e-15);
  CHECK_CLOSE_REL(_m, 1e12 * _pm, 1e-15);
  CHECK_CLOSE_REL(_m, 1e15 * _fm, 1e-15);
}

BOOST_AUTO_TEST_CASE(mass_conversions)
{
  CHECK_CLOSE_REL(_kg, 1e3 * _g, 1e-15);
  CHECK_CLOSE_REL(_kg, 1e6 * _mg, 1e-15);
  CHECK_CLOSE_REL(1.660539040e-27 * _kg, _u, 1e-15);
}

BOOST_AUTO_TEST_CASE(time_conversions)
{
  CHECK_CLOSE_REL(_s, 1e3 * _ms, 1e-15);
  CHECK_CLOSE_REL(3600 * _s, _h, 1e-15);
}

BOOST_AUTO_TEST_CASE(energy_conversions)
{
  CHECK_CLOSE_REL(_MeV, 1e-3 * _GeV, 1e-15);
  CHECK_CLOSE_REL(_MeV, 1e-6 * _TeV, 1e-15);
  CHECK_CLOSE_REL(_MeV, 1e3 * _keV, 1e-15);
  CHECK_CLOSE_REL(_MeV, 1e6 * _eV, 1e-15);
}

// clang-format off
BOOST_AUTO_TEST_CASE(si_nat_conversions)
{
  CHECK_CLOSE_REL(SI2Nat<ENERGY>(1 * _J), 1. / 1.60217733e-19 * _eV, 1e-8);
  CHECK_CLOSE_REL(SI2Nat<ENERGY>(1 * _g * _m * _m / (_s * _s)),
                  1. / 1.60217733e-16 * _eV,
                  1e-8);
  CHECK_CLOSE_REL(SI2Nat<ENERGY>(1 * _J), 1. / 1.60217733e-10 * _GeV, 1e-8);
  CHECK_CLOSE_REL(Nat2SI<ENERGY>(100 * _keV), 1.60217733e-14 * _J, 1e-8);
  CHECK_CLOSE_REL(Nat2SI<ENERGY>(SI2Nat<ENERGY>(100 * _J)), 100 * _J, 1e-8);
  CHECK_CLOSE_REL(SI2Nat<ENERGY>(Nat2SI<ENERGY>(57.3 * _MeV)),
                  57300 * _keV,
                  1e-8);

  CHECK_CLOSE_REL(SI2Nat<LENGTH>(1 * _m), 5.0677289e15 / _GeV, 1e-8);
  CHECK_CLOSE_REL(Nat2SI<LENGTH>(1. / (197.3270523 * _MeV)), 1 * _fm, 1e-8);
  CHECK_CLOSE_REL(Nat2SI<LENGTH>(SI2Nat<LENGTH>(10 * _m)), 1e-2 * _km, 1e-8);
  CHECK_CLOSE_REL(SI2Nat<LENGTH>(Nat2SI<LENGTH>(1. / _keV)),
                  1 / (1e-3 * _MeV),
                  1e-8);

  CHECK_CLOSE_REL(SI2Nat<MOMENTUM>(2 * _g * _m / _s),
                  3.742313068429198e15 * _GeV,
                  1e-8);
  CHECK_CLOSE_REL(Nat2SI<MOMENTUM>(13 * _TeV),
                  6.947574808569e-15 * _kg * _m / _s,
                  1e-8);
  CHECK_CLOSE_REL(Nat2SI<MOMENTUM>(SI2Nat<MOMENTUM>(3.6 * _mg * _km / _h)),
                  1e-6 * _kg * _m / _s,
                  1e-8);
  CHECK_CLOSE_REL(SI2Nat<MOMENTUM>(Nat2SI<MOMENTUM>(6.5 * _GeV)),
                  6.5e6 * _keV,
                  1e-8);

  CHECK_CLOSE_REL(SI2Nat<MASS>(1.67262263806e-27 * _kg),
                  938.2720813 * _MeV,
                  1e-8);
  CHECK_CLOSE_REL(Nat2SI<MASS>(125.09 * _GeV), 2.229932766466e-19 * _mg, 1e-8);
  CHECK_CLOSE_REL(Nat2SI<MASS>(SI2Nat<MASS>(3 * _mg)), 3e-3 * _g, 1e-8);
  CHECK_CLOSE_REL(SI2Nat<MASS>(Nat2SI<MASS>(1.77 * _GeV)), 1770 * _MeV, 1e-8);

  CHECK_CLOSE_REL(SI2Nat<ENERGY>(_hbar * _c),
                  0.197327052356 * _fm * _GeV,
                  1e-8);

  CHECK_CLOSE_REL(Nat2SI<MOMENTUM>(_GeV) / (_el_charge * _C * 2 * _T),
              10./6. * _m,
              1e-3);
}
// clang-format on

BOOST_AUTO_TEST_SUITE_END()
