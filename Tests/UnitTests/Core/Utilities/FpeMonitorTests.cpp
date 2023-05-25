// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/FpeMonitor.hpp"

#include <optional>

namespace utf = boost::unit_test;

namespace {

void divbyzero() {
  volatile float j = 0.0;
  volatile float r = 123 / j;
  (void)r;
}

void overflow() {
  std::cout << "PRE OVERFLOW" << std::endl;
  volatile float j = std::numeric_limits<float>::max();
  volatile float r = j * j;
  (void)r;
  std::cout << "POST OVERFLOW" << std::endl;
}

void invalid() {
  volatile float j = -1;
  volatile float r = std::sqrt(j);
  (void)r;
}

}  // namespace

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(FpeMonitorTest)

BOOST_AUTO_TEST_CASE(Invalid) {
  {
    FpeMonitor mon;
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    invalid();

    BOOST_CHECK(mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));
  }
}

BOOST_AUTO_TEST_CASE(DivByZero) {
  {
    FpeMonitor mon;
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    divbyzero();

    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(mon.encountered(FpeType::FLTDIV));
  }
}

BOOST_AUTO_TEST_CASE(Overflow) {
  {
    FpeMonitor mon;
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    overflow();
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));
  }
}

BOOST_AUTO_TEST_CASE(Combinations) {
  {
    FpeMonitor mon;
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    invalid();
    BOOST_CHECK(mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    overflow();
    BOOST_CHECK(mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    divbyzero();

    BOOST_CHECK(mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(mon.encountered(FpeType::FLTDIV));
  }
}

BOOST_AUTO_TEST_CASE(ClearOnEnter) {
  invalid();
  divbyzero();
  overflow();

  {
    FpeMonitor mon;
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));
  }
}

BOOST_AUTO_TEST_CASE(Scoping) {
  {
    FpeMonitor mon;
    BOOST_CHECK(!mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));

    invalid();

    {
      FpeMonitor mon2;
      BOOST_CHECK(!mon2.encountered(FpeType::FLTINV));
      BOOST_CHECK(!mon2.encountered(FpeType::FLTOVF));
      BOOST_CHECK(!mon2.encountered(FpeType::FLTDIV));

      overflow();

      {
        FpeMonitor mon3;
        BOOST_CHECK(!mon3.encountered(FpeType::FLTINV));
        BOOST_CHECK(!mon3.encountered(FpeType::FLTOVF));
        BOOST_CHECK(!mon3.encountered(FpeType::FLTDIV));

        divbyzero();

        BOOST_CHECK(!mon3.encountered(FpeType::FLTINV));
        BOOST_CHECK(!mon3.encountered(FpeType::FLTOVF));
        BOOST_CHECK(mon3.encountered(FpeType::FLTDIV));
      }

      BOOST_CHECK(!mon2.encountered(FpeType::FLTINV));
      BOOST_CHECK(mon2.encountered(FpeType::FLTOVF));
      BOOST_CHECK(!mon2.encountered(FpeType::FLTDIV));
    }

    BOOST_CHECK(mon.encountered(FpeType::FLTINV));
    BOOST_CHECK(!mon.encountered(FpeType::FLTOVF));
    BOOST_CHECK(!mon.encountered(FpeType::FLTDIV));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
