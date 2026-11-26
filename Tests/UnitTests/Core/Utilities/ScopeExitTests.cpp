// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/ScopeExit.hpp"

#include <stdexcept>

using Acts::detail::ScopeExit;

BOOST_AUTO_TEST_SUITE(ScopeExitTests)

BOOST_AUTO_TEST_CASE(CalledOnceAtScopeEnd) {
  int calls = 0;
  {
    const ScopeExit guard([&] { ++calls; });
    BOOST_CHECK_EQUAL(calls, 0);
  }
  BOOST_CHECK_EQUAL(calls, 1);
}

BOOST_AUTO_TEST_CASE(CalledOnEarlyReturn) {
  int calls = 0;
  auto f = [&](bool early) {
    const ScopeExit guard([&] { ++calls; });
    if (early) {
      return 1;
    }
    return 2;
  };
  BOOST_CHECK_EQUAL(f(true), 1);
  BOOST_CHECK_EQUAL(f(false), 2);
  BOOST_CHECK_EQUAL(calls, 2);
}

BOOST_AUTO_TEST_CASE(CalledOnException) {
  int calls = 0;
  BOOST_CHECK_THROW(
      {
        const ScopeExit guard([&] { ++calls; });
        throw std::runtime_error("scope");
      },
      std::runtime_error);
  BOOST_CHECK_EQUAL(calls, 1);
}

BOOST_AUTO_TEST_SUITE_END()
