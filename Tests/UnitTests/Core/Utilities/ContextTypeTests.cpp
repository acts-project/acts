// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/ContextType.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(ContextTypeTests)

BOOST_AUTO_TEST_CASE(PackUnpack) {
  ContextType ctx;

  int v = 42;
  ctx = v;

  BOOST_CHECK_EQUAL(ctx.get<int>(), 42);
  BOOST_CHECK_THROW(ctx.get<double>(), std::bad_any_cast);
}

BOOST_AUTO_TEST_CASE(MaybeUnpack) {
  ContextType ctx;

  int v = 42;
  ctx = v;

  BOOST_CHECK_EQUAL(*ctx.maybeGet<int>(), 42);
  BOOST_CHECK_EQUAL(ctx.maybeGet<double>(), nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
