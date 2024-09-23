// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Utilities/DelegateChain.hpp"

using namespace Acts;

struct AddTo {
  int value = 0;

  void add(int &x) const { x += value; }
};

BOOST_AUTO_TEST_SUITE(DelegateChainTests)

BOOST_AUTO_TEST_CASE(DelegateChainAdd) {
  DelegateChain<void(int &)> chain;
  AddTo a1{1}, a2{2}, a3{3};

  chain.connect_back<&AddTo::add>(&a1);
  chain.connect_back<&AddTo::add>(&a2);
  chain.connect_back<&AddTo::add>(&a3);

  int x = 0;
  chain(x);
  BOOST_CHECK_EQUAL(x, 6);
}

BOOST_AUTO_TEST_SUITE_END()
