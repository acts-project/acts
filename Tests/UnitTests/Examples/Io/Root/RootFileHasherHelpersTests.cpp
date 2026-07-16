// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Io/Root/detail/RootFileHasherHelpers.hpp"

#include <string>

using ActsExamples::detail::firstTemplateArg;

BOOST_AUTO_TEST_SUITE(RootFileHasherHelpers)

BOOST_AUTO_TEST_CASE(FirstTemplateArgSimple) {
  BOOST_CHECK_EQUAL(firstTemplateArg("vector<float>"), "float");
  BOOST_CHECK_EQUAL(firstTemplateArg("vector<unsigned int>"), "unsigned int");
}

BOOST_AUTO_TEST_CASE(FirstTemplateArgNested) {
  // Legacy ROOT spelling with a space before the closing angle bracket.
  BOOST_CHECK_EQUAL(firstTemplateArg("vector<vector<int> >"), "vector<int>");
  // Modern spelling without the space.
  BOOST_CHECK_EQUAL(firstTemplateArg("vector<vector<int>>"), "vector<int>");
}

BOOST_AUTO_TEST_CASE(FirstTemplateArgStopsAtTopLevelComma) {
  // Only the first top-level argument is returned; a nested comma is ignored.
  BOOST_CHECK_EQUAL(firstTemplateArg("map<int,float>"), "int");
  BOOST_CHECK_EQUAL(firstTemplateArg("map<vector<int,short>,float>"),
                    "vector<int,short>");
}

BOOST_AUTO_TEST_CASE(FirstTemplateArgTrimsWhitespace) {
  BOOST_CHECK_EQUAL(firstTemplateArg("vector< float >"), "float");
}

BOOST_AUTO_TEST_CASE(FirstTemplateArgNoTemplate) {
  BOOST_CHECK_EQUAL(firstTemplateArg("Float_t"), "");
  BOOST_CHECK_EQUAL(firstTemplateArg(""), "");
}

BOOST_AUTO_TEST_SUITE_END()
