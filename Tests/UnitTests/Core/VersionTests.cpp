// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// the purpose of these tests is to ensure that the version header is valid and
// the exported parameters are accessible. otherwise, there would is no code
// that actually includes the version header.

#include <boost/test/unit_test.hpp>
#include <string_view>

#include "Acts/ActsVersion.hpp"

BOOST_AUTO_TEST_CASE(Version) {
  // the only way to get a zero version would be zero for all components
  BOOST_TEST(0u < Acts::Version);
  // these tests are not really useful as the version components can be any
  // value. they are there so we touch all variables and ensure that they are
  // accessible.
  BOOST_TEST(0u <= Acts::VersionMajor);
  BOOST_TEST(0u <= Acts::VersionMinor);
  BOOST_TEST(0u <= Acts::VersionPatch);
}

BOOST_AUTO_TEST_CASE(CommitHash) {
  BOOST_TEST(not std::string_view(Acts::CommitHash).empty());
  BOOST_TEST(not std::string_view(Acts::CommitHashShort).empty());
}
