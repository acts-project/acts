// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"

#include <format>
#include <iostream>

#include <Eigen/Core>
#include <boost/version.hpp>

int main(void) {
  std::cout << std::format("Using Acts version {}.{}.{} commit {}\n",
                           Acts::VersionMajor, Acts::VersionMinor,
                           Acts::VersionPatch,
                           Acts::CommitHash.value_or("UNKNOWN"))
            << std::endl;
  std::cout << std::format("Using Boost version {}.{}.{}\n",
                           BOOST_VERSION / 100000, BOOST_VERSION / 100 % 1000,
                           BOOST_VERSION % 100)
            << std::endl;
  std::cout << std::format("Using Eigen version {}.{}.{}\n",
                           EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION,
                           EIGEN_MINOR_VERSION)
            << std::endl;
  return EXIT_SUCCESS;
}
