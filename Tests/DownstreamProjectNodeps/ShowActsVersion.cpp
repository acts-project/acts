// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/ActsVersion.hpp>

#include <cstdio>
#include <cstdlib>

#include <Eigen/Core>
#include <boost/version.hpp>

int main(void) {
  printf("Using Acts version %u.%u.%u commit %s\n", Acts::VersionMajor,
         Acts::VersionMinor, Acts::VersionPatch, Acts::CommitHash);
  printf("Using Boost version %u.%u.%u\n", BOOST_VERSION / 100000,
         BOOST_VERSION / 100 % 1000, BOOST_VERSION % 100);
  printf("Using Eigen version %u.%u.%u\n", EIGEN_WORLD_VERSION,
         EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
  return EXIT_SUCCESS;
}
