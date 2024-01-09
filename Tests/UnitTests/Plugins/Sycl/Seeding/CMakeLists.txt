# This file is part of the Acts project.
#
# Copyright (C) 2020-2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

add_executable(ActsUnitTestSeedFinderSycl
  SeedFinderSyclTest.cpp
  CommandLineArguments.h
  CommandLineArguments.cpp
  ATLASCuts.hpp
  SpacePoint.hpp)
target_link_libraries(ActsUnitTestSeedFinderSycl
  PRIVATE Boost::program_options ActsCore ActsPluginSycl
          vecmem::core vecmem::sycl)
