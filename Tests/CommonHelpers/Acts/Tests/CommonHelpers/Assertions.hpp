// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/unit_test.hpp>

#include <vector>

#define CHECK_NE_COLLECTIONS(col1, col2)                                       \
  do {                                                                         \
    BOOST_CHECK_EQUAL(col1.size(), col2.size());                               \
    std::vector<bool> result;                                                  \
    for (std::size_t i = 0; i < col1.size(); i++) {                            \
      result.push_back(col1[i] == col2[i]);                                    \
    }                                                                          \
    BOOST_CHECK(                                                               \
        !std::all_of(result.begin(), result.end(), [](bool r) { return r; })); \
  } while (0)
