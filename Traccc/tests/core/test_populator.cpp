// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "traccc/seeding/grids/populator.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace traccc;

GTEST_TEST(traccc_grid2, replace_populator) {
  replace_populator<> replacer;
  unsigned int stored = 3u;
  replacer(stored, 2u);
  EXPECT_EQ(stored, 2u);

  replacer(stored, 42u);
  EXPECT_EQ(stored, 42u);
}

GTEST_TEST(traccc_grid2, complete_populator) {
  using cpopulator4 = complete_populator<vecmem::vector, vecmem::jagged_vector,
                                         std::array, unsigned int, false, 4u>;
  cpopulator4 completer;

  cpopulator4::store_value stored = {completer.m_invalid, completer.m_invalid,
                                     completer.m_invalid, completer.m_invalid};

  cpopulator4::store_value test = stored;
  test[0] = 9u;
  completer(stored, 9u);
  EXPECT_EQ(stored, test);

  test[1] = 3u;
  completer(stored, 3u);
  EXPECT_EQ(stored, test);

  using sort_cpopulator4 =
      complete_populator<vecmem::vector, vecmem::jagged_vector, std::array,
                         unsigned int, true, 4u>;
  sort_cpopulator4 sort_completer;

  test = {0u, 3u, 9u, 1000u};
  sort_completer(stored, 1000u);
  sort_completer(stored, 0u);
  EXPECT_EQ(stored, test);
}

GTEST_TEST(traccc_grid2, attach_populator) {
  // Attach populator without sorting
  attach_populator<> attacher;
  attach_populator<>::store_value stored = {3u};
  attacher(stored, 2u);
  attach_populator<>::store_value test = {3u, 2u};
  EXPECT_EQ(stored, test);

  attacher(stored, 42u);
  test = {3u, 2u, 42u};
  EXPECT_EQ(stored, test);

  // Attach populator with sorting
  attach_populator<vecmem::vector, vecmem::jagged_vector, std::array,
                   unsigned int, true>
      sort_attacher;
  sort_attacher(stored, 11u);
  test = {2u, 3u, 11u, 42u};
  EXPECT_EQ(stored, test);
}
