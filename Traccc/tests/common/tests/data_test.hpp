// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/io/utils.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <string>

namespace traccc::tests {
/**
 * @brief Test fixture base class for accessing the data directory.
 *
 * We have a lot of tests which need to access data from the data directory,
 * which is set via the TRACCC_TEST_DATA_DIR environment variable. This raises
 * a lot of boilerplate code though, so we can abstract some of this away by
 * defining a fixture base class which automatically grabs the data directory.
 *
 * @warning If you need to write a fixture which subclasses this class and you
 * need a custom SetUp function, make sure to mark it as `override` and make
 * sure to call `data_test::SetUp()` at the start of it.
 */
class data_test : public ::testing::Test {
 protected:
  std::string data_directory;

  virtual void SetUp() override {
    data_directory = traccc::io::data_directory();
  }

  std::string get_datafile(std::string name) { return data_directory + name; }
};
}  // namespace traccc::tests
