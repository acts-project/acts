// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/io/write.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cstdio>
#include <filesystem>
#include <fstream>

TEST(io_json, digitization_config) {
  // Read in the ODD digitization configuration.
  const traccc::digitization_config orig = traccc::io::read_digitization_config(
      traccc::io::get_absolute_path(
          (std::filesystem::path("geometries") / std::filesystem::path("odd") /
           std::filesystem::path("odd-digi-geometric-config.json"))
              .native()),
      traccc::data_format::json);

  // Write the digitization configuration to a temporary file.
  const std::string tempfile =
      (std::filesystem::temp_directory_path() /
       std::filesystem::path("traccc-odd-digi-config-io-test.json"))
          .native();
  traccc::io::write(tempfile, traccc::data_format::json, orig);

  // Read in the digitization configuration from the temporary file.
  const traccc::digitization_config copy =
      traccc::io::read_digitization_config(tempfile, traccc::data_format::json);

  // Remove the temporary file.
  std::remove(tempfile.c_str());

  // Check that the two configurations are the same.
  ASSERT_EQ(orig.size(), copy.size());
  auto orig_it = orig.begin();
  auto copy_it = copy.begin();
  for (; orig_it != orig.end(); ++orig_it, ++copy_it) {
    EXPECT_EQ(orig_it->dimensions, copy_it->dimensions);
    EXPECT_EQ(orig_it->bin_edges.size(), copy_it->bin_edges.size());
  }
}
