// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/io/read_magnetic_field.hpp"

// GTest include(s).
#include <gtest/gtest.h>

TEST(io_bfield, read_odd_binary_field) {
  // Read in the binary ODD file.
  traccc::magnetic_field field;
  traccc::io::read_magnetic_field(field, "geometries/odd/odd-bfield.cvf",
                                  traccc::data_format::binary);

  // Check that the field is about 2 Tesla in the center of the detector.
  auto field_view =
      field.as_view<traccc::host::inhom_bfield_backend_t<traccc::scalar>>();
  EXPECT_NEAR(
      static_cast<float>(field_view.at(0, 0, 0)[0]) / traccc::unit<float>::T,
      0.f, 0.01f);
  EXPECT_NEAR(
      static_cast<float>(field_view.at(0, 0, 0)[1]) / traccc::unit<float>::T,
      0.f, 0.01f);
  EXPECT_NEAR(
      static_cast<float>(field_view.at(0, 0, 0)[2]) / traccc::unit<float>::T,
      2.f, 0.01f);
}
