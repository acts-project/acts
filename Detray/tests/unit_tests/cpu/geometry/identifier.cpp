// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s).
#include "detray/geometry/identifier.hpp"

// Google test include(s).
#include <gtest/gtest.h>

using namespace detray;

/// Test retrieval of surface from collection using brute force searching
GTEST_TEST(detray_geometry, identifier) {
  auto geo_id = geometry::identifier{};

  // Check a empty identifier
  EXPECT_EQ(geo_id.volume(), static_cast<dindex>((1UL << 12) - 1UL));
  EXPECT_EQ(geo_id.id(), static_cast<surface_id>((1UL << 4) - 1UL));
  EXPECT_EQ(geo_id.index(), static_cast<dindex>((1UL << 21) - 1UL));
  EXPECT_EQ(geo_id.transform(), static_cast<dindex>((1UL << 21) - 1UL));
  EXPECT_EQ(geo_id.extra(), static_cast<dindex>((1UL << 6) - 1UL));

  geo_id.set_volume(2UL)
      .set_id(surface_id::e_passive)
      .set_index(42UL)
      .set_transform(11UL)
      .set_extra(24UL);

  // Check the values after setting them
  EXPECT_EQ(geo_id.volume(), 2UL);
  EXPECT_EQ(geo_id.id(), surface_id::e_passive);
  EXPECT_EQ(geo_id.index(), 42UL);
  EXPECT_EQ(geo_id.transform(), 11UL);
  EXPECT_EQ(geo_id.extra(), 24UL);

  // Check invalid identifier
  EXPECT_FALSE(geo_id.is_invalid());
  geo_id.set_volume((1UL << 12) - 1UL);
  EXPECT_TRUE(geo_id.is_invalid());
  geo_id.set_volume(2UL);
  EXPECT_FALSE(geo_id.is_invalid());

  geo_id.set_id(static_cast<surface_id>((1UL << 4) - 1UL));
  EXPECT_TRUE(geo_id.is_invalid());
  geo_id.set_id(surface_id::e_passive);
  EXPECT_FALSE(geo_id.is_invalid());

  geo_id.set_index((1UL << 21) - 1UL);
  EXPECT_TRUE(geo_id.is_invalid());
  geo_id.set_index(42UL);
  EXPECT_FALSE(geo_id.is_invalid());

  geo_id.set_transform((1UL << 21) - 1UL);
  EXPECT_TRUE(geo_id.is_invalid());
  geo_id.set_transform(11UL);
  EXPECT_FALSE(geo_id.is_invalid());

  geo_id.set_extra((1UL << 6) - 1UL);
  EXPECT_FALSE(geo_id.is_invalid());
}
