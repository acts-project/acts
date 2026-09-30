/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "traccc/seeding/grids/axis.hpp"
#include "traccc/seeding/grids/serializer2.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace traccc;

GTEST_TEST(traccc_grid2, serialize) {
  vecmem::host_memory_resource resource;

  axis2::regular<> r6{6u, -3.f, 7.f, resource};
  axis2::circular<> c12{12u, -3.f, 3.f, resource};

  serializer2 ser2;

  // Serializing
  unsigned int test = ser2.serialize(r6, c12, 0u, 0u);
  EXPECT_EQ(test, 0u);
  test = ser2.serialize(r6, c12, 5u, 0u);
  EXPECT_EQ(test, 5u);
  test = ser2.serialize(r6, c12, 0u, 1u);
  EXPECT_EQ(test, 6u);
  test = ser2.serialize(r6, c12, 5u, 2u);
  EXPECT_EQ(test, 17u);
}
