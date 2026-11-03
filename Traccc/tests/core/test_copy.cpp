// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Projection include(s).
#include "traccc/edm/container.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Google test include(s).
#include <gtest/gtest.h>

TEST(ContainerCopy, HostToHost) {
  using test_container_types = traccc::container_types<int, int>;

  // Memory resources used by the application.
  vecmem::host_memory_resource host_mr;

  // Create a host container
  test_container_types::host host_orig{&host_mr};

  // Fill the host container
  for (int i = 0; i < 3; ++i) {
    host_orig.push_back(i,
                        test_container_types::host::vector_type<int>{i, i + 1});
  }

  // Copy the container to the host buffer, and back.
  auto host_buffer = traccc::get_data(host_orig);
  test_container_types::device host_copy{host_buffer};

  // Check the copied container
  ASSERT_EQ(host_copy.size(), 3u);
  const auto& host_headers = host_copy.get_headers();
  ASSERT_EQ(host_headers.size(), 3u);
  ASSERT_EQ(host_headers[0], 0);
  ASSERT_EQ(host_headers[1], 1);
  ASSERT_EQ(host_headers[2], 2);

  const auto& host_items = host_copy.get_items();
  ASSERT_EQ(host_items.size(), 3u);
  ASSERT_EQ(host_items[0].size(), 2u);
  ASSERT_EQ(host_items[0][0], 0);
  ASSERT_EQ(host_items[0][1], 1);
  ASSERT_EQ(host_items[1].size(), 2u);
  ASSERT_EQ(host_items[1][0], 1);
  ASSERT_EQ(host_items[1][1], 2);
  ASSERT_EQ(host_items[2].size(), 2u);
  ASSERT_EQ(host_items[2][0], 2);
  ASSERT_EQ(host_items[2][1], 3);
}
