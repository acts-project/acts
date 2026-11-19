// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <type_traits>

#include <gtest/gtest.h>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector_buffer.hpp"
#include "detray/detectors/telescope_metadata.hpp"
#include "detray/test/framework/types.hpp"

namespace {
using metadata = detray::test::default_metadata;
using detector = detray::host::detector<metadata>;
using telescope =
    detray::host::detector<detray::telescope_metadata<detray::test::algebra>>;
using supported = std::tuple<telescope, detector>;
}  // namespace

TEST(detray_erased_detector, ownership_and_dispatch) {
  static_assert(!std::is_copy_constructible_v<detray::host_detector>);
  static_assert(std::is_nothrow_move_constructible_v<detray::host_detector>);
  vecmem::host_memory_resource mr;
  detray::host_detector erased;
  EXPECT_FALSE(erased.is<metadata>());
  auto visit = []<typename T>(const T& concrete) -> const std::type_info& {
    return typeid(concrete);
  };
  EXPECT_THROW(detray::host_detector_visitor<supported>(erased, visit),
               std::invalid_argument);
  erased.set<metadata>(detector{mr});
  EXPECT_TRUE(erased.is<detector>());
  EXPECT_EQ(&detray::host_detector_visitor<supported>(erased, visit),
            &typeid(detector));
  EXPECT_THROW(
      detray::host_detector_visitor<std::tuple<telescope>>(erased, visit),
      std::invalid_argument);
  EXPECT_THROW(detray::host_detector_visitor<std::tuple<>>(erased, visit),
               std::invalid_argument);
  auto moved = std::move(erased);
  EXPECT_FALSE(erased.has_value());
  EXPECT_TRUE(moved.is<metadata>());
  moved.set<telescope>(telescope{mr});
  EXPECT_TRUE(moved.is<telescope>());
  EXPECT_THROW(moved.as<metadata>(), std::bad_any_cast);
  moved.reset();
  EXPECT_FALSE(moved.has_value());
}

TEST(detray_erased_detector, buffer_dispatch) {
  vecmem::host_memory_resource mr;
  vecmem::copy copy;
  detray::host_detector erased;
  erased.set<metadata>(detector{mr});
  auto buffer = detray::buffer_from_host_detector<supported>(erased, mr, copy);
  EXPECT_TRUE(buffer.is<metadata>());
  auto called = detray::detector_buffer_visitor<supported>(
      buffer, []<typename T>(const detray::detector_view_t<T>& view) {
        detray::detector_device_t<T> concrete{view};
        return concrete.volumes().empty();
      });
  EXPECT_TRUE(called);
  auto moved = std::move(buffer);
  EXPECT_FALSE(buffer.has_value());
  EXPECT_TRUE(moved.is<metadata>());
}
