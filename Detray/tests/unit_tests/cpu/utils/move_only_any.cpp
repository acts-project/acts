// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/utils/move_only_any.hpp"

#include <cstdint>
#include <memory>

#include <gtest/gtest.h>

namespace {
struct alignas(128) aligned_payload {
  std::unique_ptr<int> value = std::make_unique<int>(42);
};
struct throwing_payload {
  throwing_payload() = default;
  throwing_payload(throwing_payload&&) {
    throw std::runtime_error("move failed");
  }
};
}  // namespace

TEST(detray_move_only_any, lifetime_alignment_and_reset) {
  detray::move_only_any value{aligned_payload{}};
  EXPECT_EQ(
      reinterpret_cast<std::uintptr_t>(&value.as<aligned_payload>()) % 128, 0u);
  EXPECT_EQ(*value.as<aligned_payload>().value, 42);
  auto* address = &value.as<aligned_payload>();
  detray::move_only_any moved{std::move(value)};
  EXPECT_FALSE(value.has_value());
  EXPECT_FALSE(value.is<aligned_payload>());
  EXPECT_EQ(&moved.as<aligned_payload>(), address);
  value = std::move(moved);
  EXPECT_FALSE(moved.has_value());
  EXPECT_EQ(&value.as<aligned_payload>(), address);
  value.reset();
  value.reset();
  EXPECT_FALSE(value.has_value());
  EXPECT_THROW(value.type(), std::logic_error);
  value.set(std::make_unique<int>(7));
  EXPECT_EQ(*value.as<std::unique_ptr<int>>(), 7);
  EXPECT_THROW(value.as<int>(), std::bad_any_cast);
}

TEST(detray_move_only_any, failed_replacement_preserves_value) {
  detray::move_only_any value{std::make_unique<int>(17)};
  throwing_payload replacement;
  EXPECT_THROW(value.set(std::move(replacement)), std::runtime_error);
  EXPECT_EQ(*value.as<std::unique_ptr<int>>(), 17);
}

TEST(detray_move_only_any, releases_replaced_and_reset_values) {
  auto first = std::make_shared<int>(1);
  std::weak_ptr<int> first_lifetime = first;
  detray::move_only_any value{std::move(first)};
  auto second = std::make_shared<int>(2);
  std::weak_ptr<int> second_lifetime = second;
  value.set(std::move(second));
  EXPECT_TRUE(first_lifetime.expired());
  EXPECT_FALSE(second_lifetime.expired());
  value.reset();
  EXPECT_TRUE(second_lifetime.expired());
}
