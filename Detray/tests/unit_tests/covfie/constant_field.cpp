// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// covfie core
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// GTest include(s)
#include <gtest/gtest.h>

GTEST_TEST(Covfie, ConstantField1D) {
  using field_t =
      covfie::field<covfie::backend::constant<covfie::vector::float1,
                                              covfie::vector::float1>>;

  field_t f(
      covfie::make_parameter_pack(field_t::backend_t::configuration_t{2.f}));
  field_t::view_t v(f);

  constexpr float start{-100.f};
  constexpr float step{1.f};
  for (int i = 0; i < 200; ++i) {
    const float x{start + static_cast<float>(i) * step};

    EXPECT_EQ(v.at(x)[0], 2.f);
  }
}

GTEST_TEST(Covfie, ConstantField2D) {
  using field_t =
      covfie::field<covfie::backend::constant<covfie::vector::float2,
                                              covfie::vector::float2>>;

  field_t f(covfie::make_parameter_pack(
      field_t::backend_t::configuration_t{2.f, 5.f}));
  field_t::view_t v(f);

  constexpr float start{-100.f};
  constexpr float step{1.f};
  for (int i = 0; i < 200; ++i) {
    for (int j = 0; j < 200; ++j) {
      const float x{start + static_cast<float>(i) * step};
      const float y{start + static_cast<float>(j) * step};

      EXPECT_EQ(v.at(x, y)[0], 2.f);
      EXPECT_EQ(v.at(x, y)[1], 5.f);
    }
  }
}

GTEST_TEST(Covfie, ConstantField3D) {
  using field_t =
      covfie::field<covfie::backend::constant<covfie::vector::float3,
                                              covfie::vector::float3>>;

  field_t f(covfie::make_parameter_pack(
      field_t::backend_t::configuration_t{2.f, 5.f, -4.f}));
  field_t::view_t v(f);

  constexpr float start{-10.f};
  constexpr float step{1.f};
  for (int i = 0; i < 20; ++i) {
    for (int j = 0; j < 20; ++j) {
      for (int k = 0; k < 20; ++k) {
        const float x{start + static_cast<float>(i) * step};
        const float y{start + static_cast<float>(j) * step};
        const float z{start + static_cast<float>(k) * step};

        EXPECT_EQ(v.at(x, y, z)[0], 2.f);
        EXPECT_EQ(v.at(x, y, z)[1], 5.f);
        EXPECT_EQ(v.at(x, y, z)[2], -4.f);
      }
    }
  }
}
