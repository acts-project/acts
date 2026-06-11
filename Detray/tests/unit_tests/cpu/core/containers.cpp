// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/core/detail/tuple_container.hpp"

// Vecmem include(s)
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <tuple>
#include <vector>

using namespace detray;

constexpr float tol_single{1e-7f};
constexpr double tol_double{1e-31f};

struct test_func {
  template <typename container_t>
  auto operator()(const container_t& coll, const unsigned int /*index*/) {
    return coll.size();
  }
};

GTEST_TEST(detray_core, single_store) {
  // Vecmem memory resource
  vecmem::host_memory_resource resource;

  // Create store for double values
  single_store<double, vecmem::vector, geometry_context> store(resource);

  // Base store function check
  EXPECT_TRUE(store.empty());
  EXPECT_EQ(store.size(), 0u);

  // Add elements to the container
  geometry_context ctx{};
  store.reserve(4, ctx);
  store.emplace_back(ctx, 1.);
  store.push_back(2., ctx);
  store.insert(vecmem::vector<double>{10.5, 7.6}, ctx);

  EXPECT_FALSE(store.empty());
  EXPECT_EQ(store.size(), 4u);

  // Check  access to the data
  EXPECT_NEAR(store.at(0, ctx), 1., tol_double);
  EXPECT_NEAR(store.at(2, ctx), 10.5, tol_double);
  EXPECT_NEAR(store.at(1, ctx), 2., tol_double);
  EXPECT_NEAR(store.at(3, ctx), 7.6, tol_double);
}

GTEST_TEST(detray_core, tuple_container) {
  // Vecmem memory resource
  vecmem::host_memory_resource resource;

  using host_tuple_t =
      detail::tuple_container<std::tuple, vecmem::vector<float>,
                              vecmem::vector<std::size_t>,
                              vecmem::vector<int*>>;

  using device_tuple_t =
      detail::tuple_container<tuple, vecmem::device_vector<float>,
                              vecmem::device_vector<std::size_t>,
                              vecmem::device_vector<int*>>;

  vecmem::vector<float> vec1{};
  vecmem::vector<std::size_t> vec2{};
  vecmem::vector<int*> vec3{};

  host_tuple_t container(resource, vec1, vec2, vec3);

  static_assert(
      std::is_same_v<typename host_tuple_t::view_type,
                     dmulti_view<dvector_view<float>, dvector_view<std::size_t>,
                                 dvector_view<int*>>>,
      "View type incorrectly assembled");

  static_assert(std::is_same_v<typename host_tuple_t::const_view_type,
                               dmulti_view<dvector_view<const float>,
                                           dvector_view<const std::size_t>,
                                           dvector_view<int* const>>>,
                "Const view type incorrectly assembled");

  typename host_tuple_t::view_type view = get_data(container);
  device_tuple_t dev_container(view);

  // Base container function check
  EXPECT_EQ(container.size(), 3u);
  EXPECT_EQ(dev_container.size(), 3u);

  EXPECT_TRUE(container.get<0>().empty());
  EXPECT_TRUE(container.get<1>().empty());
  EXPECT_TRUE(detail::get<2>(container).empty());
}

GTEST_TEST(detray_core, vector_multi_store) {
  enum class type_ids : std::uint_least8_t {
    e_size = 0u,
    e_float = 1u,
    e_double = 2u,
  };

  using enum type_ids;

  using container_t =
      regular_multi_store<type_ids, empty_context, std::tuple, vecmem::vector,
                          std::size_t, float, double>;

  // Vecmem memory resource
  vecmem::host_memory_resource resource;

  // Create tuple vector container
  container_t vector_store(resource);

  // Base container function check
  EXPECT_EQ(vector_store.n_collections(), 3u);
  EXPECT_EQ(vector_store.empty<e_size>(), true);
  EXPECT_EQ(vector_store.empty<e_float>(), true);
  EXPECT_EQ(vector_store.empty<e_double>(), true);
  EXPECT_EQ(vector_store.size<e_size>(), 0u);
  EXPECT_EQ(vector_store.size<e_float>(), 0u);
  EXPECT_EQ(vector_store.size<e_double>(), 0u);

  // Add elements to the container
  vector_store.push_back<e_size>(1u);
  vector_store.emplace_back<e_size>(empty_context{}, 2u);
  vector_store.push_back<e_float>(3.1f);
  vector_store.emplace_back<e_float>(empty_context{}, 4.5f);
  vector_store.push_back<e_double>(5.5);
  vector_store.emplace_back<e_double>(empty_context{}, 6.f);

  EXPECT_EQ(vector_store.empty<e_size>(), false);
  EXPECT_EQ(vector_store.empty<e_float>(), false);
  EXPECT_EQ(vector_store.empty<e_double>(), false);
  EXPECT_EQ(vector_store.size<e_size>(), 2u);
  EXPECT_EQ(vector_store.size<e_float>(), 2u);
  EXPECT_EQ(vector_store.size<e_double>(), 2u);

  vecmem::vector<std::size_t> int_vec{3u, 4u, 5u};
  vector_store.insert(int_vec);

  vecmem::vector<float> float_vec{12.1f};
  vector_store.insert(float_vec);

  vector_store.insert(vecmem::vector<double>{10.5, 7.6});

  // int collectiont
  EXPECT_EQ(vector_store.size<e_size>(), 5u);
  EXPECT_EQ(vector_store.get<e_size>()[0], 1u);
  EXPECT_EQ(vector_store.get<e_size>()[1], 2u);
  EXPECT_EQ(vector_store.get<e_size>()[2], 3u);
  EXPECT_EQ(vector_store.get<e_size>()[3], 4u);
  EXPECT_EQ(vector_store.get<e_size>()[4], 5u);

  // float collectiont
  EXPECT_EQ(vector_store.size<e_float>(), 3u);
  EXPECT_NEAR(vector_store.get<e_float>()[0], 3.1f, tol_single);
  EXPECT_NEAR(vector_store.get<e_float>()[1], 4.5f, tol_single);
  EXPECT_NEAR(vector_store.get<e_float>()[2], 12.1f, tol_single);

  // double collectiont
  EXPECT_EQ(vector_store.size<e_double>(), 4u);
  EXPECT_NEAR(vector_store.get<e_double>()[0], 5.5, tol_double);
  EXPECT_NEAR(vector_store.get<e_double>()[1], 6., tol_double);
  EXPECT_NEAR(vector_store.get<e_double>()[2], 10.5, tol_double);
  EXPECT_NEAR(vector_store.get<e_double>()[3], 7.6, tol_double);

  // call functor
  container_t::single_link l0{e_size, 0u};
  container_t::single_link l1{e_float, 0u};
  container_t::single_link l2{e_double, 0u};
  EXPECT_EQ(vector_store.visit<test_func>(l0), 5u);
  EXPECT_EQ(vector_store.visit<test_func>(l1), 3u);
  EXPECT_EQ(vector_store.visit<test_func>(l2), 4u);
}
