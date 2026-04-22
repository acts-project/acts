// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray Core include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/core/detail/tuple_container.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace detray {

// Single store test
/// @{
using single_store_t = single_store<double, vecmem::vector, geometry_context>;
using single_store_dev_t =
    single_store<double, vecmem::device_vector, geometry_context>;
/// @}

// Tuple container test
/// @{
using tuple_cont_t = detail::tuple_container<dtuple, vecmem::vector<int>,
                                             vecmem::vector<double>>;
using tuple_cont_dev_t =
    detail::tuple_container<dtuple, vecmem::device_vector<int>,
                            vecmem::device_vector<double>>;
/// @}

// Regular multi store test (uses vectors as containers in every tuple element)
/// @{
enum class reg_type_ids : std::uint_least8_t {
  e_size = 0u,
  e_float = 1u,
  e_double = 2u,
};

using reg_multi_store_t =
    regular_multi_store<reg_type_ids, empty_context, dtuple, vecmem::vector,
                        std::size_t, float, double>;
using reg_multi_store_dev_t =
    regular_multi_store<reg_type_ids, empty_context, dtuple,
                        vecmem::device_vector, std::size_t, float, double>;
/// @}

/// Multi store test
/// @{

/// Test type that holds vecemem members and forces a hierarchical view/buffer
/// treatment

enum class type_ids : std::uint_least8_t {
  e_float = 0u,
  e_test_class = 1u,
};

template <template <typename...> class vector_t = dvector>
struct test_class {
  using view_type = dmulti_view<dvector_view<int>, dvector_view<double>>;
  using const_view_type =
      dmulti_view<dvector_view<const int>, dvector_view<const double>>;
  using buffer_type =
      dmulti_buffer<dvector_buffer<int>, dvector_buffer<double>>;

  DETRAY_HOST explicit test_class(vecmem::memory_resource* mr)
      : first(mr), second(mr) {}

  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE explicit test_class(view_t v)
      : first(detail::get<0>(v.m_view)), second(detail::get<1>(v.m_view)) {}

  DETRAY_HOST view_type get_data() {
    return view_type{vecmem::get_data(first), vecmem::get_data(second)};
  }

  vector_t<int> first;
  vector_t<double> second;
};

using multi_store_t =
    multi_store<type_ids, empty_context, dtuple, vecmem::vector<float>,
                test_class<vecmem::vector>>;
using multi_store_dev_t =
    multi_store<type_ids, empty_context, dtuple, vecmem::device_vector<float>,
                test_class<vecmem::device_vector>>;
/// @}

void test_single_store(typename single_store_t::view_type store_view,
                       vecmem::data::vector_view<double> sum_data);

void test_tuple_container(typename tuple_cont_t::view_type store_view,
                          vecmem::data::vector_view<double> sum_data);

void test_reg_multi_store(typename reg_multi_store_t::view_type store_view,
                          vecmem::data::vector_view<double> sum_data);

void test_multi_store(typename multi_store_t::view_type store_view,
                      vecmem::data::vector_view<double> sum_data);

}  // namespace detray
