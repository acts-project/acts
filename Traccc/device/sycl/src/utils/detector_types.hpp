// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/definitions/primitives.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_type_list.hpp"

namespace traccc::sycl {

/*
 * SYCL requires a little bit of extra massaging to make the kernel tags
 * work...
 */

template <typename T>
struct detector_tag_selector {};

/// Kernel tags for all available detector types
/// (see traccc::available_detector_type_list)
// clang-format off
#define TRACCC_SYCL_GENERATED_DETECTOR_TAG(NAME)            \
  struct NAME##_detector_kernel_tag {};                     \
  template <>                                               \
  struct detector_tag_selector<traccc::NAME##_detector> {   \
    using type = NAME##_detector_kernel_tag;                \
  };
DETRAY_FOR_EACH_METADATA(TRACCC_SYCL_GENERATED_DETECTOR_TAG)
#undef TRACCC_SYCL_GENERATED_DETECTOR_TAG
// clang-format on

template <typename T>
using detector_tag_selector_t = typename detector_tag_selector<T>::type;

template <typename T>
concept detector_tag_exists_for_backend =
    requires { typename detector_tag_selector_t<T>; };

template <typename>
struct detector_tag_existance_validator {};

template <typename... Ts>
struct detector_tag_existance_validator<std::tuple<Ts...>> {
  static constexpr bool value = (detector_tag_exists_for_backend<Ts> && ...);
};

static_assert(
    detector_tag_existance_validator<traccc::detector_type_list>::value,
    "Not all detector types registered for SYCL have an accompanying tag");

}  // namespace traccc::sycl
