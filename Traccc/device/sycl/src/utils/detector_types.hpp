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
struct default_detector_kernel_tag {};
struct toy_detector_kernel_tag {};
struct telescope_detector_kernel_tag {};
struct odd_detector_kernel_tag {};
struct itk_detector_kernel_tag {};

template <typename T>
struct detector_tag_selector {};

template <>
struct detector_tag_selector<traccc::default_detector> {
  using type = default_detector_kernel_tag;
};

template <>
struct detector_tag_selector<traccc::toy_detector> {
  using type = toy_detector_kernel_tag;
};

template <>
struct detector_tag_selector<traccc::telescope_detector> {
  using type = telescope_detector_kernel_tag;
};

template <>
struct detector_tag_selector<traccc::odd_detector> {
  using type = odd_detector_kernel_tag;
};

template <>
struct detector_tag_selector<traccc::itk_detector> {
  using type = itk_detector_kernel_tag;
};

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
