// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/TypeTraits.hpp"

#include <type_traits>

namespace Acts {

class Surface;

namespace concept {
  namespace detail_slc {
  template <typename T>
  using comparable_t = decltype(std::declval<T>() == std::declval<T>());

  template <typename T>
  using dereferenceable_t = decltype(*std::declval<T>());

  //~ template <typename T>
  //~ using surface_method_t = decltype(std::declval<T>().referenceSurface());

  template <typename T>
  struct SourceLinkConcept {
    constexpr static bool comparison_works =
        identical_to<bool, comparable_t, T>;
    static_assert(comparison_works,
                  "Source link does not implement equality operator");

    //~ constexpr static bool surface_method_exists =
        //~ converts_to<const Surface&, surface_method_t, T>;
    //~ static_assert(
        //~ surface_method_exists,
        //~ "Source link does not have compliant referenceSurface method");

    constexpr static bool copyable = std::is_copy_constructible_v<T>;
    static_assert(copyable, "Source link must be copy constructible");

    constexpr static bool default_constructible =
        std::is_default_constructible_v<T>;
    static_assert(default_constructible,
                  "Source link must be default-constructible");

    constexpr static bool value =
        concept ::require<comparison_works, copyable, //surface_method_exists,
                          default_constructible>;
  };
  }  // namespace detail_slc
}  // namespace concept
template <typename T>
constexpr bool SourceLinkConcept =
    concept ::detail_slc::SourceLinkConcept<T>::value;
}  // namespace Acts
