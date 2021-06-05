// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

namespace Concepts {
namespace SourceLinkAccessor {

template <typename T>
using container_t = typename T::container_type;
template <typename T>
using key_t = typename T::key_type;
template <typename T>
using value_t = typename T::value_type;
template <typename T>
using const_iterator_t = typename T::const_iterator;

METHOD_TRAIT(count_t, count);
METHOD_TRAIT(equal_range_t, equal_range);
METHOD_TRAIT(at_t, at);

// clang-format off
    template <typename S>
      struct SourceLinkAccessorConcept {
        constexpr static bool container_type_exists = exists<container_t, S>;
        static_assert(container_type_exists, "Container type not found");
        constexpr static bool key_type_exists = exists<key_t, S>;
        static_assert(key_type_exists, "Key type not found");
        constexpr static bool value_type_exists = exists<value_t, S>;
        static_assert(value_type_exists, "Value type not found");
        constexpr static bool const_iterator_exists = exists<const_iterator_t, S>;
        static_assert(const_iterator_exists, "Const iterator not found");

        constexpr static bool count_exists = has_method<const S,
          size_t, count_t, const typename S::container_type&, const typename S::key_type&>;
        static_assert(count_exists, "count method not found");
        constexpr static bool equal_range_exists = has_method<const S,
          std::pair<typename S::const_iterator, typename S::const_iterator>,
          equal_range_t, const typename S::container_type&, const typename S::key_type&>;
        static_assert(equal_range_exists, "equal_range method not found");
        constexpr static bool at_exists = has_method<const S,
          const typename S::value_type&, at_t, const typename S::const_iterator&>;
        static_assert(at_exists, "at method not found");

        constexpr static bool value = require<container_type_exists,
                                              key_type_exists,
                                              value_type_exists,
                                              const_iterator_exists,
                                              count_exists,
                                              equal_range_exists,
                                              at_exists>;
      };
// clang-format on
}  // namespace SourceLinkAccessor
}  // namespace Concepts

template <typename accessor>
constexpr bool SourceLinkAccessorConcept =
    Acts::Concepts ::SourceLinkAccessor::SourceLinkAccessorConcept<
        accessor>::value;

}  // namespace Acts
