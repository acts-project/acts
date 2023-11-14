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

class Surface;

namespace Concepts {
namespace SourceLinkAccessor {

template <typename T>
using container_t = typename T::Container;
template <typename T>
using key_t = typename T::Key;
template <typename T>
using value_t = typename T::Value;
template <typename T>
using iterator_t = typename T::Iterator;

METHOD_TRAIT(count_t, count);
METHOD_TRAIT(range_t, range);
METHOD_TRAIT(at_t, at);

// clang-format off
    template <typename S>
      struct SourceLinkAccessorConcept {
        constexpr static bool container_exists = exists<container_t, S>;
        static_assert(container_exists, "Container type not found");
        constexpr static bool key_exists = exists<key_t, S>;
        static_assert(key_exists, "Key type not found");
        constexpr static bool value_exists = exists<value_t, S>;
        static_assert(value_exists, "Value type not found");
        constexpr static bool iterator_exists = exists<iterator_t, S>;
        static_assert(iterator_exists, "Iterator type not found");

        constexpr static bool container_pointer_exists =
          std::is_same_v<std::decay_t<decltype(*(std::declval<S>().container))>, container_t<S>>;
        static_assert(container_pointer_exists, "Pointer to container not found");

        constexpr static bool count_exists = has_method<const S,
          std::size_t, count_t, const Surface&>;
        static_assert(count_exists, "count method not found");
        constexpr static bool range_exists = has_method<const S,
          std::pair<typename S::Iterator, typename S::Iterator>,
          range_t, const Surface&>;
        static_assert(range_exists, "range method not found");
        constexpr static bool at_exists = has_method<const S,
          const typename S::Value&, at_t, const typename S::Iterator&>;
        static_assert(at_exists, "at method not found");

        constexpr static bool value = require<container_exists,
                                              key_exists,
                                              value_exists,
                                              container_pointer_exists,
                                              iterator_exists,
                                              count_exists,
                                              range_exists,
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
