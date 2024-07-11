// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#if defined(__cpp_concepts)

#include "Acts/Utilities/Concepts.hpp"

namespace Acts::Utils {

template <typename value_t>
void insert(std::output_iterator<value_t> auto storage, value_t&& value) {
  storage = std::forward<value_t>(value);
}

template <typename value_t>
void insert(Acts::isCollectionThatSupportsPushBack auto& storage,
            value_t&& value) {
  storage.push_back(std::forward<value_t>(value));
}

template <std::ranges::range storage_t, typename value_t>
requires(!Acts::isCollectionThatSupportsPushBack<storage_t> &&
         Acts::isCollectionThatSupportsInsert<
             storage_t>) void insert(storage_t& storage, value_t&& value) {
  storage.insert(std::ranges::end(storage), std::forward<value_t>(value));
}

template <typename... value_t>
void emplace(Acts::isCollectionThatSupportsEmplace<value_t...> auto& storage,
             value_t&&... args) {
  storage.emplace_back(std::forward<value_t>(args)...);
}

}  // namespace Acts::Utils

#else

#include <iterator>

namespace Acts::Utils {

template <template <typename...> typename container_t, typename value_t>
void insert(std::back_insert_iterator<container_t<value_t> > storage,
            value_t&& value) {
  storage = std::forward<value_t>(value);
}

}  // namespace Acts::Utils

#endif
