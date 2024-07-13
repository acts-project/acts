// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#if defined(__cpp_concepts)

#include <ranges>

namespace Acts::details {

// Define some concepts
template <typename external_t>
concept isCollectionThatSupportsPushBack =
    std::ranges::range<external_t> && requires {
  typename external_t::value_type;
} && requires(external_t coll, typename external_t::value_type val) {
  coll.push_back(val);
};

template <typename external_t>
concept isCollectionThatSupportsInsert =
    std::ranges::range<external_t> && requires {
  typename external_t::value_type;
} && requires(external_t coll, typename external_t::value_type val) {
  coll.insert(std::ranges::end(coll), val);
};

template <typename external_t, typename... args_t>
concept isCollectionThatSupportsEmplace = std::ranges::range<external_t> &&
    requires(external_t coll, args_t... vals) {
  coll.emplace_back(vals...);
};

// Define some functions
template <typename value_t>
void insert(std::output_iterator<value_t> auto& storage, value_t&& value) {
  storage = std::forward<value_t>(value);
}

template <typename value_t>
void insert(Acts::details::isCollectionThatSupportsPushBack auto& storage,
            value_t&& value) {
  storage.push_back(std::forward<value_t>(value));
}

template <std::ranges::range storage_t, typename value_t>
requires(!Acts::details::isCollectionThatSupportsPushBack<storage_t> &&
         Acts::details::isCollectionThatSupportsInsert<
             storage_t>) void insert(storage_t& storage, value_t&& value) {
  storage.insert(std::ranges::end(storage), std::forward<value_t>(value));
}

template <typename... value_t>
void emplace(
    Acts::details::isCollectionThatSupportsEmplace<value_t...> auto& storage,
    value_t&&... args) {
  storage.emplace_back(std::forward<value_t>(args)...);
}

}  // namespace Acts::details

#endif
