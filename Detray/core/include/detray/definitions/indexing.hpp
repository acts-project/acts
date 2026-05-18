// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/indexing.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray {

/// Default index type
using dindex = unsigned int;

/// Global invalid index definition
inline constexpr dindex dindex_invalid = detail::invalid_value<dindex>();

/// Index ranges and sequences
using dindex_range = detail::index_range<dindex>;
using dsized_index_range =
    detail::index_range<dindex, detail::sized_index_range>;
using dindex_sequence = dvector<dindex>;

/// Index that consists of multiple subindices
template <typename index_t = dindex, std::size_t DIM = 3u>
using dmulti_index = detail::multi_index<index_t, DIM>;

/// Link consisting of a type ID and an index
template <typename id_t = dindex, typename index_t = dindex,
          typename value_t = std::uint_least32_t, value_t id_mask = 0xf0000000,
          value_t index_mask = ~id_mask>
using dtyped_index =
    detail::typed_index<id_t, index_t, value_t, id_mask, index_mask>;

namespace detail {

/// Stub function to get a single index
template <std::size_t ID>
dindex get(dindex idx) noexcept {
  return idx;
}

}  // namespace detail

}  // namespace detray
