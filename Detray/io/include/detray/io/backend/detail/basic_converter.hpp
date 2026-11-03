// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/utils/io_metadata.hpp"

// System include(s)
#include <string_view>

// Convert basic information like links and header data to and from payloads
namespace detray::io::detail::basic_converter {

/// @returns a link from its io payload @param link_data
inline dindex from_payload(const single_link_payload& link_data) {
  return static_cast<dindex>(link_data.link);
}

/// Convert a link @param idx into its io payload
inline single_link_payload to_payload(const std::size_t idx) {
  single_link_payload link_data;
  link_data.link = idx;

  return link_data;
}

/// Convert a typed link with a type id @param id and index @param idx into its
/// io payload
template <typename type_id>
inline typed_link_payload<type_id> to_payload(const type_id id,
                                              const std::size_t idx) {
  typed_link_payload<type_id> link_data;

  link_data.type = id;
  link_data.index = idx;

  return link_data;
}

}  // namespace detray::io::detail::basic_converter
