// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/json/json.hpp"

namespace detray::io {

inline void to_json(nlohmann::ordered_json& j, const common_header_payload& h) {
  j["version"] = h.version;
  j["detector"] = h.detector;
  j["date"] = h.date;
  j["tag"] = h.tag;
}

inline void from_json(const nlohmann::ordered_json& j,
                      common_header_payload& h) {
  h.version = j["version"];
  h.detector = j["detector"];
  h.date = j["date"];
  h.tag = j["tag"];
}

inline void to_json(nlohmann::ordered_json& j, const header_payload<bool>& h) {
  j["common"] = h.common;
  // Do write the optional subheader here, but in the dedicated serializers
}

inline void from_json(const nlohmann::ordered_json& j,
                      header_payload<bool>& h) {
  h.common = j["common"];
  // Do not look at the optional subheader here
}

/// Data links IO
/// @{
inline void to_json(nlohmann::ordered_json& j, const single_link_payload& so) {
  j = so.link;
}

inline void from_json(const nlohmann::ordered_json& j,
                      single_link_payload& so) {
  so.link = j;
}

template <typename type_id>
inline void to_json(nlohmann::ordered_json& j,
                    const typed_link_payload<type_id>& m) {
  j["type"] = static_cast<unsigned int>(m.type);
  j["index"] = m.index;
}

template <typename type_id>
inline void from_json(const nlohmann::ordered_json& j,
                      typed_link_payload<type_id>& m) {
  m.type = static_cast<type_id>(j["type"]);
  m.index = j["index"];
}
/// @}

}  // namespace detray::io
