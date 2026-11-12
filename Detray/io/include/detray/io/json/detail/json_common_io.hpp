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
#include "detray/io/json/json.hpp"

namespace detray::io {

inline void to_json(nlohmann::ordered_json& j, const common_header_payload& h) {
  j["version"] = h.version;
  j["date"] = h.date;
  j["metadata"] = h.metadata;
  j["source"] = h.source;
}

inline void from_json(const nlohmann::ordered_json& j,
                      common_header_payload& h) {
  h.version = j["version"];
  h.date = j["date"];

  if (j.find("metadata") != j.end()) {
    h.metadata = j["metadata"];
  } else {
    h.metadata = "unknown";
  }

  if (j.find("source") != j.end()) {
    h.source = j["source"];
  } else {
    h.source = "unknown";
  }
}

inline void to_json(nlohmann::ordered_json& j, const tagged_header_payload& h) {
  const common_header_payload& ch = h;
  detray::io::to_json(j, ch);
  j["detector"] = h.detector;
  j["tag"] = h.tag;
}

inline void from_json(const nlohmann::ordered_json& j,
                      tagged_header_payload& h) {
  common_header_payload& ch = h;
  detray::io::from_json(j, ch);
  h.detector = j["detector"];
  h.tag = j["tag"];
}

inline void to_json(nlohmann::ordered_json& j, const header_payload& h) {
  const common_header_payload& ch = h;
  detray::io::to_json(j, ch);
  // Put additional converters for sub headers here ...
}

inline void from_json(const nlohmann::ordered_json& j, header_payload& h) {
  common_header_payload& ch = h;
  detray::io::from_json(j, ch);
  // Put additional converters for sub headers here ...
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
