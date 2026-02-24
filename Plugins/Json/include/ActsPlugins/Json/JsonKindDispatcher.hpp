// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include <nlohmann/json.hpp>

namespace Acts {

/// JSON-specific dispatcher that routes decoder functions based on a string
/// kind tag found in the encoded payload.
template <typename return_t, typename... args_t>
class JsonKindDispatcher {
 public:
  using decoder_type = std::function<return_t(const nlohmann::json&, args_t...)>;
  using self_type = JsonKindDispatcher<return_t, args_t...>;

  explicit JsonKindDispatcher(std::string kindKey = "kind",
                              std::string context = "JSON payload")
      : m_kindKey(std::move(kindKey)), m_context(std::move(context)) {
    if (m_kindKey.empty()) {
      throw std::invalid_argument("JsonKindDispatcher kind key must be non-empty");
    }
    if (m_context.empty()) {
      m_context = "JSON payload";
    }
  }

  self_type& registerKind(std::string kind, decoder_type decoder) {
    if (kind.empty()) {
      throw std::invalid_argument("JsonKindDispatcher kind must be non-empty");
    }
    if (!decoder) {
      throw std::invalid_argument("JsonKindDispatcher decoder must be valid");
    }
    auto [_, inserted] = m_decoders.emplace(std::move(kind), std::move(decoder));
    if (!inserted) {
      throw std::invalid_argument("JsonKindDispatcher duplicate kind registration");
    }
    return *this;
  }

  return_t operator()(const nlohmann::json& encoded, args_t... args) const {
    if (!encoded.contains(m_kindKey)) {
      throw std::invalid_argument("Missing '" + m_kindKey + "' in " + m_context);
    }

    const auto& kindValue = encoded.at(m_kindKey);
    if (!kindValue.is_string()) {
      throw std::invalid_argument("Invalid '" + m_kindKey + "' type in " +
                                  m_context);
    }

    const auto kind = kindValue.template get<std::string>();
    auto decoder = m_decoders.find(kind);
    if (decoder == m_decoders.end()) {
      throw std::invalid_argument("Unsupported " + m_context + " kind: " +
                                  kind);
    }

    if constexpr (std::is_void_v<return_t>) {
      decoder->second(encoded, std::forward<args_t>(args)...);
      return;
    } else {
      return decoder->second(encoded, std::forward<args_t>(args)...);
    }
  }

  bool hasKind(std::string_view kind) const {
    return m_decoders.find(std::string{kind}) != m_decoders.end();
  }

  void clear() { m_decoders.clear(); }

  std::size_t size() const { return m_decoders.size(); }

 private:
  std::string m_kindKey;
  std::string m_context;
  std::unordered_map<std::string, decoder_type> m_decoders;
};

}  // namespace Acts
