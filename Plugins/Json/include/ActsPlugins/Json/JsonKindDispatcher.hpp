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
  /// Type of the dispatcher specialization
  using self_type = JsonKindDispatcher<return_t, args_t...>;
  /// Function signature type
  using decoder_signature = return_t(const nlohmann::json&, args_t...);
  /// Decoder callable type
  using decoder_type = std::function<decoder_signature>;
  /// Decoder callable pointer type
  using decoder_pointer_type = return_t (*)(args_t...);

  /// Explicit constructor of the dispatcher
  ///
  /// @param kindKey the key containing the type kind in json
  /// @param context context string for error signaling
  explicit JsonKindDispatcher(std::string kindKey = "kind",
                              std::string context = "JSON payload")
      : m_kindKey(std::move(kindKey)), m_context(std::move(context)) {
    if (m_kindKey.empty()) {
      throw std::invalid_argument(
          "JsonKindDispatcher kind key must be non-empty");
    }
    if (m_context.empty()) {
      m_context = "JSON payload";
    }
  }

  /// Register a kind and the corresponding decoder
  ///
  /// @param kind kind to register
  /// @param decoder corresponding decoder
  ///
  /// @return reference to this dispatcher instance
  self_type& registerKind(std::string kind, decoder_type decoder) {
    if (kind.empty()) {
      throw std::invalid_argument("JsonKindDispatcher kind must be non-empty");
    }
    if (!decoder) {
      throw std::invalid_argument("JsonKindDispatcher decoder must be valid");
    }
    auto [_, inserted] =
        m_decoders.emplace(std::move(kind), std::move(decoder));
    if (!inserted) {
      throw std::invalid_argument(
          "JsonKindDispatcher duplicate kind registration");
    }
    return *this;
  }

  /// Decode the registered kind from a json file
  ///
  /// @param encoded json file to decode
  /// @param args forwarding reference to the decoder arguments
  ///
  /// @return the object constructed from the json encoding
  template <typename... func_args_t>
  return_t operator()(const nlohmann::json& encoded,
                      func_args_t&&... args) const
    requires std::invocable<decoder_pointer_type, func_args_t...>
  {
    if (!encoded.contains(m_kindKey)) {
      throw std::invalid_argument("Missing '" + m_kindKey + "' in " +
                                  m_context);
    }

    const auto& kindValue = encoded.at(m_kindKey);
    if (!kindValue.is_string()) {
      throw std::invalid_argument("Invalid '" + m_kindKey + "' type in " +
                                  m_context);
    }

    const auto kind = kindValue.template get<std::string>();
    auto decoder = m_decoders.find(kind);
    if (decoder == m_decoders.end()) {
      throw std::invalid_argument("Unsupported " + m_context +
                                  " kind: " + kind);
    }

    if constexpr (std::is_void_v<return_t>) {
      decoder->second(encoded, std::forward<func_args_t>(args)...);
      return;
    } else {
      return decoder->second(encoded, std::forward<func_args_t>(args)...);
    }
  }

  /// Check if a certain kind is registered
  ///
  /// @param kind the kind to check for registration
  ///
  /// @return boolean showing registration status
  bool hasKind(std::string_view kind) const {
    return m_decoders.contains(std::string{kind});
  }

  /// Clear the registered decoders list
  void clear() { m_decoders.clear(); }

  /// Get the number of registered decoders
  ///
  /// @return number of registered decoders
  std::size_t size() const { return m_decoders.size(); }

 private:
  std::string m_kindKey;
  std::string m_context;
  std::unordered_map<std::string, decoder_type> m_decoders;
};

}  // namespace Acts
