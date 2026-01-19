// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Convert a geometry hierarchy map to/from Json.
///
/// @tparam value_t value type stored in the geometry hierarchy map
///
/// The value type is expected to be directly convertible to/from a Json object.
/// It has to be either a fundamental type or appropriate `to_json(json&, const
/// value_t&)` and `from_json(const json&, value_t&)` functions must be
/// available. See the relevant [nlohmann::json documentation][1] for further
/// information.
///
/// A user-defined identifier is stored in the encoded Json object that is used
/// to identify which value type is stored in the file. The identifier is
/// checked for consistency when decoding the Json object.
///
/// [1]: https://nlohmann.github.io/json/features/arbitrary_types/
template <typename value_t,
          class decorator_t = void /* e.g. ITrackingGeometryJsonDecorator */>
class GeometryHierarchyMapJsonConverter {
 public:
  using Value = value_t;
  using Container = GeometryHierarchyMap<value_t>;

  /// Construct the converter.
  ///
  /// @param valueIdentifier user-defined identifier for the stored value
  explicit GeometryHierarchyMapJsonConverter(std::string valueIdentifier)
      : m_valueIdentifier(std::move(valueIdentifier)) {
    if (m_valueIdentifier.empty()) {
      throw std::invalid_argument("Value identifier must be non-empty");
    }
  }

  /// Encode the geometry hierarchy map into a Json object.
  ///
  /// @param container Geometry hierarchy map that should be encoded
  /// @param decorator nullptr or a decorator to add extra values to the json
  /// output
  /// @return Encoded Json object
  nlohmann::json toJson(const Container& container,
                        const decorator_t* decorator) const;

  /// Decode a Json object into a geometry hierarchy map.
  ///
  /// @param encoded Json object that should be decoded
  /// @return Decoded geometry hierarchy map
  /// @throw std::invalid_argument in case of format errors
  Container fromJson(const nlohmann::json& encoded) const;

 private:
  static constexpr const char* kHeaderKey = "acts-geometry-hierarchy-map";
  static constexpr const char* kEntriesKey = "entries";
  /// The version of the encoded Json container format. This must be increased
  /// manually every time the container format changes.
  static constexpr int kFormatVersion = 0;

  std::string m_valueIdentifier;
};

// implementations

// auxiliary struct to indicate a missing specialization of a template which
// requires specialisation
template <typename T, class decorator_t>
struct missing_specialization : std::false_type {};

// methods to adapt type decorations for the given decorator
template <class T, class decorator_t>
void decorateJson([[maybe_unused]] const decorator_t* decorator,
                  [[maybe_unused]] const T& src,
                  [[maybe_unused]] nlohmann::json& dest) {
  // this needs to be specialised
  static_assert(
      missing_specialization<T, decorator_t>::value,
      "Explicit specialization needed for each decorator_t and src T");
}
template <class T, class decorator_t>
void decorateJson([[maybe_unused]] const decorator_t* decorator,
                  [[maybe_unused]] const T* src,
                  [[maybe_unused]] nlohmann::json& dest) {
  // this needs to be specialised
  static_assert(
      missing_specialization<T, decorator_t>::value,
      "Explicit specialization needed for each decorator_t and src T");
}

template <typename value_t, class decorator_t>
nlohmann::json GeometryHierarchyMapJsonConverter<value_t, decorator_t>::toJson(
    const Container& container,
    [[maybe_unused]] const decorator_t* decorator) const {
  // encode header
  nlohmann::json encoded = nlohmann::json::object();
  encoded[kHeaderKey] = nlohmann::json::object();
  encoded[kHeaderKey]["format-version"] = kFormatVersion;
  encoded[kHeaderKey]["value-identifier"] = m_valueIdentifier;
  // encode entries
  nlohmann::json entries = nlohmann::json::array();
  for (std::size_t i = 0; i < container.size(); ++i) {
    auto entry =
        GeometryIdentifierJsonConverter::encodeIdentifier(container.idAt(i));
    auto value_json = nlohmann::json(container.valueAt(i));
    if constexpr (!std::is_same_v<decorator_t, void>) {
      decorateJson(decorator, container.valueAt(i), value_json);
    }
    entry["value"] = std::move(value_json);
    entries.push_back(std::move(entry));
  }
  encoded[kEntriesKey] = std::move(entries);
  return encoded;
}

template <typename value_t, class decorator_t>
auto GeometryHierarchyMapJsonConverter<value_t, decorator_t>::fromJson(
    const nlohmann::json& encoded) const -> Container {
  // verify json format header
  auto header = encoded.find(kHeaderKey);
  if (header == encoded.end()) {
    throw std::invalid_argument(
        "Missing header entry in json geometry hierarchy map");
  }
  if (header->at("format-version").get<int>() != kFormatVersion) {
    throw std::invalid_argument(
        "Invalid format version in json geometry hierarchy map");
  }
  if (header->at("value-identifier").get<std::string>() != m_valueIdentifier) {
    throw std::invalid_argument(
        "Inconsistent value identifier in Json geometry hierarchy map");
  }
  // decode json entries
  if (!encoded.contains(kEntriesKey)) {
    throw std::invalid_argument(
        "Missing entries in json geometry hierarchy map");
  }
  std::vector<std::pair<GeometryIdentifier, Value>> elements;
  for (const auto& entry : encoded.at(kEntriesKey)) {
    auto id = GeometryIdentifierJsonConverter::decodeIdentifier(entry);
    auto value = entry.at("value").get<Value>();
    elements.emplace_back(id, std::move(value));
  }
  return Container(std::move(elements));
}

/// @}
}  // namespace Acts
