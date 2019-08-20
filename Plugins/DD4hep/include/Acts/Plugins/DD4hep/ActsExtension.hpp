// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ActsExtension.hpp, Acts project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <utility>
#include "DD4hep/Detector.h"

namespace Acts {

/// @class ActsExtension
///
/// @brief Extension of the \a %DD4hep \a DetElement needed for translation
/// into the Acts tracking geometry
///
/// This is a simple value / value container to be:
/// - filled by the Detector constructors if ACTS is needed
/// - interpreted by the DD4Hep layer builders and detector converters
class ActsExtension {
 public:
  /// Minimal constructor, sets the defail of the axis
  ActsExtension(const std::string& axes = "XYZ");

  /// Standard copy constructor
  ActsExtension(const ActsExtension& ext) = default;

  /// Copy constructor with element
  ActsExtension(const ActsExtension& ext, const dd4hep::DetElement& elem);

  /// Destructor
  ~ActsExtension() = default;

  /// Get the value
  ///
  /// @param tag the entry identifier in the value store
  /// @param type the (optional) category in the value store
  double getValue(const std::string& tag,
                  const std::string& category = "") const noexcept(false);

  /// Add the parameter to the store
  ///
  /// @param value the value to be added
  /// @param tag the entry identifier in the value store
  /// @param type the (optional) category in the value store
  void addValue(double value, const std::string& tag,
                const std::string& category = "");

  /// Check if the ActsExtension has a value (with optional category)
  ///
  /// @param type the primary identifier in the flag store
  /// @param type the (optional) category in the flag store
  bool hasValue(const std::string& tag, const std::string& category = "") const;

  /// Check if the ActsExtension has a value (with optional category)
  ///
  /// @param type the primary identifier in the flag store
  /// @param type the (optional) category in the flag store
  bool hasType(const std::string& type, const std::string& category = "") const;

  /// Add the characteristics
  ///
  /// @param type the primary identifier in the flag store
  /// @param category the (optional) category in the flag store
  /// @param the word to be stored
  void addType(const std::string& type, const std::string& category = "",
               const std::string& word = "");

  /// Get the string content
  ///
  /// @param type the primary identifier in the flag store
  /// @param category the (optional) category in the flag store
  const std::string getType(const std::string& type,
                            const std::string& category = "") const
      noexcept(false);

  /// Output to string
  std::string toString() const;

 private:
  /// Templated helper method
  template <typename T>
  void addT(std::map<std::string, T>& map, const T& val, const std::string& tag,
            const std::string& category, const T& catDeco);

  /// Templated helper method
  template <typename T>
  const T getT(const std::map<std::string, T>& map, const std::string& tag,
               const std::string& category = "") const noexcept(false);

  /// Templated helper method
  template <typename T>
  bool hasT(const std::map<std::string, T>& map, const std::string& tag,
            const std::string& category = "") const;

  /// Multiple flags to be stored, existance defines set
  std::map<std::string, std::string> m_flagStore;

  /// Unique value store for doubles
  std::map<std::string, double> m_valueStore;
};

// Templated helper method to get from the value/type store
template <typename T>
const T ActsExtension::getT(const std::map<std::string, T>& map,
                            const std::string& tag,
                            const std::string& category) const noexcept(false) {
  std::string ctag = "/";
  if (category != "") {
    ctag += category;
    ctag += "/";
  }
  ctag += tag;
  auto search = map.find(ctag);
  if (search == map.end()) {
    std::string error_message = "Acts::ActsExtension does not contain: ";
    error_message += ctag;
    throw std::runtime_error(error_message.c_str());
  }
  return search->second;
};

// Templated helper method to set from the value/type store
template <typename T>
void ActsExtension::addT(std::map<std::string, T>& map, const T& val,
                         const std::string& tag, const std::string& category,
                         const T& catDeco) {
  std::string ctag = "/";
  if (category != "") {
    ctag += category;
    map[ctag] = catDeco;
    ctag += "/";
  }
  ctag += tag;
  map[ctag] = val;
}

// Templated helper method to get from the value/type store
template <typename T>
bool ActsExtension::hasT(const std::map<std::string, T>& map,
                         const std::string& tag,
                         const std::string& category) const {
  std::string ctag = "/";
  if (category != "") {
    ctag += category;
    ctag += "/";
  }
  ctag += tag;
  auto search = map.find(ctag);
  return (search != map.end());
}

}  // namespace Acts
