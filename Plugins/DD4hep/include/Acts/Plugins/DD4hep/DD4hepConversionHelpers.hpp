// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <DD4hep/DetElement.h>
#include <DDRec/DetectorData.h>

namespace Acts {

/// Helper function to extract a parameter value from a dd4hep detector element
/// from VariantParameters
/// @tparam T The value type
/// @param key The key of the value to extract
/// @param elt The detector element instance
/// @return A copy of the value contained in the params instance
template <typename T>
T getParam(const std::string& key, dd4hep::DetElement& elt) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    throw std::runtime_error{"Detector Element has no VariantParameters"};
  }
  return params->get<T>(key);
}

/// Helper function to extract a VariantParameters instance
/// @param elt The detector element instance
/// @return The VariantParameters instance
inline dd4hep::rec::VariantParameters& getParams(dd4hep::DetElement& elt) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    throw std::runtime_error{"Detector Element has no VariantParameters"};
  }
  return *params;
}

/// Helper function to extract a VariantParameters instance, const version
/// @param elt The detector element instance
/// @return The VariantParameters instance
inline const dd4hep::rec::VariantParameters& getParams(
    const dd4hep::DetElement& elt) {
  const auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    throw std::runtime_error{"Detector Element has no VariantParameters"};
  }
  return *params;
}

/// Get a parameter value or an alternative value if either the
/// VariantParameters extension isn't set, or it doesn't contain the demanded
/// key
/// @tparam T The value type
/// @param key The key of the value to extract
/// @param elt The detector element instance
/// @param alternative The value to return if no params are set of the key doesn't exist
/// @return The value behind key, or @p alternative
template <typename T>
T getParamOr(const std::string& key, const dd4hep::DetElement& elt,
             T alternative) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    return alternative;
  }
  return params->value_or<T>(key, alternative);
}

/// Check if a detector element has a key set in its VariantParameters
/// @param key The key to check existance for
/// @param elt The detector element instance
/// @return True if the element has VariantParameters and the key exists, false if
///         either of these is not true
inline bool hasParam(const std::string& key, dd4hep::DetElement& elt) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    return false;
  }
  return params->contains(key);
}

/// Check if a detector element has VariantParameters set
/// @param elt The detector element instance
/// @return True if the VariantParameters exist, false if not
inline bool hasParams(dd4hep::DetElement& elt) {
  return elt.extension<dd4hep::rec::VariantParameters>(false) != nullptr;
}
}  // namespace Acts
