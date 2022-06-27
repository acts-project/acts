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

template <typename T>
T getParam(const std::string& key, dd4hep::DetElement& elt) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>();
  return params->get<T>(key);
}

inline dd4hep::rec::VariantParameters& getParams(dd4hep::DetElement& elt) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>();
  return *params;
}

inline const dd4hep::rec::VariantParameters& getParams(
    const dd4hep::DetElement& elt) {
  const auto* params = elt.extension<dd4hep::rec::VariantParameters>();
  return *params;
}

// template <typename T>
// T getParamOr(const std::string& key, dd4hep::DetElement& elt, T alternative)
// { auto* params = elt.extension<dd4hep::rec::VariantParameters>(false); if
// (params == nullptr) { return alternative;
// }
// return params->value_or<T>(key, alternative);
// }

template <typename T>
T getParamOr(const std::string& key, const dd4hep::DetElement& elt,
             T alternative) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    return alternative;
  }
  if (!params->contains(key)) {
    return alternative;
  }
  return params->get<T>(key);
}

inline bool hasParam(const std::string& key, dd4hep::DetElement& elt) {
  auto* params = elt.extension<dd4hep::rec::VariantParameters>(false);
  if (params == nullptr) {
    return false;
  }
  return params->contains(key);
}

inline bool hasParams(dd4hep::DetElement& elt) {
  return elt.extension<dd4hep::rec::VariantParameters>(false) != nullptr;
}
}  // namespace Acts
