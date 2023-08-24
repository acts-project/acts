// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"

namespace Acts {
namespace detail {

/// @brief Helper method to get an attribute with fallback
///
/// @note the fallback value has to be provided
///
/// @tparam value_type the primitive type allowed by variant parameters
///
/// @param node the node object from DD4hep
/// @param attrName the name of the attribute that is checked
/// @param fallbackValue the fallbackValue
///
/// @return either the gathered attribute or the fallback
template <typename value_type>
value_type getAttrValueOr(const dd4hep::xml::Component& node,
                          const std::string& attrName,
                          const value_type& fallbackValue) {
  if (node.hasAttr(dd4hep::xml::Strng_t(attrName.c_str()))) {
    return node.attr<value_type>(attrName.c_str());
  } else {
    return fallbackValue;
  }
}

}  // namespace detail
}  // namespace Acts
