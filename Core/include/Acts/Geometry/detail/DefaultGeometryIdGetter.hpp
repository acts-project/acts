// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryID.hpp"

namespace Acts {
namespace detail {

/// Uniform access to geometry identifiers for different types.
struct DefaultGeometryIdGetter {
  // support objects that implement `.geoID()`.
  template <typename object_t>
  inline auto operator()(const object_t& object) const
      -> decltype(object.geoID(), GeometryID()) {
    return object.geoID();
  }
  // support objects that implement `.geometryId()`.
  template <typename object_t>
  inline auto operator()(const object_t& object) const
      -> decltype(object.geometryId(), GeometryID()) {
    return object.geometryId();
  }
};

}  // namespace detail
}  // namespace Acts
