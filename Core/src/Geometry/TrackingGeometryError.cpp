// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometryError.hpp"

#include <string>

namespace {

class TrackingGeometryErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "TrackingGeometryError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::TrackingGeometryError;

    switch (static_cast<TrackingGeometryError>(c)) {
      case TrackingGeometryError::PositionNotOnAssociatedSurface:
        return "Position is not on the associated boundary or portal surface.";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::TrackingGeometryError e) {
  static TrackingGeometryErrorCategory c;
  return {static_cast<int>(e), c};
}
