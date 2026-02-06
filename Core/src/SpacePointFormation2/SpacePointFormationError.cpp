// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation2/SpacePointFormationError.hpp"

#include <string>

namespace {

class SpacePointFormationErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "SpacePointFormationError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::SpacePointFormationError;

    switch (static_cast<SpacePointFormationError>(c)) {
      case SpacePointFormationError::ClusterPairDistanceExceeded:
        return "Cluster pair distance exceeded";
      case SpacePointFormationError::ClusterPairThetaDistanceExceeded:
        return "Cluster pair theta distance exceeded";
      case SpacePointFormationError::ClusterPairPhiDistanceExceeded:
        return "Cluster pair phi distance exceeded";
      case SpacePointFormationError::CosmicToleranceNotMet:
        return "Cosmic tolerance not met";
      case SpacePointFormationError::OutsideLimits:
        return "Outside limits";
      case SpacePointFormationError::OutsideRelaxedLimits:
        return "Outside relaxed limits";
      case SpacePointFormationError::NoSolutionFound:
        return "No solution found";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::SpacePointFormationError e) {
  static SpacePointFormationErrorCategory c;
  return {static_cast<int>(e), c};
}
