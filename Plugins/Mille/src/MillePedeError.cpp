// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/MillePedeError.hpp"

#include <string>

namespace {

class MillePedeErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final { return "MillePedeError"; }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using ActsPlugins::ActsToMille::MillePedeError;

    switch (static_cast<MillePedeError>(c)) {
      case MillePedeError::InstallationNotFound:
        return "Installation of the 'pede' program was not found";
      case MillePedeError::SteeringNotFound:
        return "Steering file for `pede` was not found at configured location";
      case MillePedeError::SolverCrash:
        return "The solver crashed";
      case MillePedeError::InvalidSolution:
        return "The solver encountered a serious error and found no valid "
               "solution";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code ActsPlugins::ActsToMille::make_error_code(
    ActsPlugins::ActsToMille::MillePedeError e) {
  static MillePedeErrorCategory c;
  return {static_cast<int>(e), c};
}
