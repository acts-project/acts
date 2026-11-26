// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/AtlasStepperError.hpp"

#include <string>

namespace {
class AtlasStepperErrorCategory : public std::error_category {
 public:
  const char* name() const noexcept final { return "AtlasStepperError"; }

  std::string message(int error) const final {
    switch (static_cast<Acts::AtlasStepperError>(error)) {
      case Acts::AtlasStepperError::UnsupportedBoundParametrization:
        return "AtlasStepper requires global-spherical bound coordinates";
      default:
        return "unknown";
    }
  }
};
}  // namespace

std::error_code Acts::make_error_code(AtlasStepperError error) {
  static AtlasStepperErrorCategory category;
  return {static_cast<int>(error), category};
}
