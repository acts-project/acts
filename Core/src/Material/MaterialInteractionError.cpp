// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialInteractionError.hpp"

#include <string>

namespace {
class MaterialInteractionErrorCategory : public std::error_category {
 public:
  const char* name() const noexcept final { return "MaterialInteractionError"; }

  std::string message(int code) const final {
    using Acts::MaterialInteractionError;
    switch (static_cast<MaterialInteractionError>(code)) {
      case MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval:
        return "Material noise removal would produce a "
               "non-positive-semidefinite covariance";
      default:
        return "unknown";
    }
  }
};
}  // namespace

std::error_code Acts::make_error_code(MaterialInteractionError error) {
  static MaterialInteractionErrorCategory category;
  return {static_cast<int>(error), category};
}
