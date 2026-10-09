// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/BoundParametrizationError.hpp"

#include <string>

namespace {
class BoundParametrizationErrorCategory : public std::error_category {
 public:
  const char* name() const noexcept final {
    return "BoundParametrizationError";
  }

  std::string message(int code) const final {
    using Acts::BoundParametrizationError;
    switch (static_cast<BoundParametrizationError>(code)) {
      case BoundParametrizationError::SingularDirection:
        return "Direction Jacobian is singular or insufficiently conditioned";
      case BoundParametrizationError::UndefinedCanonicalization:
        return "Canonicalization derivative is undefined at a pole";
      default:
        return "unknown";
    }
  }
};
}  // namespace

std::error_code Acts::make_error_code(BoundParametrizationError error) {
  static BoundParametrizationErrorCategory category;
  return {static_cast<int>(error), category};
}
