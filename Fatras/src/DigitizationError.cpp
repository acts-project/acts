// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/DigitizationError.hpp"

const char* ActsFatras::detail::DigitizationErrorCategory::name() const
    noexcept {
  return "DigitizationError";
}

std::string ActsFatras::detail::DigitizationErrorCategory::message(
    int c) const {
  switch (static_cast<DigitizationError>(c)) {
    case DigitizationError::SmearOutOfBounds:
      return "Digitization: smeared out of surface bounds.";
    case DigitizationError::SmearError:
      return "Digitization: smearing error occured.";
    case DigitizationError::NoSurfaceDefined:
      return "Digitization: no surface for bound measurement defined.";
    default:
      return "unknown";
  }
}

const ActsFatras::detail::DigitizationErrorCategory&
ActsFatras::DigitizationErrorCategory() {
  static ActsFatras::detail::DigitizationErrorCategory c;
  return c;
}

std::error_code ActsFatras::make_error_code(ActsFatras::DigitizationError e) {
  return {static_cast<int>(e), ActsFatras::DigitizationErrorCategory()};
}
