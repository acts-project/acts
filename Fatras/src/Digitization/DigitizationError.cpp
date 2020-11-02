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
    case DigitizationError::SmearingOutOfRange:
      return "Smeared out of surface bounds.";
    case DigitizationError::SmearingError:
      return "Smearing error occured.";
    case DigitizationError::UndefinedSurface:
      return "Surface undefined for this operation.";
    case DigitizationError::MaskingError:
      return "Surface mask could not be applied.";
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
