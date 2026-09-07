// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/ISensorDesign.hpp"

#include <memory>

namespace Acts::detail {

/// @brief Mixin base for placement elements that carry a sensor design.
///        Cast target for DigitizationAlgorithm; lives in detail so it
///        can be promoted or removed in follow-up work without ABI concerns.
class ISensorDesignHolder {
 public:
  virtual ~ISensorDesignHolder() = default;

  /// @return pointer to the attached sensor design, or nullptr
  virtual const ISensorDesign* sensorDesign() const = 0;

  /// @brief Attach a sensor design to this element
  virtual void assignSensorDesign(
      const std::shared_ptr<const ISensorDesign>& design) const = 0;
};

}  // namespace Acts::detail
