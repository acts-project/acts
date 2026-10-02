// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string_view>

namespace Acts {

/// @class ISensorDesign
///
/// Abstract base class for the physical design description of a detector
/// element.
///
class ISensorDesign {
 public:
  virtual ~ISensorDesign() = default;

  /// Name identification for logging and debugging
  /// @return string view of the design name
  virtual std::string_view name() const = 0;
};

}  // namespace Acts
