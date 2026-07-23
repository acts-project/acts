// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string_view>

namespace Acts {

/// @class IDesign
///
/// Abstract base class for the physical design description of a detector
/// element.
///
class IDesign {
 public:
  virtual ~IDesign() = default;

  /// Name identification for logging and debugging
  virtual std::string_view name() const = 0;

  const IDesign* design() const { return m_design.get(); }

  void assignDesign(std::shared_ptr<const IDesign> design) {
    m_design = std::move(design);
  }

 protected:
  std::shared_ptr<const IDesign> m_design{nullptr};
};

}  // namespace Acts
