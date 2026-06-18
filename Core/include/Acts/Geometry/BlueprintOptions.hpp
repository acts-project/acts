// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/NavigationPolicyFactory.hpp"

#include <memory>

namespace Acts::Experimental {

/// Options controlling blueprint navigation policies.
struct BlueprintOptions {
  /// Default navigation policy factory
  std::shared_ptr<NavigationPolicyFactory> defaultNavigationPolicyFactory{
      makeDefaultNavigationPolicyFactory()};

  /// If set, material designated on a portal face that must be merged during
  /// container stacking does not abort construction. Instead the offending
  /// material is discarded, the merged surface is tagged with a
  /// @ref MergedMaterialMarker, and a warning is emitted. This is lossy and
  /// intended as a debugging aid: the material designation should be moved to a
  /// face that is not merged (e.g. the enclosing container's face).
  bool keepGoingOnMaterialMergeFailure = false;

  /// Validates the blueprint options
  void validate() const;

 private:
  static std::unique_ptr<NavigationPolicyFactory>
  makeDefaultNavigationPolicyFactory();
};

}  // namespace Acts::Experimental
