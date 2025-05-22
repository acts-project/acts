// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/AlignmentDelegate.hpp"

#include "Acts/Geometry/DetectorElementBase.hpp"

#include <iostream>

const Acts::Transform3* Acts::contextualTransform(
    const GeometryContext& gctx, const DetectorElementBase& detElement) {
  // Treating non-empty context => contextual alignment present
  if (gctx.hasValue()) {
    const AlignmentDelegate* alignment = gctx.maybeGet<AlignmentDelegate>();
    if (alignment != nullptr) {
      // Check if the alignment delegate is valid
      if (!alignment->connected()) {
        return nullptr;  // No alignment available
      }
    }
    // Check if a contextual transform is available for this detector element
    auto aTransform = (alignment != nullptr && alignment->connected())
                          ? (*alignment)(detElement.surface())
                          : nullptr;
    return aTransform;  // If nullptr - caller needs to handle this
  }
  // Empty context, return nullptr - caller needs to handle this
  return nullptr;
}
