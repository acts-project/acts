// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"

#include <any>

namespace ActsExamples {

/// A detector element that is aligned with the ACTS framework
/// This class is a wrapper around the DetectorElement class
/// for different detectors and can be used with the same alignment
/// showcase infrastructure.
template <typename DetectorElement>
class Aligned : public DetectorElement {
 public:
  /// Using the constructor from the base class
  using DetectorElement::DetectorElement;
  /// An alignment aware transform call
  /// @param gctx the geometry context which is - if possible - unpacked to an AlignementContext
  /// @return The alignment-corrected transform if available, otherwise the nominal transform.
  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const final {
    if (gctx.hasValue()) {
      const auto* alignmentContext = gctx.maybeGet<AlignmentContext>();
      if (alignmentContext != nullptr && alignmentContext->store != nullptr) {
        const auto* transform =
            alignmentContext->store->contextualTransform(this->surface());
        // The store may only have a subset of surfaces registered
        if (transform != nullptr) {
          return *transform;
        }
      }
    }
    // If no alignment context is available, return the nominal transform
    return DetectorElement::nominalTransform();
  }
};

}  // namespace ActsExamples
