// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>

#include "ACTFW/GenericDetector/GenericDetectorElement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace FW {

namespace Contextual {

/// @class PayloadDetectorElement extends GenericDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
///
/// The PayloadDetectorElement demonstrates how a GeometryContext
/// can be used if it carries the entire set of Transforms through
/// the program flow.
///
/// The nominal transform is only used to once create the alignment
/// store and then in a contextual call the actual detector element
/// position is taken from the alignment Store.
///
/// In this simple implementation, it does rely on the Identifier
/// to be orderded from 0 to N-1, as the identifier is simply taken
/// as a vector index for the alignment store
class PayloadDetectorElement : public Generic::GenericDetectorElement {
 public:
  /// @class ContextType
  /// convention: nested to the Detector element
  struct ContextType {
    // The alignment store of this event
    // not the fastest, but good enough for a demonstrator
    std::vector<Acts::Transform3D> alignmentStore;
  };

  /// Constructor for an alignable surface
  ///
  /// @note see Generic::GenericDetectorElement for documentation
  template <typename... Args>
  PayloadDetectorElement(Args&&... args)
      : Generic::GenericDetectorElement(std::forward<Args>(args)...) {}

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform(gctx)
  const Acts::Transform3D& transform(
      const Acts::GeometryContext& gctx) const final override;
};

inline const Acts::Transform3D& PayloadDetectorElement::transform(
    const Acts::GeometryContext& gctx) const {
  // cast into the right context object
  auto alignContext = std::any_cast<ContextType>(gctx);
  identifier_type idValue = identifier_type(identifier());

  // check if we have the right alignment parameter in hand
  if (idValue < alignContext.alignmentStore.size()) {
    return alignContext.alignmentStore[idValue];
  }
  // Return the standard transform if not found
  return GenericDetectorElement::transform(gctx);
}

}  // end of namespace Contextual
}  // end of namespace FW
