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
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

class DD4hepDetectorElement;

/// @class DD4hepGeometryContext
///
/// @brief Geometry context for DD4hep
///
/// This class is used to store the geometry context for DD4hep
/// detector elements.
class DD4hepGeometryContext {
 public:
  using AlignmentStore =
      Delegate<const Acts::Transform3*(const DD4hepDetectorElement&)>;

  DD4hepGeometryContext() {}

  /// Check if a contextual transform is available for this detector element
  ///
  /// @param detElem The DD4hep DetElement which should be associated to
  ///
  /// If none is found, a null pointer is returned, which trigers the
  /// detector element to use its nominal transform.
  ///
  /// @return a pointer to the transform if found, otherwise nullptr
  const Transform3* contextualTransform(
      const DD4hepDetectorElement& detElem) const {
  
    if (m_aligmentStore.connected()) {
      return m_aligmentStore(detElem);
    }
    return nullptr;
  }

 private:
  AlignmentStore m_aligmentStore;
};
}  // namespace Acts
