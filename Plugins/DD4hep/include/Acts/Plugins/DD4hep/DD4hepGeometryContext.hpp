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

namespace Acts {

class DD4hepDetectorElement;

/// @class GeometryContext
///
/// @brief DD4hep specific geometry context for alignment handling
///
/// Extends the base GeometryContext to provide DD4hep-specific alignment
/// capabilities. The context can be active or inactive, controlling whether
/// alignment corrections should be applied.
///
/// @note This context is specifically designed to work with DD4hepDetectorElement
/// and provides contextual transformations for alignment purposes.
///
class DD4hepGeometryContext : public GeometryContext{
    public:
    /// Explicit inactive constructor
    /// Creates an inactive geometry context that skips alignment corrections
    /// @return DD4hepGeometryContext instance with inactive state
    static DD4hepGeometryContext inactive() { return DD4hepGeometryContext(false); }

    /// The transform of this detector element within the given context
    ///
    /// @param dElement The detector element
    ///
    /// @return The transform of the detector element
    const Transform3& contextualTransform( const DD4hepDetectorElement& dElement) const;

    /// @brief  Return the active status of the context
    /// @return boolean that indicates if the context is active
    bool isActive() const { return m_active; }

    private :
    /// Constructor
    explicit DD4hepGeometryContext(bool isContextActive) : m_active(isContextActive) {}

    const bool m_active = true;

};


}  // namespace Acts


