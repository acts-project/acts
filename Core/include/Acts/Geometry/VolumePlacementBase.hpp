// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"


namespace Acts{
    class GeometryContext;

    /// @brief Interface class to define the backend for alignable volumes 
    ///        that move coherently with the snesitive surfaces inside
    ///        The interface provides the transform from local -> global
    ///        coordinates as well as its way back.
    ///        To move the oriented surfaces along with the volume itself
    ///        the interface needs to provide a factory that creates the 
    ///        detector elements that are then passed to the oriented surfaces
    class VolumePlacementBase {
        public:
            /// @brief Default desctructor
            virtual ~VolumePlacementBase() = default;
            /// @brief Returns the transformation from the local volume coordinates to
            ///        the experiment's global coordinate system
            /// @param gctx The current geometry context object, e.g. alignment
            virtual const Transform3& localToGlobal(const GeometryContext& gctx) const = 0;
            /// @brief Returns the transformation from the experiment's gloabl frame to the
            ///        local volume coordinate system
            /// @param gctx The current geometry context object, e.g. alignment
            virtual const Transform3& globalToLocal(const GeometryContext& gctx) const = 0;

            /// @brief  Create internally a new detector element to align the portal surface
            ///         accordingly with the central volume transform
            /// @param faceIdx: Index of the oriented surface 
            /// @param internalTrf: Transform from the portal's frame to the volume's frame 
            virtual DetectorElementBase& createPortalSurfaceElement(const unsigned faceIdx,
                                                                    Transform3&& internalTrf) = 0;
            /// @brief Dispatch function that the client 
            /// @param portalElement 
            /// @param portalSurface 
            virtual void attachSurfaceToElement(DetectorElementBase& portalElement,
                                                std::shared_ptr<Surface> portalSurface) const = 0;

    };
}