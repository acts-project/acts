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

    class VolumePlacementBase {
        public:
            virtual ~VolumePlacementBase() = default;

            virtual const Transform3& localToGlobal(const GeometryContext& gctx) const = 0;

            virtual const Transform3& globalToLocal(const GeometryContext& gctx) const = 0;

            virtual DetectorElementBase& createPortalSurfaceElement(const unsigned faceIdx,
                                                                    Transform3&& internalTrf) = 0;

            virtual void attachSurfaceToElement(DetectorElementBase& portalElement,
                                                std::shared_ptr<Surface> portalSurface) const = 0;

    };
}