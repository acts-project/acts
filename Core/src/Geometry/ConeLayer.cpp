// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ConeLayer.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

namespace Acts {

std::shared_ptr<ConeLayer> ConeLayer::create(
    const Transform3& transform, std::shared_ptr<const ConeBounds> cbounds,
    std::unique_ptr<SurfaceArray> surfaceArray, double thickness,
    std::unique_ptr<ApproachDescriptor> ad, LayerType laytyp) {
  return std::shared_ptr<ConeLayer>(
      new ConeLayer(transform, std::move(cbounds), std::move(surfaceArray),
                    thickness, std::move(ad), laytyp));
}

ConeLayer::ConeLayer(const Transform3& transform,
                     std::shared_ptr<const ConeBounds> cbounds,
                     std::unique_ptr<SurfaceArray> surfaceArray,
                     double thickness, std::unique_ptr<ApproachDescriptor> ade,
                     LayerType laytyp)
    : ConeSurface(transform, std::move(cbounds)),
      Layer(std::move(surfaceArray), thickness, std::move(ade), laytyp) {}

const ConeSurface& ConeLayer::surfaceRepresentation() const {
  return (*this);
}

ConeSurface& ConeLayer::surfaceRepresentation() {
  return (*this);
}

}  // namespace Acts
