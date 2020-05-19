// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ConeLayer.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <utility>

Acts::ConeLayer::ConeLayer(std::shared_ptr<const Transform3D> transform,
                           std::shared_ptr<const ConeBounds> cbounds,
                           std::unique_ptr<SurfaceArray> surfaceArray,
                           double thickness,
                           std::unique_ptr<ApproachDescriptor> ade,
                           LayerType laytyp)
    : ConeSurface(std::move(transform), std::move(cbounds)),
      Layer(std::move(surfaceArray), thickness, std::move(ade), laytyp) {}

const Acts::ConeSurface& Acts::ConeLayer::surfaceRepresentation() const {
  return (*this);
}

Acts::ConeSurface& Acts::ConeLayer::surfaceRepresentation() {
  return (*this);
}
