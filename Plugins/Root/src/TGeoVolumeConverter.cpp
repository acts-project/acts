// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/TGeoVolumeConverter.hpp"

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"

#include "TGeoMatrix.h"
#include "TGeoShape.h"

using namespace Acts;

namespace ActsPlugins {

// @TODO: This currently uses @ref cylinderComponents to create a volume. This is probably not what we want.
std::unique_ptr<Acts::Volume> TGeoVolumeConverter::cylinderVolume(
    const TGeoShape& tgShape, const TGeoMatrix& tgTransform,
    double lengthScale) {
  auto [bounds, transform, thickness] =
      ActsPlugins::TGeoSurfaceConverter::cylinderComponents(
          tgShape, tgTransform.GetRotationMatrix(),
          tgTransform.GetTranslation(), "XY", lengthScale);

  if (bounds == nullptr) {
    throw std::invalid_argument(
        "TGeoShape -> CylinderVolume: could not convert TGeoShape to "
        "CylinderBounds.");
  }

  double medR = bounds->get(CylinderBounds::eR);
  double minR = medR - 0.5 * thickness;
  double maxR = medR + 0.5 * thickness;
  double hlZ = bounds->get(CylinderBounds::eHalfLengthZ);

  auto volBounds =
      std::make_shared<Acts::CylinderVolumeBounds>(minR, maxR, hlZ);

  return std::make_unique<Acts::Volume>(transform, volBounds);
}

}  // namespace ActsPlugins
