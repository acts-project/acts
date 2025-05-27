// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "ActsExamples/DetectorCommons/Aligned.hpp"

#include <memory>

class PVConstLink;

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsExamples {
/// Define the aligned DD4hep detector element and factory type
using AlignedGeoModelDetectorElement =
    ActsExamples::Aligned<Acts::GeoModelDetectorElement>;

/// Factory function to create an aligned GeoModel detector element
/// It forwards the arguments to the GeoModelDetectorElement constructor
using AlignedGeoModelDetectorElementFactory =
    [](PVConstLink geoPhysVol, std::shared_ptr<Surface> surface,
       const Transform3& sfTransform, double thickness) {
      return std::make_shared<AlignedGeoModelDetectorElement>(
          geoPhysVol, std::move(surface), sfTransform, thickness);
    };

}  // namespace ActsExamples
