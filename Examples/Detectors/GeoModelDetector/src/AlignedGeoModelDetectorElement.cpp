// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/AlignedGeoModelDetectorElement.hpp"

std::shared_ptr<ActsExamples::AlignedGeoModelDetectorElement>
ActsExamples::alignedGeoModelDetectorElementFactory(
    const PVConstLink& geoPhysVol, std::shared_ptr<Acts::Surface> surface,
    const Acts::Transform3& sfTransform, double thickness) {
  return std::make_shared<AlignedGeoModelDetectorElement>(
      geoPhysVol, std::move(surface), sfTransform, thickness);
}
