// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/AlignedTGeoDetectorElement.hpp"

std::shared_ptr<ActsExamples::AlignedTGeoDetectorElement>
ActsExamples::alignedTGeoDetectorElementFactory(
    const Acts::TGeoDetectorElement::Identifier& identifier,
    const TGeoNode& tGeoNode, const TGeoMatrix& tGeoMatrix,
    const std::string& axes, double scalor,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<AlignedTGeoDetectorElement>(
      identifier, tGeoNode, tGeoMatrix, axes, scalor, std::move(material));
};
