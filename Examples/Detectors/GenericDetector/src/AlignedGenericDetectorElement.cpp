// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/AlignedGenericDetectorElement.hpp"

std::shared_ptr<ActsExamples::AlignedGenericDetectorElement>
ActsExamples::alignedGenericDetectorElementFactory(
    GenericDetectorElement::Identifier identifier,
    const Acts::Transform3& transform,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<AlignedGenericDetectorElement>(
      identifier, std::move(transform), std::move(pBounds), thickness,
      std::move(material));
}
