// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/DetectorCommons/AlignedDetectorElement.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

namespace Acts {
    class PlanarBounds;
    class ISurfaceMaterial;
} // namespace Acts

namespace ActsExamples {

/// Define the aligned generic detector element and factory type
using AlignedGenericDetectorElement =
    ActsExamples::AlignedDetectorElement<ActsExamples::GenericDetectorElement>;

/// Factory function to create an aligned generic detector element
std::shared_ptr<AlignedGenericDetectorElement>
alignedGenericDetectorElementFactory(
    GenericDetectorElement::Identifier identifier,
    std::shared_ptr<const Acts::Transform3> transfomr,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material);

}  // namespace ActsExamples
