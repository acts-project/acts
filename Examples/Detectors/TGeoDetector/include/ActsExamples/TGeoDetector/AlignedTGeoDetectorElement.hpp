// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Root/TGeoDetectorElement.hpp"
#include "ActsExamples/DetectorCommons/Aligned.hpp"

#include <memory>
#include <string>

class TGeoNode;
class TGeoMatrix;

namespace Acts {
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsExamples {
/// Define the aligned TGeo detector element and factory type
using AlignedTGeoDetectorElement =
    ActsExamples::Aligned<Acts::TGeoDetectorElement>;

/// @brief The factory for creating an aligned TGeo detector element
std::shared_ptr<AlignedTGeoDetectorElement> alignedTGeoDetectorElementFactory(
    const Acts::TGeoDetectorElement::Identifier& identifier,
    const TGeoNode& tGeoNode, const TGeoMatrix& tGeoMatrix,
    const std::string& axes, double scalor,
    std::shared_ptr<const Acts::ISurfaceMaterial> material);

}  // namespace ActsExamples
