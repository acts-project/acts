// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DetectorCommons/Aligned.hpp"
#include "ActsPlugins/Root/TGeoDetectorElement.hpp"

#include <memory>
#include <string>

class TGeoNode;
class TGeoMatrix;

namespace Acts {
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsExamples {
/// Define the aligned TGeo detector element and factory type
using AlignedTGeoDetectorElement = Aligned<ActsPlugins::TGeoDetectorElement>;

/// @brief The factory for creating an aligned TGeo detector element
std::shared_ptr<AlignedTGeoDetectorElement> alignedTGeoDetectorElementFactory(
    const ActsPlugins::TGeoDetectorElement::Identifier& identifier,
    const TGeoNode& tGeoNode, const TGeoMatrix& tGeoMatrix,
    const std::string& axes, double scalor,
    std::shared_ptr<const Acts::ISurfaceMaterial> material);

}  // namespace ActsExamples
