// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "ActsExamples/DetectorCommons/Aligned.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"

#include <memory>
#include <string>

/// @cond
namespace dd4hep {
class DetElement;
}  // namespace dd4hep
/// @endcond

namespace ActsExamples {
/// Define the aligned DD4hep detector element and factory type
using AlignedDD4hepDetectorElement =
    Aligned<ActsPlugins::DD4hepDetectorElement>;

/// Factory function to create an aligned DD4hep detector element
/// It forwards the arguments to the DD4hepDetectorElement constructor
std::shared_ptr<AlignedDD4hepDetectorElement>
alignedDD4hepDetectorElementFactory(
    const dd4hep::DetElement detElement, const std::string& axes, double scalor,
    bool isDisc,
    std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

}  // namespace ActsExamples
