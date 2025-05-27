// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsExamples/DetectorCommons/AlignedDetectorElement.hpp"

#include <memory>
#include <string>

namespace Acts {
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsExamples {
/// Define the aligned DD4hep detector element and factory type
using AlignedDD4hepDetectorElement =
    ActsExamples::AlignedDetectorElement<Acts::DD4hepDetectorElement>;

/// Factory function to create an aligned DD4hep detector element
/// It forwards the arguments to the DD4hepDetectorElement constructor
std::shared_ptr<AlignedDD4hepDetectorElement>
alignedDD4hepDetectorElementFactory(
    const dd4hep::DetElement detElement, const std::string& axes, double scalor,
    bool isDisc,
    std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

}  // namespace ActsExamples