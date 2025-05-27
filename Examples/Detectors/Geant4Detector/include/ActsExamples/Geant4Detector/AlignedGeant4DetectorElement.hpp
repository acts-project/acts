// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "ActsExamples/DetectorCommons/Aligned.hpp"

#include <memory>

class G4VPhysicalVolume;

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsExamples {
/// Define the aligned DD4hep detector element and factory type
using AlignedGeant4DetectorElement =
    ActsExamples::Aligned<Acts::Geant4DetectorElement>;

/// Factory function to create an aligned Geant4 detector element
/// It forwards the arguments to the Geant4DetectorElement constructor
std::shared_ptr<AlignedGeant4DetectorElement>
alignedGeant4DetectorElementFactory(std::shared_ptr<Acts::Surface> surface,
                                    const G4VPhysicalVolume& g4physVol,
                                    const Acts::Transform3& toGlobal,
                                    double thickness);

}  // namespace ActsExamples
