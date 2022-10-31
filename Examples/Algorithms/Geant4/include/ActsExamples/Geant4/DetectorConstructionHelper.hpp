// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Detector/IBaseDetector.hpp"

#include <G4VUserDetectorConstruction.hh>

namespace ActsExamples {

G4VUserDetectorConstruction *getG4DetectorContruction(
    const IBaseDetector &detector);

}
