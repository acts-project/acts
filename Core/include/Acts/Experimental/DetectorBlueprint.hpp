// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Experimental/VolumeBlueprint.hpp"

#include <vector>
#include <tuple>

namespace Acts {

    /// This is the detector blueprint class, it contains a 
    /// building plan for building a detector from scratch
    /// using InternalBlueprint and VolumeBlueprint objects.
    class DetectorBlueprint {

        std::tuple<std::vector<VolumeBlueprint>, VolumeBlueprint::ContainerBuilder> manual;

    };

}

