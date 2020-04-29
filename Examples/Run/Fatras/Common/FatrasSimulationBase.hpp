// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "ACTFW/Utilities/OptionsFwd.hpp"

namespace Acts {
class TrackingGeometry;
}

namespace FW {

class Sequencer;
class RandomNumbers;

/// Setup the simulation algorithm and related algorithms.
///
/// @param variables user configuration variables map
/// @param sequencer the framework sequencer
/// @param randomNumbers the random numbers tools
/// @param trackingGeometry the tracking geometry
void setupSimulation(
    const Options::Variables& variables, Sequencer& sequencer,
    std::shared_ptr<const RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry);

}  // namespace FW
