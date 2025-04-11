// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <span>

namespace HepMC3 {
class GenEvent;
}

namespace Acts {
class Logger;
}

namespace ActsExamples::HepMC3Util {

std::shared_ptr<HepMC3::GenEvent> mergeEvents(
    std::span<const HepMC3::GenEvent*> genEvents, const Acts::Logger& logger);

}
