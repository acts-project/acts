// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <stdexcept>

class G4MTRunManager;

namespace ActsTests::NativeWorkers {

class Exchange;

class TransportFailure final : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// Register the test telescope, transport-only physics and worker-local
/// actions. A geantino starts at x = request.x + the first seeded uniform draw
/// (mm). The recorded (x,z) positions are silicon entries at z = 29.75
/// and 59.75 mm. The callback runs on the native worker at the first silicon
/// entry.
void configure(G4MTRunManager& manager, Exchange& exchange,
               std::function<void(int)> onTransportStep);

}  // namespace ActsTests::NativeWorkers
