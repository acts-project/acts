// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {

class Sequencer;

namespace Options {

/// Add magnetic field options with a `bf-` prefix.
void addMagneticFieldOptions(Description& desc);

/// Setup any additional services required for the magnetic field.
void setupMagneticFieldServices(const Variables& vars, Sequencer& seq);

/// Read and create the magnetic field from the given user variables.
std::shared_ptr<Acts::MagneticFieldProvider> readMagneticField(
    const Variables& vars);

}  // namespace Options
}  // namespace ActsExamples
