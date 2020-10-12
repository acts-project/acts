// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include "HepMC3Writer.hpp"
#include "HepMC3Reader.hpp"

namespace ActsExamples {

namespace Options {

/// @brief HepMC3 specific writer options
///
/// @param desc The option description forward
void addHepMC3WriterOptions(Description& desc);

/// Read the HepMC3 writer options and @return a HepMC3WriterAscii::Config
///
/// @param variables is the parameter map for the options
///
/// @returns a Config object for the HepMC3WriterAscii
HepMC3WriterAscii::Config readHepMC3WriterOptions(
    const Variables& variables);
    
/// @brief HepMC3 specific reader options
///
/// @param desc The option description forward
void addHepMC3ReaderOptions(Description& desc);

/// Read the HepMC3 reader options and @return a HepMC3ReaderAscii::Config
///
/// @param variables is the parameter map for the options
///
/// @returns a Config object for the HepMC3ReaderAscii
HepMC3ReaderAscii::Config readHepMC3ReaderOptions(
    const Variables& variables);
}  // namespace Options
}  // namespace ActsExamples
