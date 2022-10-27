// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

/// Add HepMC3 specific writer options.
///
/// @param desc The option description forward
void addHepMC3WriterOptions(Description& desc);

/// Read the HepMC3 writer options and @return a HepMC3WriterAscii::Config
///
/// @param variables is the parameter map for the options
///
/// @returns a Config object for the HepMC3WriterAscii
HepMC3AsciiWriter::Config readHepMC3WriterOptions(const Variables& variables);

/// Add HepMC3 specific reader options.
///
/// @param desc The option description forward
void addHepMC3ReaderOptions(Description& desc);

/// Read the HepMC3 reader options and @return a HepMC3ReaderAscii::Config
///
/// @param variables is the parameter map for the options
///
/// @returns a Config object for the HepMC3ReaderAscii
HepMC3AsciiReader::Config readHepMC3ReaderOptions(const Variables& variables);

}  // namespace Options
}  // namespace ActsExamples
