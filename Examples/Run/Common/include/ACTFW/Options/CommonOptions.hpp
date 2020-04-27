// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <boost/program_options.hpp>
#include <string>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"

namespace FW {
namespace Options {

/// Construct the options description with minimal default options.
///
/// @param caption Optional help text caption
boost::program_options::options_description makeDefaultOptions(
    std::string caption = std::string());

/// Add sequencer options, e.g. number of events
void addSequencerOptions(boost::program_options::options_description& opt);

/// Add random number options such as the global seed.
void addRandomNumbersOptions(boost::program_options::options_description& opt);

/// Add common geometry-related options.
void addGeometryOptions(boost::program_options::options_description& opt);

/// Add common material-related options.
void addMaterialOptions(boost::program_options::options_description& opt);

/// Add common input-related options.
void addInputOptions(boost::program_options::options_description& opt);

/// Add common output-related options.
void addOutputOptions(boost::program_options::options_description& opt);

/// Parse options and return the resulting variables map.
///
/// Automatically prints the help text if requested.
///
/// @returns Empty variables map if help text was shown.
boost::program_options::variables_map parse(
    const boost::program_options::options_description& opt, int argc,
    char* argv[]);

/// Read the log level.
Acts::Logging::Level readLogLevel(
    const boost::program_options::variables_map& vm);

/// Read the sequencer config.
Sequencer::Config readSequencerConfig(
    const boost::program_options::variables_map& vm);

// Read the random numbers config.
RandomNumbers::Config readRandomNumbersConfig(
    const boost::program_options::variables_map& vm);

}  // namespace Options
}  // namespace FW
