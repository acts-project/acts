// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Options {

/// Add common fitted related options
void addFittingOptions(boost::program_options::options_description& opt);

}  // namespace Options
}  // namespace ActsExamples
