// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace boost {
namespace program_options {
class options_description;
class variables_map;
}  // namespace program_options
}  // namespace boost

namespace FW {
namespace Options {
using Description = ::boost::program_options::options_description;
using Variables = ::boost::program_options::variables_map;
}  // namespace Options
}  // namespace FW
