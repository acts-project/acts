// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Where the descriptions, material and configurations that ship with Fatras
/// are.

#include <filesystem>
#include <string_view>

namespace ActsFatras {

/// The full path of one of the files shipped in `Fatras/data`.
///
/// The source tree is searched first and the installed copy second, so that a
/// build tree uses the files as they are checked out and an installed library
/// finds them where they were put. Everything that reads one of these also
/// takes a path outright, so nothing is forced through this.
///
/// @param name the file name, e.g. `itk-description.json`
/// @return the path to it
/// @throws std::runtime_error if there is no such file in either place, saying
///         where it looked
std::filesystem::path dataPath(std::string_view name);

}  // namespace ActsFatras
