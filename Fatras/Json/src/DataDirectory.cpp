// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Json/DataDirectory.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace ActsFatras {

std::filesystem::path dataPath(const std::string_view name) {
  const std::array<std::filesystem::path, 2> directories{
      ACTS_FATRAS_DATA_DIR, ACTS_FATRAS_INSTALLED_DATA_DIR};

  std::string looked;
  for (const std::filesystem::path& directory : directories) {
    const std::filesystem::path candidate = directory / name;
    if (std::filesystem::exists(candidate)) {
      return std::filesystem::weakly_canonical(candidate);
    }
    looked += looked.empty() ? "" : " nor ";
    looked += candidate.string();
  }
  throw std::runtime_error("ActsFatras::dataPath: there is no " +
                           std::string(name) + "; looked in " + looked);
}

}  // namespace ActsFatras
