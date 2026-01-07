// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <string>
#include <vector>

namespace Acts::Experimental {

class GbtsLutParser {
 public:
  explicit GbtsLutParser() = default;

  void parseLutFile(std::string& lutInputFile, bool useML);

  const std::vector<std::array<float, 5>>& getParsedLut() const {
    return m_mlLUT;
  }

 private:
  std::vector<std::array<float, 5>> m_mlLUT{};
};

}  // namespace Acts::Experimental
