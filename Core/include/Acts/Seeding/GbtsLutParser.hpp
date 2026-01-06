// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Acts::Experimental {

class GbtsLutParser {
 public:
  explicit GbtsLutParser(std::string& lutInputFile) {
    if (lutInputFile.empty()) {
      throw std::runtime_error("Cannot find ML predictor LUT file");
    } else {
      m_mlLUT.reserve(100);
      std::ifstream ifs(std::string(lutInputFile).c_str());

      if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open LUT file");
      }

      float cl_width{}, min1{}, max1{}, min2{}, max2{};

      while (ifs >> cl_width >> min1 >> max1 >> min2 >> max2) {
        std::array<float, 5> lut_line = {cl_width, min1, max1, min2, max2};
        m_mlLUT.emplace_back(lut_line);
      }
      if (!ifs.eof()) {
        // ended if parse error present, not clean EOF

        throw std::runtime_error("Stopped reading LUT file due to parse error");
      }

      ifs.close();
    }
  }

  const std::vector<std::array<float, 5>>& getParsedLut() const {
    return m_mlLUT;
  }

 private:
  std::vector<std::array<float, 5>> m_mlLUT{};
};

}  // namespace Acts::Experimental
