// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsLutParser.hpp"

#include <fstream>

namespace Acts::Experimental {

void GbtsLutParser::parseLutFile(std::string& lutInputFile, bool useML) {
  if (useML) {
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
}
}  // namespace Acts::Experimental
