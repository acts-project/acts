// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"

#include <array>
#include <string>
#include <vector>
#include <fstream>

#include "actsvg/meta.hpp"

namespace Acts {

namespace Svg {

struct Style {
  // Fill parameters
  std::array<int, 3> fillColor = {0, 0, 0};
  ActsScalar fillOpacity = 1.;

  // Highlight parameters
  std::array<int, 3> highlightColor = {0, 0, 0};
  std::vector<std::string> highlights = {};

  ActsScalar strokeWidth = 0.5;
  std::array<int, 3> strokeColor = {0, 0, 0};

  unsigned int nSegments = 72u;
};


/// Helper method to write to file
///
/// @param objects to be written out
/// @param fileName the file name is to be given
///
inline static void toFile(const std::vector<actsvg::svg::object>& objects, const std::string& fileName) {

    actsvg::svg::file foutFile;
    
    for (const auto& o : objects){
        foutFile.add_object(o);
    }

    std::ofstream fout;
    fout.open(fileName);
    fout << foutFile;
    fout.close();

}


}  // namespace Svg

}  // namespace Acts