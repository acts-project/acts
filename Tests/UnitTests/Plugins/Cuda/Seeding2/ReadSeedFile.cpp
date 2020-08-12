// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Local include(s).
#include "ReadSeedFile.hpp"

// System include(s).
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

/// Difference allowed on floating point numbers to still be treated equal
static constexpr float allowedDiff = std::numeric_limits<float>::epsilon() * 4;

std::vector<std::unique_ptr<TestSpacePoint> > readSeedFile(
    const std::string& fileName, bool filterDuplicates) {
  // The result object.
  std::vector<std::unique_ptr<TestSpacePoint> > result;

  // Open the input file.
  std::ifstream spFile(fileName);
  if (!spFile.is_open()) {
    throw std::runtime_error("Could not open file: " + fileName);
  }

  // Read the file's lines one by one, and create spacepoint objects out of
  // them.
  std::size_t duplicatesFound = 0;
  while (!spFile.eof()) {
    std::string line;
    std::getline(spFile, line);
    std::stringstream ss(line);
    std::string linetype;

    ss >> linetype;
    if (linetype != "lxyz") {
      continue;
    }

    int layer;
    float x, y, z, varianceR, varianceZ;
    ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
    const float r = std::sqrt(x * x + y * y);
    const float f22 = varianceR;
    const float wid = varianceZ;
    float cov = wid * wid * .08333;

    if (cov < f22) {
      cov = f22;
    }
    if (std::abs(z) > 450.) {
      varianceZ = 9. * cov;
      varianceR = .06;
    } else {
      varianceR = 9. * cov;
      varianceZ = .06;
    }

    // Create a new spacepoint object.
    std::unique_ptr<TestSpacePoint> sp(
        new TestSpacePoint{x, y, z, r, layer, varianceR, varianceZ});

    // Check if we already have another spacepoint with the same coordinates.
    if (filterDuplicates) {
      for (const auto& otherSP : result) {
        if ((std::abs(sp->x() - otherSP->x()) < allowedDiff) &&
            (std::abs(sp->y() - otherSP->y()) < allowedDiff) &&
            (std::abs(sp->z() - otherSP->z()) < allowedDiff)) {
          ++duplicatesFound;
          continue;
        }
      }
    }

    // Store the new spacepoint.
    result.push_back(std::move(sp));
  }
  // Tell the user how many duplicates were found.
  if (duplicatesFound) {
    std::cerr << duplicatesFound << " duplicates found in: " << fileName
              << std::endl;
  }

  return result;
}
