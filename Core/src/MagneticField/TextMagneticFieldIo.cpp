// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/TextMagneticFieldIo.hpp"

#include "Acts/MagneticField/BFieldMapUtils.hpp"

#include <fstream>
#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>

namespace {
constexpr std::size_t kDefaultSize = 1 << 15;

bool ignoreLine(const std::string& line) {
  return line.empty() || line[0] == '%' || line[0] == '#' ||
         std::isalpha(line[0]) != 0 ||
         line.find_first_not_of(' ') == std::string::npos;
}

}  // namespace

Acts::InterpolatedBFieldMap<
    Acts::Grid<Acts::Vector2, Acts::Axis<Acts::AxisType::Equidistant>,
               Acts::Axis<Acts::AxisType::Equidistant>>>
Acts::makeMagneticFieldMapRzFromText(
    const std::function<std::size_t(std::array<std::size_t, 2> binsRZ,
                                    std::array<std::size_t, 2> nBinsRZ)>&
        localToGlobalBin,
    const std::string& fieldMapFile, double lengthUnit, double BFieldUnit,
    bool firstQuadrant, const std::string& delimiter) {
  /// [1] Read in field map file
  // Grid position points in r and z
  std::vector<double> rPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Vector2> bField;
  // reserve estimated size
  rPos.reserve(kDefaultSize);
  zPos.reserve(kDefaultSize);
  bField.reserve(kDefaultSize);
  // [1] Read in file and fill values
  std::ifstream map_file(fieldMapFile.c_str(), std::ios::in);
  std::string line;
  double r = 0., z = 0.;
  double br = 0., bz = 0.;
  while (std::getline(map_file, line)) {
    if (ignoreLine(line)) {
      continue;
    }

    std::istringstream tmp(line);
    if (delimiter.empty()) {
      tmp >> r >> z >> br >> bz;
    } else {
      std::vector<std::string> tokens;
      boost::split(tokens, line, boost::is_any_of(delimiter));
      if (tokens.size() == 4) {
        r = std::stod(tokens[0]);
        z = std::stod(tokens[1]);
        br = std::stod(tokens[2]);
        bz = std::stod(tokens[3]);
      } else {
        throw std::runtime_error("Invalid field map format");
      }
    }
    rPos.push_back(r);
    zPos.push_back(z);
    bField.push_back(Vector2(br, bz));
  }
  map_file.close();
  rPos.shrink_to_fit();
  zPos.shrink_to_fit();
  bField.shrink_to_fit();

  /// [2] use helper function in core
  return fieldMapRZ(localToGlobalBin, rPos, zPos, bField, lengthUnit,
                    BFieldUnit, firstQuadrant);
}

Acts::InterpolatedBFieldMap<
    Acts::Grid<Acts::Vector3, Acts::Axis<Acts::AxisType::Equidistant>,
               Acts::Axis<Acts::AxisType::Equidistant>,
               Acts::Axis<Acts::AxisType::Equidistant>>>
Acts::makeMagneticFieldMapXyzFromText(
    const std::function<std::size_t(std::array<std::size_t, 3> binsXYZ,
                                    std::array<std::size_t, 3> nBinsXYZ)>&
        localToGlobalBin,
    const std::string& fieldMapFile, double lengthUnit, double BFieldUnit,
    bool firstOctant, const std::string& delimiter) {
  /// [1] Read in field map file
  // Grid position points in x, y and z
  std::vector<double> xPos;
  std::vector<double> yPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Vector3> bField;
  // reserve estimated size
  xPos.reserve(kDefaultSize);
  yPos.reserve(kDefaultSize);
  zPos.reserve(kDefaultSize);
  bField.reserve(kDefaultSize);
  // [1] Read in file and fill values
  std::ifstream map_file(fieldMapFile.c_str(), std::ios::in);
  std::string line;
  double x = 0., y = 0., z = 0.;
  double bx = 0., by = 0., bz = 0.;
  while (std::getline(map_file, line)) {
    if (ignoreLine(line)) {
      continue;
    }

    std::istringstream tmp(line);
    if (delimiter.empty()) {
      tmp >> x >> y >> z >> bx >> by >> bz;
    } else {
      std::vector<std::string> tokens;
      boost::split(tokens, line, boost::is_any_of(delimiter));
      if (tokens.size() == 6) {
        x = std::stod(tokens[0]);
        y = std::stod(tokens[1]);
        z = std::stod(tokens[2]);
        bx = std::stod(tokens[3]);
        by = std::stod(tokens[4]);
        bz = std::stod(tokens[5]);
      } else {
        throw std::runtime_error("Invalid field map format");
      }
    }
    xPos.push_back(x);
    yPos.push_back(y);
    zPos.push_back(z);
    bField.push_back(Vector3(bx, by, bz));
  }
  map_file.close();
  xPos.shrink_to_fit();
  yPos.shrink_to_fit();
  zPos.shrink_to_fit();
  bField.shrink_to_fit();

  return fieldMapXYZ(localToGlobalBin, xPos, yPos, zPos, bField, lengthUnit,
                     BFieldUnit, firstOctant);
}
