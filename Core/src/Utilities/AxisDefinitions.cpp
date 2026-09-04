// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/AxisDefinitions.hpp"

#include <algorithm>
#include <ostream>
#include <stdexcept>

using namespace Acts;

namespace {

// Legacy binning value names
//
// NOTE: this should be removed once the BinUtility is removed
const std::vector<std::string> s_legacyBinningValueNames = {
    "binX",    "binY", "binZ",   "binR",  "binPhi",
    "binRPhi", "binH", "binEta", "binMag"};
// end of legacy binning values

const std::vector<std::string> s_axisDirectionNames = {
    "AxisX",    "AxisY",     "AxisZ",   "AxisR",  "AxisPhi",
    "AxisRPhi", "AxisTheta", "AxisEta", "AxisMag"};

const std::vector<AxisDirection> s_axisDirections = {
    AxisDirection::AxisX,     AxisDirection::AxisY,   AxisDirection::AxisZ,
    AxisDirection::AxisR,     AxisDirection::AxisPhi, AxisDirection::AxisRPhi,
    AxisDirection::AxisTheta, AxisDirection::AxisEta, AxisDirection::AxisMag};

}  // namespace

const std::vector<AxisDirection>& Acts::allAxisDirections() {
  return s_axisDirections;
}

AxisDirection Acts::axisDirectionFromName(const std::string& name) {
  auto it = std::ranges::find(s_axisDirectionNames, name);
  if (it == s_axisDirectionNames.end()) {
    // Legacy binning check - this should be removed once BinUtility is gone
    it = std::ranges::find(s_legacyBinningValueNames, name);
    if (it == s_legacyBinningValueNames.end()) {
      // both legacy and current failed
      throw std::invalid_argument("Unknown AxisDirection value name: " + name);
    }
    return static_cast<AxisDirection>(
        std::distance(s_legacyBinningValueNames.begin(), it));
  }
  return static_cast<AxisDirection>(
      std::distance(s_axisDirectionNames.begin(), it));
}

const std::string& Acts::axisDirectionName(AxisDirection aDir) {
  return s_axisDirectionNames.at(
      static_cast<std::underlying_type_t<AxisDirection>>(aDir));
}

std::string Acts::axesDirectionName(const std::vector<AxisDirection>& manyDir) {
  std::string name = "{";
  for (const auto& dir : manyDir) {
    name += axisDirectionName(dir);
    if (&dir != &manyDir.back()) {
      name += ", ";
    }
  }
  name += "}";
  return name;
}

std::ostream& Acts::operator<<(std::ostream& os, AxisDirection aDir) {
  os << axisDirectionName(aDir);
  return os;
}

std::ostream& Acts::operator<<(std::ostream& os,
                               const std::vector<AxisDirection>& manyDir) {
  os << axesDirectionName(manyDir);
  return os;
}

std::ostream& Acts::operator<<(std::ostream& os, AxisBoundaryType bdt) {
  using enum AxisBoundaryType;
  switch (bdt) {
    case Open:
      os << "Open";
      break;
    case Bound:
      os << "Bound";
      break;
    case Closed:
      os << "Closed";
      break;
  }
  return os;
}

std::ostream& Acts::operator<<(std::ostream& os, AxisType type) {
  switch (type) {
    case AxisType::Equidistant:
      os << "Equidistant";
      break;
    case AxisType::Variable:
      os << "Variable";
      break;
  }
  return os;
}
