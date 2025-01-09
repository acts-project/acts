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

namespace Acts {

namespace {

static const std::vector<std::string> s_axisDirectionNames = {
    "AxisX",    "AxisY",     "AxisZ",   "AxisR",  "AxisPhi",
    "AxisRPhi", "AxisTheta", "AxisEta", "AxisMag"};

static const std::vector<AxisDirection> s_axisDirections = {
    AxisDirection::AxisX,     AxisDirection::AxisY,   AxisDirection::AxisZ,
    AxisDirection::AxisR,     AxisDirection::AxisPhi, AxisDirection::AxisRPhi,
    AxisDirection::AxisTheta, AxisDirection::AxisEta, AxisDirection::AxisMag};

}  // namespace

const std::vector<AxisDirection>& allAxisDirections() {
  return s_axisDirections;
}

AxisDirection axisDirectionFromName(const std::string& name) {
  auto it = std::ranges::find(s_axisDirectionNames, name);
  if (it == s_axisDirectionNames.end()) {
    throw std::invalid_argument("Unknown Axisning value name: " + name);
  }
  return static_cast<AxisDirection>(
      std::distance(s_axisDirectionNames.begin(), it));
}

const std::string& axisDirectionName(AxisDirection aDir) {
  return s_axisDirectionNames.at(
      static_cast<std::underlying_type_t<AxisDirection>>(aDir));
}

std::ostream& operator<<(std::ostream& os, AxisDirection aDir) {
  os << axisDirectionName(aDir);
  return os;
}

}  // namespace Acts
