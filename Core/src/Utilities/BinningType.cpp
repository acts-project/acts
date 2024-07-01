// This file is part of the Acts project.
//
// Copyright (C) 2016-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/BinningType.hpp"

#include <algorithm>
#include <stdexcept>

namespace {
const std::vector<std::string>& binningValueNames() {
  static const std::vector<std::string> _binningValueNames = {
      "binX",    "binY", "binZ",   "binR",  "binPhi",
      "binRPhi", "binH", "binEta", "binMag"};
  return _binningValueNames;
}
}  // namespace

namespace Acts {

const std::vector<BinningValue>& allBinningValues() {
  static const std::vector<BinningValue> s_binningValues = {
      BinningValue::binX, BinningValue::binY,   BinningValue::binZ,
      BinningValue::binR, BinningValue::binPhi, BinningValue::binRPhi,
      BinningValue::binH, BinningValue::binEta, BinningValue::binMag};
  return s_binningValues;
}

BinningValue binningValueFromName(const std::string& name) {
  auto it =
      std::find(binningValueNames().begin(), binningValueNames().end(), name);
  if (it == binningValueNames().end()) {
    throw std::invalid_argument("Unknown binning value name: " + name);
  }
  return static_cast<BinningValue>(
      std::distance(binningValueNames().begin(), it));
}

const std::string& binningValueName(BinningValue bValue) {
  return binningValueNames().at(
      static_cast<std::underlying_type_t<BinningValue>>(bValue));
}

std::ostream& operator<<(std::ostream& os, BinningValue bValue) {
  os << binningValueName(bValue);
  return os;
}

}  // namespace Acts
