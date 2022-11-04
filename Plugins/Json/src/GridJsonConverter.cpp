// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/GridJsonConverter.hpp"

#include "Acts/Utilities/IAxis.hpp"

void Acts::to_json(nlohmann::json& j, const IAxis& ia) {
  j["type"] = int(ia.getBoundaryType());
  if (ia.isEquidistant()) {
    j["spacing"] = "equidistant";
    j["range"] = std::array<ActsScalar, 2u>({ia.getMin(), ia.getMax()});
    j["bins"] = ia.getNBins();
  } else {
    j["spacing"] = "variable";
    j["boundaries"] = ia.getBinEdges();
  }
}
