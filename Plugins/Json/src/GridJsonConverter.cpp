// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/GridJsonConverter.hpp"

#include "Acts/Utilities/IAxis.hpp"

nlohmann::json Acts::AxisJsonConverter::toJson(const IAxis& ia) {
  nlohmann::json jAxis;
  jAxis["boundary_type"] = ia.getBoundaryType();
  // type, range, bins or boundaries
  if (ia.isEquidistant()) {
    jAxis["type"] = detail::AxisType::Equidistant;
    jAxis["range"] = std::array<ActsScalar, 2u>({ia.getMin(), ia.getMax()});
    jAxis["bins"] = ia.getNBins();
  } else {
    jAxis["type"] = detail::AxisType::Variable;
    jAxis["boundaries"] = ia.getBinEdges();
  }
  return jAxis;
}
