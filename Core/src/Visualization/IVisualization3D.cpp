// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/IVisualization3D.hpp"

bool Acts::IVisualization3D::hasExtension(const std::string& path) const {
  return (path.find(".") != std::string::npos);
}

void Acts::IVisualization3D::replaceExtension(std::string& path,
                                              const std::string& suffix) const {
  auto ppoint = path.find_last_of(".");
  if (ppoint != std::string::npos) {
    path.replace(ppoint, path.length(), suffix);
  } else {
    path += suffix;
  }
}
