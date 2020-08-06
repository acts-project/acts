// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/IVisualization3D.hpp"

void Acts::IVisualization3D::vertex(const Vector3F& vtx, ColorRGB color) {
  Vector3D vtxd = vtx.template cast<double>();
  vertex(vtxd, color);
}

void Acts::IVisualization3D::face(const std::vector<Vector3F>& vtxs,
                                ColorRGB color) {
  std::vector<Vector3D> vtxsd;
  std::transform(vtxs.begin(), vtxs.end(), std::back_inserter(vtxsd),
                 [](auto& v) { return v.template cast<double>(); });
  face(vtxsd, color);
}

void Acts::IVisualization3D::line(const Vector3F& a, const Vector3F& b,
                                ColorRGB color) {
  Vector3D ad = a.template cast<double>();
  Vector3D bd = b.template cast<double>();
  line(ad, bd, color);
}

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