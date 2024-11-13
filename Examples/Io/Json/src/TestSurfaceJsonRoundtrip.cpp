// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Json/JsonSurfacesReader.hpp"
#include <Acts/Plugins/Json/SurfaceJsonConverter.hpp>

#include <fstream>

#include <nlohmann/json.hpp>

void writeSurfaces(
    std::ostream &os,
    const std::vector<std::shared_ptr<Acts::Surface>> &jSurfaces) {
  nlohmann::json j = nlohmann::json::array();

  for (const auto &jSurface : jSurfaces) {
    j.push_back(
        Acts::SurfaceJsonConverter::toJson(Acts::GeometryContext{}, *jSurface));
  }

  os << j;
}

int main(int argc, char **argv) {
  std::vector<std::string> args(argv, argv + argc);

  auto surfaces =
      ActsExamples::JsonSurfacesReader::readVector({args.at(1), {}});

  std::ofstream outfile("output_surfaces.json");
  writeSurfaces(outfile, surfaces);
}
