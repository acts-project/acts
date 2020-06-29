// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>

#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

/// Common obj writing options
///
/// @tparam aopt_t Type of the options object (from BOOST)
///
/// @param opt The options object, where string based options are attached
template <typename aopt_t>
void addObjWriterOptions(aopt_t& opt) {
  opt.add_options()("obj-precision", po::value<int>()->default_value(6),
                    "Floating number output precission.")(
      "obj-scalor", po::value<double>()->default_value(1.),
      "Optional scaling from Acts units to ouput units.")(
      "obj-container-view",
      po::value<read_series>()->multitoken()->default_value(
          {0, 220, 220, 220, 0}),
      "View configuration of container volumes (vis/novis, r, g, b, trimesh).")(
      "obj-volume-view",
      po::value<read_series>()->multitoken()->default_value(
          {1, 220, 220, 0, 0}),
      "View configuration of navigation volumes (vis/novis, r, g, b, "
      "trimesh).")(
      "obj-layer-view",
      po::value<read_series>()->multitoken()->default_value(
          {1, 100, 180, 240, 0}),
      "View configuration of layer structures (vis/novis, r, g, b, trimesh).")(
      "obj-sensitive-view",
      po::value<read_series>()->multitoken()->default_value(
          {1, 0, 180, 240, 0}),
      "View configuration of sensitive surfaces (vis/novis, r, g, b, "
      "trimesh).")("obj-passive-view",
                   po::value<read_series>()->multitoken()->default_value(
                       {1, 240, 280, 0, 0}),
                   "View configuration of sensitive surfaces (vis/novis, r, g, "
                   "b, trimesh).")(
      "obj-grid-view",
      po::value<read_series>()->multitoken()->default_value({1, 220, 0, 0, 0}),
      "View configuration of grid structures (vis/novis, r, g, b, trimesh).")(
      "obj-grid-offset", po::value<double>()->default_value(0.),
      "View offset of grid values.")("obj-grid-thickness",
                                     po::value<double>()->default_value(0.5),
                                     "Thickness of grid objects.");
}

/// read the evgen options and return a Config file
template <class amap_t>
FW::Obj::ObjTrackingGeometryWriter::Config readObjTrackingGeometryWriterConfig(
    const amap_t& vm, const std::string& name,
    Acts::Logging::Level loglevel = Acts::Logging::INFO) {
  FW::Obj::ObjTrackingGeometryWriter::Config objTgConfig(name, loglevel);

  objTgConfig.outputPrecision = vm["obj-precision"].template as<int>();
  objTgConfig.outputScalor = vm["obj-scalor"].template as<double>();

  auto setView = [&](const std::string& vname,
                     Acts::ViewConfig& viewCfg) -> void {
    read_series cview = vm[vname].template as<read_series>();
    if (not cview.empty()) {
      if (cview[0] == 0) {
        viewCfg.visible = false;
      } else if (cview.size() > 3) {
        viewCfg.color = {cview[1], cview[2], cview[3]};
        if (cview.size() > 4 and cview[4] != 0) {
          viewCfg.triangulate = true;
        }
      }
    }
  };

  setView("obj-container-view", objTgConfig.containerView);
  setView("obj-volume-view", objTgConfig.containerView);
  setView("obj-sensitive-view", objTgConfig.sensitiveView);
  setView("obj-passive-view", objTgConfig.passiveView);
  setView("obj-grid-view", objTgConfig.gridView);
  objTgConfig.gridView.offset = vm["obj-grid-offset"].template as<double>();
  objTgConfig.gridView.lineThickness =
      vm["obj-grid-thickness"].template as<double>();

  return objTgConfig;
}

}  // namespace Options
}  // namespace FW
