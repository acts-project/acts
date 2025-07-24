// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Utilities/Context.hpp"

#include <pybind11/pybind11.h>

/// This adds the plugins module entries to the Python module
/// There is a certain order necessary as py::class_ definitions
/// need to be registered before they can be used in other modules.
namespace ActsPython {

void addCovfie(Context& ctx);
void addDetray(Context& ctx);
void addObj(Context& ctx);
void addJson(Context& ctx);
void addSvg(Context& ctx);
void addFpeMonitoring(Context& ctx);
// Geometry extensions
void addGeoModel(Context& ctx);
void addTGeo(Context& ctx);

}  // namespace ActsPython

PYBIND11_MODULE(ActsPluginsPythonBindings, p) {
  using namespace ActsPython;
  Context ctx;
  ctx.modules["plugins"] = p;

  addCovfie(ctx);
  addDetray(ctx);
  addObj(ctx);
  addJson(ctx);
  addSvg(ctx);
  addFpeMonitoring(ctx);
  // Geometry extensions
  addGeoModel(ctx);
  addTGeo(ctx);
}
