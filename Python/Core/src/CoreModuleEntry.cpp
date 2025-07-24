// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "ActsPython/Utilities/Context.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

/// This adds the core module entries to the Python module
/// There is a certain order necessary as py::class_ definitions
/// need to be registered before they can be used in other modules.
namespace ActsPython {

void addDefinitions(Context& ctx);
void addMagneticField(Context& ctx);
void addUtilities(Context& ctx);
void addVisualization(Context& ctx);

void addMaterial(Context& ctx);
void addSurfaces(Context& ctx);
void addGeometry(Context& ctx);
void addGeometryGen1(Context& ctx);
void addGeometryGen2(Context& ctx);
void addGeometryGen3(Context& ctx);
void addNavigation(Context& ctx);
void addPropagation(Context& ctx);
void addSeeding(Context& ctxt);
void addTrackFinding(Context& ctx);
}  // namespace ActsPython

namespace py = pybind11;

PYBIND11_MODULE(ActsPythonBindings, m) {
  using namespace ActsPython;

  Context ctx;
  ctx.modules["main"] = m;
  m.doc() = "Acts";

  m.attr("__version__") =
      std::tuple{Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

  {
    auto mv = m.def_submodule("version");

    mv.attr("major") = Acts::VersionMajor;
    mv.attr("minor") = Acts::VersionMinor;
    mv.attr("patch") = Acts::VersionPatch;

    mv.attr("commit_hash") = Acts::CommitHash;
    mv.attr("commit_hash_short") = Acts::CommitHashShort;
  }

  addDefinitions(ctx);
  addMagneticField(ctx);
  addMaterial(ctx);
  addUtilities(ctx);
  addVisualization(ctx);

  addSurfaces(ctx);
  addGeometry(ctx);
  addGeometryGen1(ctx);
  addGeometryGen2(ctx);
  addGeometryGen3(ctx);
  addNavigation(ctx);
  addPropagation(ctx);
  addSeeding(ctx);
  addTrackFinding(ctx);
}
