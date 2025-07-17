// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Module/Entries.hpp"
#include "ActsPython/Utilities/Context.hpp"

namespace ActsPython {
void addDefinitions(Context& ctx);
void addMagneticField(Context& ctx);
void addMaterial(Context& ctx);
void addSurfaces(Context& ctx);
void addGeometry(Context& ctx);
void addGeometryGen1(Context& ctx);
void addGeometryGen2(Context& ctx);
void addGeometryGen3(Context& ctx);
void addNavigation(Context& ctx);
void addUtilities(Context& ctx);
void addSeeding(Context& ctxt);
void addTrackFinding(Context& ctx);
}  // namespace ActsPython

void ActsPython::addCoreModule(Context& ctx) {
  addDefinitions(ctx);
  addMagneticField(ctx);
  addMaterial(ctx);
  addSurfaces(ctx);
  addGeometry(ctx);
  addGeometryGen1(ctx);
  addGeometryGen2(ctx);
  addGeometryGen3(ctx);
  addNavigation(ctx);
  addUtilities(ctx);
  addSeeding(ctx);
  addTrackFinding(ctx);
}
