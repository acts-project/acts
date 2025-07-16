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
    void addUtilities(Context& ctx);
}

void ActsPython::addCoreModule(Context& ctx) {
  addDefinitions(ctx);
  addMagneticField(ctx);
  addMaterial(ctx);
  addSurfaces(ctx);
  addGeometry(ctx);
  addUtilities(ctx);
}
