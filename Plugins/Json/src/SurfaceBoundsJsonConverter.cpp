// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/SurfaceBoundsJsonConverter.hpp"

#include "Acts/Surfaces/SurfaceBounds.hpp"

void Acts::to_json(nlohmann::json& j, const Acts::SurfaceBounds& bounds) {
  j["type"] = boundTypes[bounds.type()];
  j["values"] = bounds.values();
}
