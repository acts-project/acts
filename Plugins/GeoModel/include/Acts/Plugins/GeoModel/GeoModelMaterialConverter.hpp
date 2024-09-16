// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Material/Material.hpp"

class GeoMaterial;

namespace Acts::GeoModel {

/// @brief Convert GeoMaterial to Acts::Material
///
/// @param gm The GeoMaterial to be converted
/// @return the Acts::Material
Material geoMaterialConverter(const GeoMaterial& gm,
                              bool useMolarDensity = true);

}  // namespace Acts::GeoModel
