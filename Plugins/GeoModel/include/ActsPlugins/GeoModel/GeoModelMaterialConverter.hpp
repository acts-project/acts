// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Material/Material.hpp"

class GeoMaterial;

namespace ActsPlugins::GeoModel {

/// @addtogroup geomodel_plugin
/// @{

/// @brief Convert GeoMaterial to Acts::Material
///
/// @param gm The GeoMaterial to be converted
/// @param useMolarDensity Flag to indicate whether to use molar density
/// @return the Acts::Material
Acts::Material geoMaterialConverter(const GeoMaterial& gm,
                                    bool useMolarDensity = true);

/// @}

}  // namespace ActsPlugins::GeoModel
