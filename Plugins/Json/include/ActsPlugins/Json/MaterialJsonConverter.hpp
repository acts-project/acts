// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

class Surface;
class ISurfaceMaterial;
class IVolumeMaterial;
class BinUtility;

/// Pointer to a constant volume material object
using volumeMaterialPointer = const Acts::IVolumeMaterial*;
/// Pointer to a constant surface material object
using surfaceMaterialPointer = const Acts::ISurfaceMaterial*;

/// Convert Material to JSON
/// @param j Destination JSON object
/// @param t Source Material to convert
void to_json(nlohmann::json& j, const Material& t);

/// Convert JSON to Material
/// @param j Source JSON object
/// @param t Destination Material to populate
void from_json(const nlohmann::json& j, Material& t);

/// Convert MaterialSlab to JSON
/// @param j Destination JSON object
/// @param t Source MaterialSlab to convert
void to_json(nlohmann::json& j, const MaterialSlab& t);

/// Convert JSON to MaterialSlab
/// @param j Source JSON object
/// @param t Destination MaterialSlab to populate
void from_json(const nlohmann::json& j, MaterialSlab& t);

/// Convert JSON to MaterialSlabMatrix
/// @param j Source JSON object
/// @param t Destination MaterialSlabMatrix to populate
void from_json(const nlohmann::json& j, MaterialSlabMatrix& t);

/// Convert volumeMaterialPointer to JSON
/// @param j Destination JSON object
/// @param material Source volumeMaterialPointer to convert
void to_json(nlohmann::json& j, const volumeMaterialPointer& material);

/// Convert JSON to volumeMaterialPointer
/// @param j Source JSON object
/// @param material Destination volumeMaterialPointer to populate
void from_json(const nlohmann::json& j, volumeMaterialPointer& material);

/// Convert surfaceMaterialPointer to JSON
/// @param j Destination JSON object
/// @param material Source surfaceMaterialPointer to convert
void to_json(nlohmann::json& j, const surfaceMaterialPointer& material);

/// Convert JSON to surfaceMaterialPointer
/// @param j Source JSON object
/// @param material Destination surfaceMaterialPointer to populate
void from_json(const nlohmann::json& j, surfaceMaterialPointer& material);

/// JSON serialization mapping for MappingType enum
// This macro create a conversion for the mapping type enum
NLOHMANN_JSON_SERIALIZE_ENUM(Acts::MappingType,
                             {
                                 {Acts::MappingType::PreMapping, "PreMapping"},
                                 {Acts::MappingType::Default, "Default"},
                                 {Acts::MappingType::PostMapping,
                                  "PostMapping"},
                                 {Acts::MappingType::Sensor, "Sensor"},
                             })

/// @}
}  // namespace Acts
