// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// A synthetic detector description and its material, as files.
///
/// The two are separate on purpose: a description says where a detector's
/// layers are and a decoration says what they are made of, they are generated
/// by different scripts from different inputs, and either can be replaced
/// without touching the other. `ActsFatras::Synthetic::decorate` puts them
/// together.
///
/// The JSON follows the structs one for one, so that a file can be read as the
/// description it is. Two conventions are worth knowing:
///
/// - An infinite bound is written as the string `"inf"`, `nlohmann::json`
///   having no number for it. A service surface spanning the whole detector is
///   the usual case.
/// - The per-band material numbers may each be a single value or one per band,
///   so a surface of one composition banded only in how much of it there is
///   stays as short to read as it is to state.

#include "ActsFatras/Synthetic/DetectorLayout.hpp"

#include <filesystem>

#include <nlohmann/json.hpp>

namespace ActsFatras::Synthetic {

/// @defgroup fatras_synthetic_json Synthetic detector description as JSON
/// @{

/// Convert SurfaceShape to JSON
/// @param j Destination JSON object
/// @param shape Source SurfaceShape to convert
void to_json(nlohmann::json& j, const SurfaceShape& shape);
/// Convert JSON to SurfaceShape
/// @param j Source JSON object
/// @param shape Destination SurfaceShape to populate
void from_json(const nlohmann::json& j, SurfaceShape& shape);

/// Convert EndcapPlacement to JSON
/// @param j Destination JSON object
/// @param placement Source EndcapPlacement to convert
void to_json(nlohmann::json& j, const EndcapPlacement& placement);
/// Convert JSON to EndcapPlacement
/// @param j Source JSON object
/// @param placement Destination EndcapPlacement to populate
void from_json(const nlohmann::json& j, EndcapPlacement& placement);

/// Convert LayerKind to JSON
/// @param j Destination JSON object
/// @param kind Source LayerKind to convert
void to_json(nlohmann::json& j, const LayerKind& kind);
/// Convert JSON to LayerKind
/// @param j Source JSON object
/// @param kind Destination LayerKind to populate
void from_json(const nlohmann::json& j, LayerKind& kind);

/// Convert RingBounds to JSON
/// @param j Destination JSON object
/// @param ring Source RingBounds to convert
void to_json(nlohmann::json& j, const RingBounds& ring);
/// Convert JSON to RingBounds
/// @param j Source JSON object
/// @param ring Destination RingBounds to populate
void from_json(const nlohmann::json& j, RingBounds& ring);

/// Convert SurfaceMaterial to JSON
/// @param j Destination JSON object
/// @param material Source SurfaceMaterial to convert
void to_json(nlohmann::json& j, const SurfaceMaterial& material);
/// Convert JSON to SurfaceMaterial
/// @param j Source JSON object
/// @param material Destination SurfaceMaterial to populate
void from_json(const nlohmann::json& j, SurfaceMaterial& material);

/// Convert CylinderDescription to JSON
/// @param j Destination JSON object
/// @param cylinder Source CylinderDescription to convert
void to_json(nlohmann::json& j, const CylinderDescription& cylinder);
/// Convert JSON to CylinderDescription
/// @param j Source JSON object
/// @param cylinder Destination CylinderDescription to populate
void from_json(const nlohmann::json& j, CylinderDescription& cylinder);

/// Convert DiscDescription to JSON
/// @param j Destination JSON object
/// @param disc Source DiscDescription to convert
void to_json(nlohmann::json& j, const DiscDescription& disc);
/// Convert JSON to DiscDescription
/// @param j Source JSON object
/// @param disc Destination DiscDescription to populate
void from_json(const nlohmann::json& j, DiscDescription& disc);

/// Convert PassiveSurfaceDescription to JSON
/// @param j Destination JSON object
/// @param passive Source PassiveSurfaceDescription to convert
void to_json(nlohmann::json& j, const PassiveSurfaceDescription& passive);
/// Convert JSON to PassiveSurfaceDescription
/// @param j Source JSON object
/// @param passive Destination PassiveSurfaceDescription to populate
void from_json(const nlohmann::json& j, PassiveSurfaceDescription& passive);

/// Convert BarrelDescription to JSON
/// @param j Destination JSON object
/// @param barrel Source BarrelDescription to convert
void to_json(nlohmann::json& j, const BarrelDescription& barrel);
/// Convert JSON to BarrelDescription
/// @param j Source JSON object
/// @param barrel Destination BarrelDescription to populate
void from_json(const nlohmann::json& j, BarrelDescription& barrel);

/// Convert EndcapDescription to JSON
/// @param j Destination JSON object
/// @param endcap Source EndcapDescription to convert
void to_json(nlohmann::json& j, const EndcapDescription& endcap);
/// Convert JSON to EndcapDescription
/// @param j Source JSON object
/// @param endcap Destination EndcapDescription to populate
void from_json(const nlohmann::json& j, EndcapDescription& endcap);

/// Convert SubsystemDescription to JSON
/// @param j Destination JSON object
/// @param subsystem Source SubsystemDescription to convert
void to_json(nlohmann::json& j, const SubsystemDescription& subsystem);
/// Convert JSON to SubsystemDescription
/// @param j Source JSON object
/// @param subsystem Destination SubsystemDescription to populate
void from_json(const nlohmann::json& j, SubsystemDescription& subsystem);

/// Convert DetectorDescription to JSON
/// @param j Destination JSON object
/// @param description Source DetectorDescription to convert
void to_json(nlohmann::json& j, const DetectorDescription& description);
/// Convert JSON to DetectorDescription
/// @param j Source JSON object
/// @param description Destination DetectorDescription to populate
void from_json(const nlohmann::json& j, DetectorDescription& description);

/// Convert LayerId to JSON
/// @param j Destination JSON object
/// @param layer Source LayerId to convert
void to_json(nlohmann::json& j, const LayerId& layer);
/// Convert JSON to LayerId
/// @param j Source JSON object
/// @param layer Destination LayerId to populate
void from_json(const nlohmann::json& j, LayerId& layer);

/// Convert MaterialEntry to JSON
/// @param j Destination JSON object
/// @param entry Source MaterialEntry to convert
void to_json(nlohmann::json& j, const MaterialEntry& entry);
/// Convert JSON to MaterialEntry
/// @param j Source JSON object
/// @param entry Destination MaterialEntry to populate
void from_json(const nlohmann::json& j, MaterialEntry& entry);

/// Read a detector description from file.
/// @param path the file to read
/// @return the description
/// @throws std::runtime_error if the file is not there or is not one of these
DetectorDescription readDetectorDescription(const std::filesystem::path& path);

/// Write a detector description to file, material included where its layers
/// carry any.
/// @param path the file to write
/// @param description the description to write
void writeDetectorDescription(const std::filesystem::path& path,
                              const DetectorDescription& description);

/// Read the material of a detector from file.
/// @param path the file to read
/// @return the decoration, to be applied with `decorate`
/// @throws std::runtime_error if the file is not there or is not one of these
MaterialDecoration readMaterialDecoration(const std::filesystem::path& path);

/// Write the material of a detector to file.
/// @param path the file to write
/// @param decoration the material to write
void writeMaterialDecoration(const std::filesystem::path& path,
                             const MaterialDecoration& decoration);

/// @}

}  // namespace ActsFatras::Synthetic
