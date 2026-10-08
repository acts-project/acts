// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// A fitted synthetic event configuration, as a file. This is where the numbers
/// a fit against a full simulation produces belong: they are measurements of
/// one detector and one layout of it, not constants.
///
/// Every field is written and every field is required on reading. A
/// configuration is fitted as a whole, so a missing one silently taking a
/// default would retune the rest of it -- which is exactly what the shipped
/// presets used to spell out every field in C++ to avoid.

#include "ActsFatras/Synthetic/EventConfig.hpp"

#include <filesystem>

#include <nlohmann/json.hpp>

namespace ActsFatras::Synthetic {

/// @addtogroup fatras_synthetic_json
/// @{

/// Convert EnergyLossModel to JSON
/// @param j Destination JSON object
/// @param model Source EnergyLossModel to convert
void to_json(nlohmann::json& j, const EnergyLossModel& model);
/// Convert JSON to EnergyLossModel
/// @param j Source JSON object
/// @param model Destination EnergyLossModel to populate
void from_json(const nlohmann::json& j, EnergyLossModel& model);

/// Convert GenerationConfig to JSON
/// @param j Destination JSON object
/// @param config Source GenerationConfig to convert
void to_json(nlohmann::json& j, const GenerationConfig& config);
/// Convert JSON to GenerationConfig
/// @param j Source JSON object
/// @param config Destination GenerationConfig to populate
void from_json(const nlohmann::json& j, GenerationConfig& config);

/// Convert PropagationConfig to JSON
/// @param j Destination JSON object
/// @param config Source PropagationConfig to convert
void to_json(nlohmann::json& j, const PropagationConfig& config);
/// Convert JSON to PropagationConfig
/// @param j Source JSON object
/// @param config Destination PropagationConfig to populate
void from_json(const nlohmann::json& j, PropagationConfig& config);

/// Convert MaterialConfig to JSON
/// @param j Destination JSON object
/// @param config Source MaterialConfig to convert
void to_json(nlohmann::json& j, const MaterialConfig& config);
/// Convert JSON to MaterialConfig
/// @param j Source JSON object
/// @param config Destination MaterialConfig to populate
void from_json(const nlohmann::json& j, MaterialConfig& config);

/// Convert MeasurementConfig to JSON
/// @param j Destination JSON object
/// @param config Source MeasurementConfig to convert
void to_json(nlohmann::json& j, const MeasurementConfig& config);
/// Convert JSON to MeasurementConfig
/// @param j Source JSON object
/// @param config Destination MeasurementConfig to populate
void from_json(const nlohmann::json& j, MeasurementConfig& config);

/// Convert SecondarySamplingConfig to JSON
/// @param j Destination JSON object
/// @param config Source SecondarySamplingConfig to convert
void to_json(nlohmann::json& j, const SecondarySamplingConfig& config);
/// Convert JSON to SecondarySamplingConfig
/// @param j Source JSON object
/// @param config Destination SecondarySamplingConfig to populate
void from_json(const nlohmann::json& j, SecondarySamplingConfig& config);

/// Convert SecondaryConfig to JSON
/// @param j Destination JSON object
/// @param config Source SecondaryConfig to convert
void to_json(nlohmann::json& j, const SecondaryConfig& config);
/// Convert JSON to SecondaryConfig
/// @param j Source JSON object
/// @param config Destination SecondaryConfig to populate
void from_json(const nlohmann::json& j, SecondaryConfig& config);

/// Convert SimulationConfig to JSON
/// @param j Destination JSON object
/// @param config Source SimulationConfig to convert
void to_json(nlohmann::json& j, const SimulationConfig& config);
/// Convert JSON to SimulationConfig
/// @param j Source JSON object
/// @param config Destination SimulationConfig to populate
void from_json(const nlohmann::json& j, SimulationConfig& config);

/// Convert EventConfig to JSON
/// @param j Destination JSON object
/// @param config Source EventConfig to convert
void to_json(nlohmann::json& j, const EventConfig& config);
/// Convert JSON to EventConfig
/// @param j Source JSON object
/// @param config Destination EventConfig to populate
void from_json(const nlohmann::json& j, EventConfig& config);

/// Read an event configuration from file.
/// @param path the file to read
/// @return the configuration
/// @throws std::runtime_error if the file is not there or is not one of these
/// @throws nlohmann::json::out_of_range if it leaves a field out
EventConfig readEventConfig(const std::filesystem::path& path);

/// Write an event configuration to file.
/// @param path the file to write
/// @param config the configuration to write
void writeEventConfig(const std::filesystem::path& path,
                      const EventConfig& config);

/// @}

}  // namespace ActsFatras::Synthetic
