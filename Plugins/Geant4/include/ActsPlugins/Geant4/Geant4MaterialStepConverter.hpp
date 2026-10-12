// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/MaterialInteraction.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

class G4Step;
class G4Track;

namespace ActsPlugins {
/// @addtogroup geant4_plugin
/// @{

/// @brief Converts Geant4 steps into recorded material interactions
///
/// This is the per-step part of a material recording with geantinos: every
/// step through a material becomes one @c Acts::MaterialInteraction, which
/// are collected per track into an @c Acts::RecordedMaterialTrack.
class Geant4MaterialStepConverter {
 public:
  /// Nested configuration struct
  struct Config {
    /// Names of materials whose steps are not recorded
    std::vector<std::string> excludeMaterials = {};
    /// Whether to record the element composition of the stepped material
    bool recordElementFractions = false;
  };

  /// Constructor
  ///
  /// @param cfg the configuration
  explicit Geant4MaterialStepConverter(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Convert a single step
  ///
  /// @param step the Geant4 step
  ///
  /// @return the material interaction, or nothing if the step has no
  ///         material or the material is excluded
  std::optional<Acts::MaterialInteraction> convert(const G4Step& step) const;

  /// @brief The start of a recorded material track
  ///
  /// @param track the Geant4 track
  ///
  /// @return the vertex position and the initial direction of the track
  static std::pair<Acts::Vector3, Acts::Vector3> trackStart(
      const G4Track& track);

  /// Readonly access to the configuration
  /// @return the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins
