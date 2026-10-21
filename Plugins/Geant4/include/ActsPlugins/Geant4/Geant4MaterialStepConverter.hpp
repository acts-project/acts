// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialInteraction.hpp"

#include <optional>
#include <string>
#include <vector>

class G4Step;
class G4Track;

namespace ActsPlugins {
/// @addtogroup geant4_plugin
/// @{

/// @brief Converts Geant4 steps into recorded material tracks
///
/// This is the Geant4 side of a material recording with geantinos: an
/// @ref Acts::RecordedMaterialTrack is started from a Geant4 track with
/// @ref startTrack, and every step of it through a material is added to it
/// with @ref addStep as one @ref Acts::MaterialInteraction. The caller owns
/// the tracks.
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

  /// @brief Start a recorded material track
  ///
  /// @param track the Geant4 track
  ///
  /// @return an empty material track, starting at the vertex position and
  ///         the initial direction of @p track
  Acts::RecordedMaterialTrack startTrack(const G4Track& track) const;

  /// @brief Add the material of a step to a recorded material track
  ///
  /// The material is appended as one interaction, and added to the
  /// accumulated material of the track.
  ///
  /// @param step the Geant4 step
  /// @param track the material track the step belongs to
  ///
  /// @return whether the step was recorded, it is not if it has no material
  ///         or the material is excluded
  bool addStep(const G4Step& step, Acts::RecordedMaterialTrack& track) const;

  /// Readonly access to the configuration
  /// @return the configuration
  const Config& config() const { return m_cfg; }

 private:
  /// Convert a step, nothing if it has no or an excluded material
  std::optional<Acts::MaterialInteraction> materialInteraction(
      const G4Step& step) const;

  Config m_cfg;
};

/// @}
}  // namespace ActsPlugins
