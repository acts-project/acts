// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include <memory>
#include <string>

#include <G4UserSteppingAction.hh>

class G4VPhysicalVolume;
class G4Step;

namespace Acts {
class Surface;
}

namespace ActsExamples::Geant4 {

///   @brief G4 user stepping action to record the step of a particle through the
///          the G4 volumes. The information is dumped in the `ActsFatras::Hit`
///  format containing the 4-momentum before & after the pre-step, the
///  particle's
///          barcode & pdgId as well as the GeometryIdentifer from the
///          associated
///  Acts::Surface provided by an external translation of the geometry to the
///  Acts::TrackingGeometry. Using the configuration object, the action can be
///  configured to record the steps of primary / secondary (un)charged
///  particles. By default, only steps in volumes that are marked by
///  `SensitiveSurfaceMapper` as sensitive are recorded, but the
///  `SensitiveSteppingAction` may also record every single G4 step.
class SensitiveSteppingAction : public G4UserSteppingAction {
 public:
  /// Configuration of the Stepping action
  struct Config {
    std::shared_ptr<EventStore> eventStore;

    /// @brief Record charged particles
    bool charged = true;
    /// @brief Record neutral particles
    bool neutral = false;
    /// @brief Record the particles produced in the primary interactions
    bool primary = true;
    /// @brief Record secondary particles produced by the traversing
    ///        primary particles
    bool secondary = true;
    /// @brief Record every single step made by the particles. Otherwise,
    ///        steps in volumes where the G4VPhysicalVolume's name has the
    ///        SensitiveSurfaceMapper::mappingPrefix prepended
    bool stepLogging = false;
  };

  /// Construct the stepping action
  ///
  /// @param cfg the configuration struct
  /// @param logger the ACTS logging instance
  explicit SensitiveSteppingAction(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "SensitiveSteppingAction", Acts::Logging::INFO));
  ~SensitiveSteppingAction() override = default;

  /// @brief Interface Method doing the step and records the data
  /// @param step is the Geant4 step of the particle
  void UserSteppingAction(const G4Step* step) override;

  /// Set the multimap that correlates G4VPhysicalVolumes to Acts::Surfaces
  ///
  /// @param surfaceMapping the multimap of physical volumes to surfaces
  using VolumeToSurfAssocMap_t = SensitiveSurfaceMapper::VolumeToSurfAssocMap_t;
  void assignSurfaceMapping(const VolumeToSurfAssocMap_t& surfaceMapping) {
    m_surfaceMapping = surfaceMapping;
  }

 protected:
  Config m_cfg;

 private:
  /// Private access method to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// Private access method to the event store
  EventStore& eventStore() const { return *m_cfg.eventStore; }

  /// The looging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  VolumeToSurfAssocMap_t m_surfaceMapping{};
};

}  // namespace ActsExamples::Geant4
