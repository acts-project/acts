// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

class G4VPhysicalVolume;

namespace Acts {
class Surface;
}

namespace ActsExamples::Geant4 {

///   @brief Interface class to query sensitive surfaces from the TrackingGeometry based on the
///          external global position. It serves  the `SensitiveSurfaceMapper`
///          to
///  map the Acts::Surfaces to the G4PhysicalVolumes. Based on the global
///  position of the G4PhyicalVolume or on the position of the volume bound
///  vertices, the tracking geometry is queried to provide the candidate
///  sensitive surfaces. Due to the concurrent TrackingGeometry schemes
///  currently in use, concrete implementations need to follow
struct SensitiveCandidatesBase {
  /// Get the sensitive surfaces for a given position
  ///
  /// @param gctx the geometry context
  /// @param position the position to look for sensitive surfaces
  ///
  /// @return a vector of sensitive surfaces
  virtual std::vector<const Acts::Surface*> queryPosition(
      const Acts::GeometryContext& gctx,
      const Acts::Vector3& position) const = 0;

  /// Get all sensitive surfaces
  ///
  /// @param gctx the geometry context
  /// @param position the position to look for sensitive surfaces
  ///
  /// @return a vector of sensitive surfaces
  virtual std::vector<const Acts::Surface*> queryAll() const = 0;

  virtual ~SensitiveCandidatesBase() = default;
};

/// Implementation of the SensitiveCandidates for Gen1 geometry
struct SensitiveCandidates : public SensitiveCandidatesBase {
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;

  std::vector<const Acts::Surface*> queryPosition(
      const Acts::GeometryContext& gctx,
      const Acts::Vector3& position) const override;

  std::vector<const Acts::Surface*> queryAll() const override;
};

///   @brief The SensitiveSurfaceMapper connects the Geant 4 geometry with the Acts::TrackingGeometry.
///          The Geant 4 geometry tree is traversed and logical volumes where
///  either its name or the name of the volume's material satisfies the
///  user-configured tokens, the name of the associated physical volume is
///  prepended with the `mappingPrefix` to mark the volume as sensitive for the
///  `SensitiveSteppingAction`. In the next step, the global position of the
///  volume or of the boundary vertices is used to query the tracking geometry
///  for candidate Acts::Surfaces. In case of multiple candidates, the closest
///  one is chosen. The surface's geometry identifier is later used to associate
///  the recorded G4 hits to the surfaces.
///
///  The constructed pair of G4PhysicalVolume* and Surface* is then
///  inserted into the connection map.
class SensitiveSurfaceMapper {
 public:
  /// This prefix is used to indicate a sensitive volume that is matched
  constexpr static std::string_view mappingPrefix = "ActsSensitive#";
  /// @brief Abrivation of the association between G4 volumes and surfaces
  using VolumeToSurfAssocMap_t =
      std::multimap<const G4VPhysicalVolume*, const Acts::Surface*>;

  /// Configuration struct for the surface mapper
  struct Config {
    /// For which G4 material names we try to find a mapping
    std::vector<std::string> materialMappings;

    /// For which G4 volume names we try to find a mapping
    std::vector<std::string> volumeMappings;

    /// The sensitive surfaces that are being mapped to
    std::shared_ptr<SensitiveCandidatesBase> candidateSurfaces;
  };

  /// State object that coutns the assignments and makes
  /// a replica save copy association map
  struct State {
    /// The map of G4 physical volumes to the mapped surfaces (can be many as
    /// there can be replicas)
    VolumeToSurfAssocMap_t g4VolumeToSurfaces{};

    /// Record of the missing volumes
    std::vector<std::pair<const G4VPhysicalVolume*, Acts::Transform3>>
        missingVolumes;
  };

  /// Constructor with:
  ///
  /// @param cfg the configuration struct
  /// @param logger the logging instance
  explicit SensitiveSurfaceMapper(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "SensitiveSurfaceMapper", Acts::Logging::INFO));
  ~SensitiveSurfaceMapper() = default;

  /// Recursive mapping function that walks through the Geant4
  /// hierarchy and applies name remapping to the Physical volumes
  /// of the Geant4 geometry.
  ///
  /// @param state is the state object (for caching)
  /// @param gctx the geometry context
  /// @param g4PhysicalVolume the current physical volume in process
  /// @param motherPosition the absolute position of the mother
  void remapSensitiveNames(State& state, const Acts::GeometryContext& gctx,
                           G4VPhysicalVolume* g4PhysicalVolume,
                           const Acts::Transform3& motherTransform) const;

  /// Function that checks the success of the mapping, and exposes
  /// some additional information for debugging
  ///
  /// @param state state object after a call to remapSensitiveNames
  /// @param gctx the geometry context
  /// @param writeMissingG4VolsAsObj write the Geant4 volumes that are
  /// not mapped to 'missing_g4_volumes.obj' in the working directory
  /// @param writeMissingSurfacesAsObj write the sensitive surfaces that
  /// where not mapped to 'missing_acts_surfaces.obj' in the working directory
  /// @return Returns true only if all sensitive surfaces where mapped
  bool checkMapping(const State& state, const Acts::GeometryContext& gctx,
                    bool writeMissingG4VolsAsObj = false,
                    bool writeMissingSurfacesAsObj = false) const;

 protected:
  /// Configuration object
  Config m_cfg;

 private:
  /// Private access method to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The looging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples::Geant4
