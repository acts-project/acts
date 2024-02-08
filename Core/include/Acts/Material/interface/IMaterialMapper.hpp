// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

class ISurfaceMaterial;
class IVolumeMaterial;
class TrackingGeometry;

namespace Experimental {
class Detector;
}

struct MaterialMappingState {
  virtual ~MaterialMappingState() = default;
};

struct MaterialMappingResult {
  std::map<GeometryIdentifier, std::unique_ptr<const ISurfaceMaterial>>
      surfaceMaterial;
  std::map<GeometryIdentifier, std::unique_ptr<IVolumeMaterial>> volumeMaterial;
};

using MaterialInteractionVeto = std::function<bool(const MaterialInteraction&)>;

struct NoVeto {
  bool operator()(const MaterialInteraction&) const { return false; }
};

/// @brief Interface definition of Material mappers, it defines how to
/// create a mapping state either from a TrackingGeometry or a Detector
/// object, which is then handled by the :
/// - mapMaterialTrack method
/// - finalizeMaps method
///
/// @note IMaterialMapper instances are designed to work in a multi-threaded
/// environment via mutex-protected state objects.
///
/// @note the individual material mapping instances, are processing a
/// recorded material track and return the unused interaction records
/// as a new material track for subsequent mappers to process.
///
/// The mapping function thus returns a pair of Recorded Material tracks:
/// - the ones used for mapping
/// - the unprocessed ones
class IMaterialMapper {
 public:
  /// Virtual destructor
  virtual ~IMaterialMapper() = default;

  /// @brief Create the state object from a TrackingGeometry
  ///
  /// @param geoContext is the geometry context
  /// @param magFieldContext is the magnetic field context
  /// @param trackingGeometry is the tracking geometry
  ///
  /// @return the state object
  virtual std::unique_ptr<MaterialMappingState> createState(
      const GeometryContext& geoContext,
      const MagneticFieldContext& magFieldContext,
      const TrackingGeometry& trackingGeometry) const = 0;

  /// @brief Create the state object from a Detector object
  ///
  /// @param geoContext is the geometry context
  /// @param magFieldContext is the magnetic field context
  /// @param detector is the tracking geometry
  ///
  /// @return the state object
  virtual std::unique_ptr<MaterialMappingState> createState(
      const GeometryContext& geoContext,
      const MagneticFieldContext& magFieldContext,
      const Experimental::Detector& detector) const = 0;

  /// @brief Mapping interface method
  ///
  /// The concrete implementation of the Material Mapper is supposed
  /// to remove the mapped material records from the provided material
  /// track and to return the remaining material track.
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note the RecordedMaterialSlab of the track are assumed
  /// to be ordered from the starting position along the starting direction
  ///
  /// @return the mapped part and unused part of the material track
  virtual std::array<RecordedMaterialTrack, 2u> mapMaterialTrack(
      MaterialMappingState& mState,
      const RecordedMaterialTrack& mTrack) const = 0;

  /// @brief Finalize the maps, which usually involves the merging of
  /// event-wise collected information
  ///
  /// @param state is the state object
  ///
  /// @return a Mapping result
  virtual MaterialMappingResult finalizeMaps(
      MaterialMappingState& state) const = 0;
};

}  // namespace Acts
