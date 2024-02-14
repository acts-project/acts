// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/interface/IMaterialMapper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/VolumeCollector.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {

class IVolumeMaterial;
class TrackingGeometry;
class VolumeBounds;

namespace Experimental {
class Detector;
class DetectorVolume;
}  // namespace Experimental

/// @brief VolumeMaterialMapper
///
/// This is the main feature tool to map material information
/// from a 3D geometry onto the TrackingGeometry with its surface
/// material description.
///
/// The process runs as such:
///
///  1) Input geometry is parsed and for each Volume with
///     ProtoVolumeMaterial a local store is initialized
///     the identification is done hereby through the
///     'Volume'::GeometryIdentifier
///
/// @note that this mapper is designed to work with both the `TrackingGeometry` and
/// the `Detector` schema as long as they co-exist.
///
///  2) A number of N material tracks is read in, each track has :
///       origin, direction, material steps (< position, step length, x0, l0, a,
///       z, rho >, thichness)
///
///       for each track:
///          volume along the origin/direction path are collected.
///          the step are then associated to volume inside which they are.
///          Additional step are created along the track direction.
///
///  3) Each 'hit' bin per event is counted and averaged at the end of the run

class VolumeMaterialMapper final : public IMaterialMapper {
 public:
  using StraightLineTGPropagator = Propagator<StraightLineStepper, Navigator>;
  using StraightLineDetPropagator =
      Propagator<StraightLineStepper, Experimental::DetectorNavigator>;

  /// @struct Config
  ///
  /// Nested Configuration struct for the material mapper
  struct Config {
    /// Size of the step for the step extrapolation
    ActsScalar mappingStep = 1.;
    /// A potential veto
    MaterialInteractionVeto veto = NoVeto{};
  };

  /// @struct State
  ///
  /// Nested State struct which is used for the mapping prococess
  struct State final : public MaterialMappingState {
    /// Constructor of the State with contexts
    State(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : geoContext(gctx), magFieldContext(mctx) {}

    /// The recorded material per geometry ID
    std::map<const GeometryIdentifier, Acts::AccumulatedVolumeMaterial>
        homogeneousGrid;

    /// The recorded 2D transform associated the grid for each geometry ID
    std::map<const GeometryIdentifier,
             std::function<Acts::Vector2(Acts::Vector3)>>
        transform2D;

    /// The 2D material grid for each geometry ID
    std::map<const GeometryIdentifier, Grid2D> grid2D;

    /// Recorded 3D transform associated the material grid for each geometry ID
    std::map<const GeometryIdentifier,
             std::function<Acts::Vector3(Acts::Vector3)>>
        transform3D;

    /// The 3D material grid for each geometry ID
    std::map<const GeometryIdentifier, Grid3D> grid3D;

    /// The binning for each geometry ID
    std::map<const GeometryIdentifier, BinUtility> materialBin;

    /// Reference to the geometry context for the mapping
    std::reference_wrapper<const GeometryContext> geoContext;

    /// Reference to the magnetic field context
    std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  };

  /// Delete the Default constructor
  VolumeMaterialMapper() = delete;

  /// Constructor with config object
  ///
  /// @param cfg Configuration struct
  /// @param propagator The straight line propagator with the TrackingGeometry navigation
  /// @param slogger The logger
  VolumeMaterialMapper(const Config& cfg, StraightLineTGPropagator& propagator,
                       std::unique_ptr<const Logger> slogger = getDefaultLogger(
                           "VolumeMaterialMapper", Logging::INFO));

  /// Constructor with config object
  ///
  /// @param cfg Configuration struct
  /// @param propagator The straight line propagator with the Detector navigation
  /// @param slogger The logger
  VolumeMaterialMapper(const Config& cfg, StraightLineDetPropagator& propagator,
                       std::unique_ptr<const Logger> slogger = getDefaultLogger(
                           "VolumeMaterialMapper", Logging::INFO));

  /// @brief helper method that creates the cache for the mapping
  ///
  /// @param[in] gctx The geometry context to use
  /// @param[in] mctx The magnetic field context to use
  /// @param[in] tGeometry The geometry which should be mapped
  ///
  /// This method takes a TrackingGeometry,
  /// finds all volume material proxes and fills them into a State
  /// object to be used in the mapping process
  std::unique_ptr<MaterialMappingState> createState(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const TrackingGeometry& tGeometry) const final;

  /// @brief helper method that creates the cache for the mapping
  ///
  /// @param[in] gctx The geometry context to use
  /// @param[in] mctx The magnetic field context to use
  /// @param[in] detector The detector which should be mapped
  ///
  /// This method takes a Detector object,
  /// finds all volume material proxes and fills them into a State
  /// object to be used in the mapping process
  std::unique_ptr<MaterialMappingState> createState(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const Experimental::Detector& detector) const final;

  /// @brief Method to finalize the maps
  ///
  /// It calls the final run averaging and then transforms
  /// the Homogeneous material into HomogeneousVolumeMaterial and
  /// the 2D and 3D grid into a InterpolatedMaterialMap
  ///
  /// @param mState
  MaterialMappingResult finalizeMaps(MaterialMappingState& mState) const final;

  /// Process/map a single track
  ///
  /// @param mState The current state map
  /// @param mTrack The material track to be mapped
  ///
  /// @note the RecordedMaterialSlab of the track are assumed
  /// to be ordered from the starting position along the starting direction
  ///
  /// @note it will @return the mapped and unmapped part of the material track
  std::array<RecordedMaterialTrack, 2u> mapMaterialTrack(
      MaterialMappingState& mState, const RecordedMaterialTrack& mTrack) const final;

 private:
  /// selector for finding
  struct MaterialVolumeSelector {
    template <typename volume_type>
    bool operator()(const volume_type& vf) const {
      return (vf.volumeMaterial() != nullptr);
    }
  };

  /// @brief finds all TrackingVolume objects with ProtoVolumeMaterial
  ///
  /// @param mState The state to be filled
  /// @param tVolume is current TrackingVolume
  void resolveMaterialVolume(State& mState,
                             const TrackingVolume& tVolume) const;

  /// @brief finds all DetectorVolume objects with ProtoVolumeMaterial
  ///
  /// @param mState The state to be filled
  /// @param dVolume is current TrackingVolume
  void resolveMaterialVolume(State& mState,
                             const Experimental::DetectorVolume& dVolume) const;

  /// @brief check and insert
  ///
  /// @param mState is the map to be filled
  /// @param volumeMaterial is the material found for this volume
  /// @param volumeBounds the bounds of the volume
  /// @param transform the transform of the volume
  /// @param geoID is the volume geometry identifier which is used to store
  void checkAndInsert(State& mState, const IVolumeMaterial& volumeMaterial,
                      const VolumeBounds& volumeBounds,
                      const Transform3& transform,
                      const GeometryIdentifier& geoID) const;

  /// Create extra material point for the mapping and add them to the grid
  ///
  /// @param mState The state to be filled
  /// @param currentBinning a pair containing the current geometry ID and the current binning
  /// @param properties material properties of the original hit
  /// @param position position of the original hit
  /// @param direction direction of the track
  void createExtraHits(
      State& mState,
      std::pair<const GeometryIdentifier, BinUtility>& currentBinning,
      Acts::MaterialSlab properties, const Vector3& position,
      Vector3 direction) const;

  /// Standard logger method
  const Logger& logger() const { return *m_logger; }

  /// The configuration object
  Config m_cfg;

  /// The straight line propagator
  std::shared_ptr<const StraightLineTGPropagator> m_tgPropagator = nullptr;
  std::shared_ptr<const StraightLineDetPropagator> m_detPropagator = nullptr;

  /// The logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
