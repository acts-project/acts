// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialMapper.h, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include <map>
#include <utility>
#include "Acts/Extrapolation/IExtrapolationEngine.hpp"
#include "Acts/Plugins/MaterialMapping/AssignedMaterialSteps.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialStep.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialTrack.hpp"
#include "Acts/Plugins/MaterialMapping/SurfaceMaterialRecord.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class Surface;
class TrackingGeometry;
class MaterialProperties;

/// @class MaterialMapper
///
/// @brief Class for material mapping
///
/// This class should be used to map material from the full and detailed
/// detector geometry onto the simplified Acts geometry. It offers options to
/// map, average and finalize the material.
///
/// Preconditions are:
/// - material steps recorded from the detailed geometry (e.g. from Geant4)
/// - a prepared Acts::TrackingGeometry with Acts::SurfaceMaterialProxy ob
///   surfaces when the mapping should be done
///
/// All surfaces of the TrackingGeometry with a material proxy
/// are identified and SurfaceMaterialRecords are created.
///
/// @todo update the following description
///
/// One MaterialTrack (containing all the MaterialSteps along a Track) is
/// mapped by using the function Acts::MaterialMapper::mapMaterial(). The
/// mapping process starts by extrapolating from the same starting
/// point and direction as the MaterialTrack through the Acts geometry.
/// The extrapolation engine then finds  the closest surface marked to carry
/// material (by carrying a SurfaceMaterialProxy).
/// The material steps are then assigned to the corresponding surfaces
/// (and the according bin) at the assigned position.
///
/// Along one track in one bin of a layer the material is averaged:
/// \image html MaterialAveraging.jpeg
///
/// When the material mapping is done many MaterialTracks will be mapped.
/// Everytime the same bin is hit, the material parameters are summed up.
/// This information is cached in the corresponding
/// SurfaceMaterialRecord object.
///
/// In a finalization step, the SurfaceMaterialRecord bins are averaged
/// by the number of hits per bin and the final BinnedSufaceMaterial
/// are created.

class MaterialMapper
{
public:
  /// @struct Config
  ///
  /// Configuration for the MaterialMapper
  struct Config
  {
    /// ignore events with eta bigger than the cutoff value
    double etaCutoff;
    /// extrapolation engine
    std::shared_ptr<const IExtrapolationEngine> extrapolationEngine = nullptr;
  };

  /// @struct Cache
  ///
  /// This is the cache object used for calling the mapping method
  struct Cache
  {

    /// object which connects the layer with its SurfaceMaterialRecord
    std::map<GeometryID, SurfaceMaterialRecord> surfaceMaterialRecords;

    // counter in case one wants to combine output from several jobs
    size_t materialTrackCounter = 0;

    /// Constructor from a new map
    Cache(std::map<GeometryID, SurfaceMaterialRecord> smr)
      : surfaceMaterialRecords(std::move(smr)), materialTrackCounter(0)
    {
    }
  };

  /// @brief default constructor
  ///
  /// @param cfg the internal configuration object
  /// @param log the logging instance
  MaterialMapper(const Config&                 cfg,
                 std::unique_ptr<const Logger> log
                 = getDefaultLogger("MaterialMapper", Logging::INFO));

  /// @brief destructor
  ~MaterialMapper();

  /// @brief helper method that creates the cache for the mapping
  ///
  /// This method takes a TrackingGeometry,
  /// finds all surfaces with material proxis
  /// and returns you a Cache object ot be used
  Cache
  materialMappingCache(const TrackingGeometry& tGeometry) const;

  /// maps the material for the given direction(eta,phi) onto the layers of the
  /// given tracking geometry
  ///
  /// @param materialTrack the MaterialTrack to be mapped
  ///
  /// @return is the mapped material track, i.e. it is collapsed
  ///      onto the available
  MaterialTrack
  mapMaterialTrack(Cache&               mappingCache,
                   const MaterialTrack& materialTrack) const;

  /// finds the TrackingGeometry steps associated to the material steps
  ///
  /// @param materialSteps are the full geometry steps
  /// @param assignedSteps are the Tracking geometry associated points
  ///
  /// @note this method is currently public for Testing
  void
  assignSteps(const std::vector<MaterialStep>&    materialSteps,
              std::vector<AssignedMaterialSteps>& assignedSteps) const;

  /// creates the final surface material records
  ///
  /// @param mappingCache
  std::map<GeometryID, SurfaceMaterial*>
  createSurfaceMaterial(Cache& mappingCache) const;

  /// set logging instance
  ///
  /// @param newLogger is the unique logger instance
  void
  setLogger(std::unique_ptr<const Logger> newLogger);

private:
  /// finds all surfaces with SurfaceMaterialProxy of a volume
  ///
  /// @param sMap is the map to be filled
  /// @param tVolume is current TrackingVolume
  void
  resolveMaterialSurfaces(std::map<GeometryID, SurfaceMaterialRecord>& sMap,
                          const TrackingVolume& tVolume) const;

  /// check and insert
  ///
  /// @param sMap is the map to be filled
  /// @param surface is the surface to be checked for a Proxy
  void
  checkAndInsert(std::map<GeometryID, SurfaceMaterialRecord>& sMap,
                 const Surface& surface) const;

  /// standard logger method
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the configuration object
  Config m_cfg;

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;
};
}