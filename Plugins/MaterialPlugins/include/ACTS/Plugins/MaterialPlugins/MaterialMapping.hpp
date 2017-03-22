// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialMapping.h, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIALPLUGINS_MATERIALMAPPIN_H
#define ACTS_MATERIALPLUGINS_MATERIALMAPPIN_H 1

#include <map>
#include <utility>
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialTrackRecord.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class Surface;
class TrackingGeometry;
class MaterialProperties;
class SurfaceMaterialRecord;

/// @class MaterialMapping
///
/// @brief Class for material mapping
///
/// This class should be used to map material from the full and detailed
/// detector geometry onto the simplified ACTS geometry. It offers options to
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
/// One MaterialTrackRecord (containing all the MaterialSteps along a Track) is
/// mapped by using the function Acts::MaterialMapping::mapMaterial(). The
/// mapping process starts by extrapolating from the same starting
/// point and direction as the MaterialTrackRecord through the ACTS geometry.
/// The extrapolation engine then finds  the closest surface marked to carry  
/// material (by carrying a SurfaceMaterialProxy).
/// The material steps are then assigned to the corresponding surfaces 
/// (and the according bin) at the assigned position.
///
/// Along one track in one bin of a layer the material is averaged:
/// \image html MaterialAveraging.jpeg
///
/// When the material mapping is done many MaterialTrackRecords will be mapped.
/// Everytime the same bin is hit, the material parameters are summed up.
/// This information is cached in the corresponding
/// SurfaceMaterialRecord object.
/// 
/// In a finalization step, the SurfaceMaterialRecord bins are averaged
/// by the number of hits per bin and the final BinnedSufaceMaterial
/// are created.

class MaterialMapping
{
public:
  /// @struct Config
  ///
  /// Configuration for the MaterialMapping
  struct Config
  {
    /// ignore events with eta bigger than the cutoff value
    double etaCutoff;
    /// extrapolation engine
    std::shared_ptr<IExtrapolationEngine> extrapolationEngine = nullptr;
  };

  /// @struct Cache 
  /// 
  /// This is the cache object used for calling the mapping method
  struct Cache {
    /// object which connects the layer with its SurfaceMaterialRecord
    std::map<GeometryID, SurfaceMaterialRecord> surfaceMaterialRecords;
    
  };

  /// @brief default constructor
  ///
  /// @param cfg the internal configuration object
  /// @param logger the logging instance
  MaterialMapping(const Config&           cfg,
                  std::unique_ptr<Logger> logger
                  = getDefaultLogger("MaterialMapping", Logging::INFO));

  /// @brief destructor
  ~MaterialMapping();

  /// @brief helper method that creates the cache for the mapping
  ///
  /// This method takes a TrackingGeometry, finds all surfaces with
  /// material proxis 
  std::map<GeometryID, SurfaceMaterialRecord >
  materialMappingCache(const TrackingGeometry& tGeometry) const;
  
  /// maps the material for the given direction(eta,phi) onto the layers of the
  /// given tracking geometry
  /// 
  /// @param matTrackRec the MaterialTrackRecord to be mapped
  ///
  /// @return if the mapping was successful
  bool
  mapMaterialTrackRecord(Cache& mappingCache, 
                         const MaterialTrackRecord& matTrackRec) const;

  /// finds the TrackingGeometry steps associated to the material steps
  ///
  /// @param materialSteps are the full geometry steps
  /// @param assignedSteps are the Tracking geometry associated points
  void
  assignSteps(const std::vector<MaterialStep>& materialSteps,
              std::vector< AssignedSteps >& assignedSteps) const;
  
  /// set logging instance
  ///
  /// @param logger is the unique logger instance
  void
  setLogger(std::unique_ptr<Logger> logger);

private:
  
  
  /// finds all surfaces with SurfaceMaterialProxy of a volume
  ///
  /// @param sMap is the map to be filled
  /// @param tVolume is current TrackingVolume
  void
  collectMaterialSurfaces(std::map<GeometryID, SurfaceMaterialRecord>& sMap,
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
  std::unique_ptr<Logger> m_logger;
  
};
}

#endif  // ACTS_MATERIALPLUGINS_MATERIALMAPPING_H
