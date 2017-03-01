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
#define ACTS_MATERIALPLUGINS_MATERIALMAPPIN_H

#include <map>
#include <utility>
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialTrackRecord.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class Layer;
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
/// - material steps recorded from the detailed geometry
/// - a prepared Acts::TrackingGeometry with Acts::SurfaceMaterialProxy ob
///   surfaces when the mapping should be done  
/// 
/// In a first step all surfaces in the TrackingGeometry with a material proxy
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
    // ignore events with eta bigger than the cutoff value @todo add later
    //         double etaCutoff;
    // needed for debugging: -1 negative | 0 all | 1 positive @todo add later
    //          int etaSide;
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
  /// @param cnf the internal configuration object
  /// @param logger the logging instance
  MaterialMapping(const Config&           cnf,
                  std::unique_ptr<Logger> logger
                  = getDefaultLogger("MaterialMapping", Logging::INFO));

  /// @brief destructor
  ~MaterialMapping();

  /// maps the material for the given direction(eta,phi) onto the layers of the
  /// given tracking geometry
  /// 
  /// @param matTrackRec the Acts::MaterialTrackRecord to be mapped
  ///
  /// @return if the mapping was successful
  bool
  mapMaterialTrackRecord(Cache& mappingCache, 
                         const MaterialTrackRecord& matTrackRec) const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

private:
  /// internally used method to collect hits on their corresponding layers
  /// received when extrapolating through the tracking geometry
  bool
  collectLayersAndHits(
      const MaterialTrackRecord& matTrackRec,
      std::vector<std::pair<GeometryID, Acts::Vector3D>>&
          surfacesAndHits);
      
  /// internally used method to associate the material to the right layer in the
  /// tracking geometry
  void
  associateLayerMaterial(
      const MaterialTrackRecord& matTrackRec,
      std::vector<std::pair<GeometryID, Vector3D>>&
          surfacesAndHits);
      
  /// internally used method to associate a hit to a given layer by recording it
  /// in the layer records map
  void
  associateHit(const Layer*                     layer,
               const Vector3D&                  position,
               const std::vector<MaterialStep>& layerMaterialSteps);
  /// configuration object

  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the configuration object
  Config m_cnf;
  
  /// the logging instance
  std::unique_ptr<Logger> m_logger;
  
};
}

inline const std::map<const Acts::Layer*, Acts::SurfaceMaterialRecord>
Acts::MaterialMapping::layerRecords() const
{
  return m_layerRecords;
}

#endif  // ACTS_MATERIALPLUGINS_MATERIALMAPPIN_Hr
