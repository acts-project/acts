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
class LayerMaterialRecord;

/// @class MaterialMapping
///
/// @brief Class for material mapping
///
/// This class should be used to map material from the full and detailed
/// detector geometry onto the simplified ACTS geometry. It offers options to
/// map, average and finalize the material.
///

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
  /// @param matTrackRec the Acts::MaterialTrackRecord to be mapped
  void
  mapMaterial(const MaterialTrackRecord& matTrackRec);
  /// averages the material of the layer records collected so far (for each bin)
  void
  averageLayerMaterial();
  /// after all step collections have been mapped this method needs to be called
  /// it sets the created material to the layers
  void
  finalizeLayerMaterial();
  /// @return hands back all layers carrying material with their corresponding
  /// Acts::LayerMaterialRecord
  const std::map<const Layer*, LayerMaterialRecord>
  layerRecords() const;
  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

private:
  /// internally used method to collect hits on their corresponding layers
  /// received when extrapolating through the tracking geometry
  bool
  collectLayersAndHits(
      const MaterialTrackRecord& matTrackRec,
      std::map<const Acts::Layer*, Acts::Vector3D>& layersAndHits);
  /// internally used method to associate the material to the right layer in the
  /// tracking geometry
  void
  associateLayerMaterial(
      const MaterialTrackRecord& matTrackRec,
      std::map<const Acts::Layer*, Acts::Vector3D>& layersAndHits);
  /// internally used method to associate a hit to a given layer by recording it
  /// in the layer records map
  void
  associateHit(const Layer*                           layer,
               const Acts::Vector3D&                  position,
               const std::vector<Acts::MaterialStep>& layerMaterialSteps);
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
  /// object which connects the layer with its LayerMaterialRecord
  std::map<const Layer*, LayerMaterialRecord> m_layerRecords;
};
}

inline const std::map<const Acts::Layer*, Acts::LayerMaterialRecord>
Acts::MaterialMapping::layerRecords() const
{
  return m_layerRecords;
}

#endif  // ACTS_MATERIALPLUGINS_MATERIALMAPPIN_Hr
