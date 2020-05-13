// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include "Acts/Geometry/HierarchicalGeometryContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <nlohmann/json.hpp>

#include <map>

namespace Acts {

/// Convert a tracking geometry to json to initialise a Hierarchical Container
template <typename object_t>
struct JsonGeometricalHierarchyConverter {
 public:
  /// Configuration of the Reader/Writer
  class Config {
   public:
    /// The object key
    std::string datakey = "Object";
    /// The default logger
    std::shared_ptr<const Logger> logger;
    /// The name of the writer
    std::string name = "";

    /// Steering to handle sensitive data
    bool processSensitives = true;
    /// Steering to handle approach data
    bool processApproaches = true;
    /// Steering to handle representing data
    bool processRepresenting = true;
    /// Steering to handle boundary data
    bool processBoundaries = true;
    /// Steering to handle volume data
    bool processVolumes = true;

    /// Constructor
    ///
    /// @param lname Name of the writer tool
    /// @param lvl The output logging level
    Config(const std::string& lname = "JsonGeometricalHierarchyConverter",
           Logging::Level lvl = Logging::INFO)
        : logger(getDefaultLogger(lname, lvl)), name(lname) {}
  };

  using Representation = std::map<GeometryID::Value, object_t>;

  /// @brief Layer representation for Json writing
  struct LayerRep {
    // the layer id
    GeometryID layerID;
    Representation sensitives;
    Representation approaches;
    object_t representing;
  };

  /// @brief Volume representation for Json writing
  struct VolumeRep {
    // The geometry id
    GeometryID volumeID;
    /// The namne
    std::string volumeName;

    std::map<GeometryID::Value, LayerRep> layers;
    Representation boundaries;
    object_t volume;
  };

  /// @brief Detector representation for Json writing
  struct DetectorRep {
    std::map<GeometryID::Value, VolumeRep> volumes;
  };

  /// Constructor
  ///
  /// @param cfg configuration struct for the reader
  JsonGeometricalHierarchyConverter(const Config& cfg);

  /// Destructor
  ~JsonGeometricalHierarchyConverter() = default;

  /// Write json map from Tracking Geometry
  ///
  /// @param tGeometry is the tracking geometry
  nlohmann::json trackingGeometryToJson(
      const TrackingGeometry& tGeometry,
      std::function<nlohmann::json(const object_t&)> toJson,
      std::function<object_t(const GeometryID&)> initialise) const;

 private:
  /// Convert to internal representation method, recursive call
  ///
  /// @param detRep is the representation of the detector
  /// @param tVolume is the tracking volume
  /// @param initialise is a function that return an initialised Hierarchical
  /// Object
  void convertToRep(
      DetectorRep& detRep, const TrackingVolume& tVolume,
      std::function<object_t(const GeometryID&)> initialise) const;

  /// Convert to internal representation method
  ///
  /// @param tLayer is a layer
  /// @param initialise is a function that return an initialised Hierarchical
  /// Object
  LayerRep convertToRep(
      const Layer& tLayer,
      std::function<object_t(const GeometryID&)> initialise) const;

  /// Create Json from a detector represenation
  ///
  /// @param detRep is the representation of the detector
  /// @param toJson Function that convert Hierarchical Object to json
  nlohmann::json detectorRepToJson(
      const DetectorRep& detRep,
      std::function<nlohmann::json(const object_t&)> toJson) const;

  /// The detector tag
  std::string m_detkey = "detector";
  /// The volume identification string
  std::string m_volkey = "volumes";
  /// The boundary surface string
  std::string m_boukey = "boundaries";
  /// The layer identification string
  std::string m_laykey = "layers";
  /// The approach identification string
  std::string m_appkey = "approach";
  /// The sensitive identification string
  std::string m_senkey = "sensitive";
  /// The representing idntification string
  std::string m_repkey = "representing";
  /// The name identification
  std::string m_namekey = "Name";

  /// The config class
  Config m_cfg;

  /// Private access to the logging instance
  const Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace Acts
#include "Acts/Plugins/Json/JsonGeometricalHierarchyConverter.ipp"
