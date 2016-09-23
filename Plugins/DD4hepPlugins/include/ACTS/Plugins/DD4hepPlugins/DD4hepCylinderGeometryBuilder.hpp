// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepCylinderGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H 1

#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Tools/ISurfaceArrayCreator.hpp"
#include "ACTS/Tools/ITrackingGeometryBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"
#include "DD4hep/Detector.h"

namespace Acts {
class TrackingGeometry;
}
namespace Acts {

/// @class DD4hepCylinderGeometryBuilder
///
/// @brief builds the tracking geometry from DD4hep input
///
/// This class receives the DD4hep geometry over the configuration
/// as the world DD4hep::DetElement. It walks through the DD4hep detector
/// hierarchy
/// and translates it to the ACTS tracking geometry. It walks through the
/// possible subdetectors,
/// which can contain layers. The layers themselves can contain modules and/or
/// are support structure.
/// It returns the world tracking geometry element.

class DD4hepCylinderGeometryBuilder
    : virtual public Acts::ITrackingGeometryBuilder
{
public:
  /// @struct Config
  /// Configuration for the DD4hepCylinderGeometryBuilder
  struct Config
  {
    /// Building the contained sub detectors
    std::shared_ptr<ITrackingVolumeBuilder> volumeBuilder = nullptr;
    /// Helper tool needed for volume building
    std::shared_ptr<ITrackingVolumeHelper> volumeHelper = nullptr;
    /// Surface array creator needed for building binned surface arrays
    std::shared_ptr<ISurfaceArrayCreator> surfaceArrayCreator = nullptr;
    /// World detector element of the DD4hep geometry
    DD4hep::Geometry::DetElement detWorld;
  };

  /// Constructor
  /// @param dbgConfig The conficuration struct for this builder
  /// @param logger logging instance
  DD4hepCylinderGeometryBuilder(
      const Config            dgbConfig,
      std::unique_ptr<Logger> logger
      = getDefaultLogger("DD4hepCylinderGeometryBuilder", Logging::INFO));

  /// Destructor
  virtual ~DD4hepCylinderGeometryBuilder() = default;

  /// setting the builders and helpers with the configuration object
  /// @param dbgConfig The conficuration struct for this builder
  void
  setConfiguration(const Config dgbConfig);

  /// get the configuration object
  /// @return The conficuration struct for this builder
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

  /// TrackingGeometry Interface method - translates the DD4hep geometry into
  /// the tracking geometry
  /// @return The world TrackingGeometry object
  std::unique_ptr<TrackingGeometry>
  trackingGeometry() const override;

  /// helper to convert the TGeo transformation matrix into a ACTS
  /// transformation matrix
  /// @param tGeoTrans TGeo transformation matrix which should be converted
  std::shared_ptr<Acts::Transform3D>
  convertTransform(const TGeoMatrix* tGeoTrans) const;
  /// helper method to extract the ACTS volume boundaries of a cylindrical
  /// volume
  /// @note Only cylindrical volume boundary implemented
  /// @param detElement DD4hep detector element of which the volume boundaries
  /// should be returned
  std::shared_ptr<const Acts::VolumeBounds>
  extractVolumeBounds(DD4hep::Geometry::DetElement& detElement) const;

private:
  /// Creates a triple of volumes of a possible barrel-endcap configuration and
  /// of
  /// all the three possible Layer types of the given volume detector element.
  /// constructs all subvolumes contained by this volume (motherDetELement) with
  /// its layers and modules, if present
  /// @param motherDetElement The DD4hep detector element containing the sub
  /// volumes
  /// @param layerTriple The layer triple which should get filled
  /// @param volumeTriple The volume triple which whould get filled
  void
  createSubVolumes(DD4hep::Geometry::DetElement& motherDetElement,
                   LayerTriple&                  layerTriple) const;

  /// Creates the cylindrical shaped layers
  /// @param motherDetElement the DD4hep detector element containing the layers
  /// @param layers a vector of layers which should get filled
  /// @param motherTransform the global transformation matrix of the
  /// motherDetElement
  void
  createCylinderLayers(DD4hep::Geometry::DetElement& motherDetElement,
                       Acts::LayerVector&            layers,
                       const TGeoMatrix*             motherTransform) const;
  /// Creates the disc shaped layers
  /// @param motherDetElement The DD4hep detector element containing the layers
  /// @param centralLayers A vector of layers which should get filled
  /// @param motherTransform The global transformation matrix of the
  /// motherDetElement
  void
  createDiscLayers(DD4hep::Geometry::DetElement& motherDetElement,
                   Acts::LayerVector&            layers,
                   const TGeoMatrix*             motherTransform) const;

  /// Creates a std::vector of Acts::Surface out of std::vector of DD4hep
  /// DetELements
  /// @param modules Vector of DD4hep detector modules which should get
  /// translated into surfaces
  /// @param motherTransform The global transformation matrix of the mother
  /// @param axes the orientation of the modules
  std::vector<const Acts::Surface*>
  createSurfaceVector(std::vector<DD4hep::Geometry::DetElement>& modules,
                      const TGeoMatrix*  motherTransform,
                      const std::string& axes = "xyz") const;

  /// Private access method to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// configuration object
  Config m_cfg;

  /// The logging instance
  std::unique_ptr<Logger> m_logger;
};

inline DD4hepCylinderGeometryBuilder::Config
DD4hepCylinderGeometryBuilder::getConfiguration() const
{
  return m_cfg;
}
}  // end of namespace

#endif  // ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
