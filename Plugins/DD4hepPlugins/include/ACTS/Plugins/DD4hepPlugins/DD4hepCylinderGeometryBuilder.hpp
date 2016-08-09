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
    /// The logging instance
    std::shared_ptr<Logger> logger;
    /// Building the contained sub detectors
    std::shared_ptr<ITrackingVolumeBuilder> volumeBuilder;
    /// Helper tool needed for volume building
    std::shared_ptr<ITrackingVolumeHelper> volumeHelper;
    /// World detector element of the DD4hep geometry
    DD4hep::Geometry::DetElement detWorld;

    Config()
      : logger(getDefaultLogger("DD4hepCylinderGeometryBuilder", Logging::INFO))
      , volumeBuilder(nullptr)
      , volumeHelper(nullptr)
      , detWorld()
    {
    }
  };

  /// Constructor
  /// @param dbgConfig The conficuration struct for this builder
  DD4hepCylinderGeometryBuilder(const Config dgbConfig);

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

  /// TrackingGeometry Interface method - translates the DD4hep geometry into
  /// the tracking geometry
  /// @return The world TrackingGeometry object
  std::unique_ptr<TrackingGeometry>
  trackingGeometry() const override;

  /// helper to convert the TGeo transformation matrix into a ACTS transformation matrix
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
                   LayerTriple&                  layerTriple,
                   VolumeTriple&                 volumeTriple) const;

  /// Creates the cylindrical shaped layers
  /// @param motherDetElement the DD4hep detector element containing the layers
  /// @param layers a vector of layers which should get filled
  /// @param motherTransform the global transformation matrix of the
  /// motherDetElement
  void
  createCylinderLayers(DD4hep::Geometry::DetElement&      motherDetElement,
                       Acts::LayerVector&                 layers,
                       const TGeoMatrix* motherTransform) const;
  /// Creates the disc shaped layers
  /// @param motherDetElement The DD4hep detector element containing the layers
  /// @param centralLayers A vector of layers which should get filled
  /// @param motherTransform The global transformation matrix of the
  /// motherDetElement
  void
  createDiscLayers(DD4hep::Geometry::DetElement&      motherDetElement,
                   Acts::LayerVector&                 layers,
                   const TGeoMatrix* motherTransform) const;

  /// Creates a binned array of Acts::Surfaces out of vector of DD4hep detector
  /// modules
  /// @note the BinnedArray is two dimensional and one binning value will alway
  /// be phi
  /// @param modules Vector of DD4hep detector modules which should get
  /// translated into surfaces
  /// @param lValue The second binning value which can be either r or z
  /// @param motherTransform The global transformation matrix of the mother
  /// detector element containing the modules/surfaces(the layer)
  std::unique_ptr<Acts::SurfaceArray>
  createSurfaceArray(std::vector<DD4hep::Geometry::DetElement>& modules,
                     Acts::BinningValue                         lValue,
                     const TGeoMatrix*   motherTransform,
                     const std::string&  axes = "xyz") const;
  /// Creates a two dimensional surface array binned in phi and a longitudinal
  /// direction which
  /// can either be z or r
  /// @param surfaces The vector of pointers to surfaces which should be binned
  /// @param lValue The second binning value which can be either r or z
  std::unique_ptr<Acts::SurfaceArray>
  binnedSurfaceArray2DPhiL(const std::vector<const Acts::Surface*> surfaces,
                           Acts::BinningValue lValue) const;
  /// Helper method needed to get the bin values for a binned array out of
  /// overlapping modules in r
  std::vector<float>
  createBinValues(std::vector<std::pair<float, float>> old) const;
  // Helper method to sort pairs of floats
  static bool
  sortFloatPairs(std::pair<float, float> ap, std::pair<float, float> bp);

  /// Private access method to the logger
  const Logger&
  logger() const
  {
    return *m_cfg.logger;
  }

  /// configuration object
  Config m_cfg;
};

inline DD4hepCylinderGeometryBuilder::Config
DD4hepCylinderGeometryBuilder::getConfiguration() const
{
  return m_cfg;
}
}  // end of namespace

#endif  // ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
