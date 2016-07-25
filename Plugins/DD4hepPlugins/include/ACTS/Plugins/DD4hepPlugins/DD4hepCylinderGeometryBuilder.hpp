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

// Geometry module
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Tools/ITrackingGeometryBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Logger.hpp"
// DD4hep
#include "DD4hep/Detector.h"

namespace Acts {
class TrackingGeometry;
}
namespace Acts {

/** @ class DD4hepCylinderGeometryBuilder

This class receives the DD4hep geometry from the given implementation of the
IDD4hepGeomtrySvc, walks through the subvolumes and initiates their translation
into the tracking geometry.
    It returns the world tracking geometry element.
 @TODO find replacement for Gaudi exeption and message stream

 */

class DD4hepCylinderGeometryBuilder
    : virtual public Acts::ITrackingGeometryBuilder
{
public:
  /** @struct Config
   Configuration for the DD4hepCylinderGeometryBuilder */
  struct Config
  {
      /// the logging instance
      std::shared_ptr<Logger>                            logger;
    std::shared_ptr<ITrackingVolumeBuilder>
        volumeBuilder;  //!< building the contained sub detectors
    std::shared_ptr<ITrackingVolumeHelper>
        volumeHelper;  //!< helper tool needed for volume building
    DD4hep::Geometry::DetElement
        detWorld;  //!< world detector element of the DD4hep geometry

    Config() : logger(getDefaultLogger("DD4hepCylinderGeometryBuilder", Logging::INFO)), volumeBuilder(nullptr), volumeHelper(nullptr), detWorld() {}
  };

  /** Constructor */
  DD4hepCylinderGeometryBuilder(const Config dgbConfig);

  /** Destructor */
  virtual ~DD4hepCylinderGeometryBuilder();

  /** setting the builders and helpers with the configuration object*/
  void
  setConfiguration(const Config dgbConfig);

  /** get the configuration object */
  Config
  getConfiguration() const;

  /** TrackingGeometry Interface method - translates the DD4hep geometry into
   * the tracking geometry*/
  std::unique_ptr<TrackingGeometry>
  trackingGeometry() const override;
    
    /**helper method to extract the transformation matrix from a DD4hep
     * DetElement*/
    std::shared_ptr<Acts::Transform3D>
    extractTransform(DD4hep::Geometry::DetElement& detElement) const;
    /**helper method to extract the volume boundaries of a cylindrical volume*/
    std::shared_ptr<const Acts::VolumeBounds>
    extractVolumeBounds(DD4hep::Geometry::DetElement& detElement) const;

private:
    
    /** Creates a triple of volumes a possible barrel-endcap configuration and of
     * all the three possible Layer types of the given volume detector element*/
    /** constructs all subvolumes contained by this volume (motherDetELement) with
     * its layers and modules, if present */
    void
    createSubVolumes(DD4hep::Geometry::DetElement& motherDetElement,
                     LayerTriple&                  layerTriple,
                     VolumeTriple&                 volumeTriple) const;
    
    /**creates the cylindrical shaped layers*/
    void
    createCylinderLayers(DD4hep::Geometry::DetElement&      motherDetElement,
                         Acts::LayerVector&                 centralLayers,
                         std::shared_ptr<Acts::Transform3D> motherTransform
                         = nullptr) const;
    /**creates disc shaped layers*/
    void
    createDiscLayers(DD4hep::Geometry::DetElement&      motherDetElement,
                     Acts::LayerVector&                 layers,
                     std::shared_ptr<Acts::Transform3D> motherTransform
                     = nullptr) const;
    
    /**creates a binned array of Acts::Surfaces out of vector of DD4hep detector
     * modules*/
    std::unique_ptr<Acts::SurfaceArray>
    createSurfaceArray(std::vector<DD4hep::Geometry::DetElement>& modules,
                       Acts::BinningValue                         lValue,
                       std::shared_ptr<const Acts::Transform3D>   motherTransform
                       = nullptr) const;
    /**creating a surface array binned in phi and a longitudinal direction which
     * can either be z or r*/
    std::unique_ptr<Acts::SurfaceArray>
    binnedSurfaceArray2DPhiL(const std::vector<const Acts::Surface*> surfaces,
                             Acts::BinningValue                      lValue) const;
    /**helper method to get the bin values for a binned array out of overlapping
     * modules*/
    std::vector<float>
    createBinValues(std::vector<std::pair<float, float>> old) const;
    /**helper method to sort pairs of doubles*/
    static bool
    sortFloatPairs(std::pair<float, float> ap, std::pair<float, float> bp);
    
    /// Private access method to the logger
    const Logger&
    logger() const
    {
        return *m_cfg.logger;
    }
    
    
  /** configuration object */
  Config m_cfg;
};

inline DD4hepCylinderGeometryBuilder::Config
DD4hepCylinderGeometryBuilder::getConfiguration() const
{
  return m_cfg;
}
}  // end of namespace

#endif  //#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPCYLINDERGEOMETRYBUILDER_H
