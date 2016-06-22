// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H
#define ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H 1

// Geometry module
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Logger.hpp"

#ifndef ATAS_GEOMETRYTOOLS_TAKESMALLERBIGGER
#define ATAS_GEOMETRYTOOLS_TAKESMALLERBIGGER
#define takeSmaller(current, test) current = current < test ? current : test
#define takeBigger(current, test) current  = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test)                           \
  takeSmaller(cSmallest, test);                                                \
  takeBigger(cBiggest, test)
#endif

namespace Acts {

class TrackingVolume;
class VolumeBounds;

/// LayerSetup struct to understand the layer setup
struct LayerSetup
{
  bool         present;       ///< layers are present
  BinningValue binningValue;  ///< in what way they are binned

  std::pair<double, double> rBoundaries;  //!< raidal boundaries
  std::pair<double, double> zBoundaries;  //!< zBoundaries

  /// std::vector<double>      ringBoundaries;      //!< ring boundaries if
  /// present //!< @TODO insert ring layout
  
  LayerSetup()
    : present(false)
    , binningValue(binR)
    , rBoundaries(std::pair<double, double>(10e10, -10e10))
    , zBoundaries(std::pair<double, double>(10e10, -10e10))
  {
  }

  /// Conversion operator to bool 
  operator bool() const { return present; }
};

/// @class CylinderVolumeBuilder
///
///  A simple cylindrical volume builder to be used for building a concentrical
/// cylindrical volume
///  - a) configured volume
///  - b) wrapping around a cylindrical/disk layer setup
///
///  All are optionally wrapped around a given volume which has to by a cylinder
/// volume
///  and which has to be center at z == 0
///
///  To receive the tracking volume it is possible to also hand over a triple of
/// layers, which is a C++ tuple of three pointers to layer vectors (defined in
/// the ITrackingVolumeBuilder). This functionality is needed for a possible
/// translation of an geometry existing in another format. The first entry
/// represents the layers of the negative endcap, the second the layers of the
/// barrel and the third the layers of the positive endcap. If the one of these
/// pointers is a nullptr no layers will be created for this volume
///  Another functionality needed to translate an already existing geometry is to
/// hand over a volume triple, which is a triple of shared pointers of volumes
/// (defined in the ITrackingVolumeBuilder). The first entry contains the
/// negative endcap volume, the second the barrel volume and the third one the
/// positive endcap volume. This volumes are then used to get the internal
/// boundaries of the current hierarchy.
///
class CylinderVolumeBuilder : public ITrackingVolumeBuilder
{
public:
  /// @struct Config
  /// Nested configuration struct for this CylinderVolumeBuilder
  struct Config
  {
    /// the logging instance
    std::shared_ptr<Logger>                logger;        
    /// the trackign volume helper for construction
    std::shared_ptr<ITrackingVolumeHelper> trackingVolumeHelper;  
    /// the string based indenfication
    std::string                            volumeName;  
    /// The dimensions of the manually created world
    std::vector<double>                    volumeDimension; 
    /// the world material 
    std::shared_ptr<Material>              volumeMaterial;  
    /// build the volume to the beam pipe
    bool                                    volumeToBeamPipe; 
    /// needed to build layers within the volume 
    std::shared_ptr<ILayerBuilder>          layerBuilder;
    /// the envelope covering the potential layers     
    double                                  layerEnvelopeR;
    /// the envelope covering the potential layers   
    double                                  layerEnvelopeZ; 
    /// the volume signature  
    int                                     volumeSignature;  

    Config() : logger(getDefaultLogger("CylinderVolumeBuilder", Logging::INFO))
    {
    }
  };

  /// Constructor 
  /// @param cvbConfig is the configuraiton struct to steer the builder
  CylinderVolumeBuilder(const Config& cvbConfig);

  /// Destructor 
  virtual ~CylinderVolumeBuilder();

  /// CylinderVolumeBuilder interface method  
  /// @param insideVolume is an (optional) volume to be wrapped
  /// @param outsideBounds is an (optional) outside confinement
  /// @param layerTriple is an (optional) triplet of layers 
  /// @param volumeTriple is an (optional) triplet of volumes
  TrackingVolumePtr
  trackingVolume(TrackingVolumePtr   insideVolume  = nullptr,
                 VolumeBoundsPtr     outsideBounds = nullptr,
                 const LayerTriple*  layerTriple   = nullptr,
                 const VolumeTriple* volumeTriple  = nullptr) const override;

  /// Set configuration method 
  /// @param cvbConfig is the new configuration to be set               
  void
  setConfiguration(const Config& cvbConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

private:
  /// Configuration struct 
  Config m_cfg;

  /// Private access to the logger
  const Logger&
  logger() const
  {
    return *m_cfg.logger;
  }
  
  /// analyse the layer setup 
  LayerSetup
  analyzeLayerSetup(const LayerVector lVector) const;
};

/// Return the configuration object 
inline CylinderVolumeBuilder::Config
CylinderVolumeBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H
