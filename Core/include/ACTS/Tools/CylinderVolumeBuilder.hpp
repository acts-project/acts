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
#define ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H

#include <string>
#include <limits>
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeBuilder.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/Units.hpp"

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

/// @enum WrappingCondition
enum WrappingCondition {
  SynchronizationError = 0, ///< inconsistency detected
  Undefined            = 0, ///< inconsistency detected
  Attaching            = 1, ///< attach the volumes
  Inserting            = 2, ///< insert the new volume
  Wrapping             = 3  ///< wrap the new volume
};

/// VolumeConfig struct to understand the layer config
struct VolumeConfig
{
  bool        present;   ///< layers are present
  bool        wrapping;  ///< in what way they are binned
  double      rMin;      ///< min parameter r
  double      rMax;      ///< max parameter r
  double      zMin;      ///< min parameter z
  double      zMax;      ///< max parameter z
  LayerVector layers;    ///< the layers you have

  /// Default constructor
  VolumeConfig()
    : present(false)
    , wrapping(false)
    , rMin(std::numeric_limits<double>::max())
    , rMax(std::numeric_limits<double>::lowest())
    , zMin(std::numeric_limits<double>::max())
    , zMax(std::numeric_limits<double>::lowest())
    , layers()
  {}

  /// Adapt to the dimensions of another config in Z
  /// it will take the maximum/minimum values and just overwrite them 
  ///
  /// @param [in] lConfig is the config to which it should be adapded
  void
  adaptZ(const VolumeConfig& lConfig)
  {
    if (lConfig){
      takeSmaller(zMin, lConfig.zMin);
      takeBigger(zMax, lConfig.zMax);
    }
  }


  /// Adapt to the dimensions of another config in R
  /// it will take the maximum/minimum values and just overwrite them 
  ///
  /// @param [in] lConfig is the config to which it should be adapded
  void
  adaptR(const VolumeConfig& lConfig)
  {
    if (lConfig){
      takeSmaller(rMin, lConfig.rMin);
      takeBigger(rMax, lConfig.rMax);
    }
  }

  /// Adapt to the dimensions of another config
  /// it will take the maximum/minimum values and just overwrite them 
  ///
  /// @param [in] lConfig is the config to which it should be adapded
  void
  adapt(const VolumeConfig& lConfig)
  {
    adaptZ(lConfig);
    adaptR(lConfig);
  }

  /// Attach method - non-const
  /// it attaches the one volume config to the other one in Z
  /// this is the non-cost method, i.e. the mit point is is used
  /// 
  /// @param [in] lConfig is the config to which it should be attached
  /// @note lConfig will be changed 
  void midPointAttachZ(VolumeConfig& lConfig){
    if (lConfig.zMin >= zMax){
      double zMid = 0.5*(lConfig.zMin+zMax);
      lConfig.zMin = zMid;
      zMax = zMid;
    } else {
      double zMid = 0.5*(zMin+lConfig.zMax);
      lConfig.zMax = zMid;
      zMin = zMid;  
    }
  }

  /// Attach method - const
  /// it attaches the one volume config to the other one
  /// 
  /// @param [in] lConfig is the confit to which it should be attached
  void attachZ(const VolumeConfig& lConfig){
    if (lConfig.zMin >= zMax)
       zMax = lConfig.zMin;
     else 
      zMin = lConfig.zMax;  
  }

  /// Overlap check radially
  ///
  /// @param [in] vConfig is the config against which is checked
  /// @return boolean if the overlap in r exists
  bool
  overlapsInR(const VolumeConfig& vConfig) const
  {
    if (!present) return false;
    return std::max(rMin, vConfig.rMin) <= std::min(rMax, vConfig.rMax);
  }

  /// Overlap check longitudinally
  ///
  /// @param [in] vConfig is the config against which is checked
  /// @return boolean if the overlap in z exists
  bool
  overlapsInZ(const VolumeConfig& vConfig) const
  {
    if (!present) return false;
    return std::max(zMin, vConfig.zMin) <= std::min(zMax, vConfig.zMax);
  }

  /// Compatibility check full set
  ///
  /// @param [in] vConfig is the config against which is checked
  /// @return boolean if the current volume wraps the vConfig fully
  bool
  wraps(const VolumeConfig& vConfig) const
  {
    if ((zMax <= vConfig.zMin) || (zMin >= vConfig.zMax)) return true;
    return containsInR(vConfig);
  }

  /// Check if contained full set
  ///
  /// @param [in] vConfig is the config against which is checked
  bool
  contains(const VolumeConfig& vConfig) const
  {
    return (containsInR(vConfig) && containsInZ(vConfig));
  }

  /// Check if contained radially
  ///
  /// @param [in] vConfig is the config against which is checked
  bool
  containsInR(const VolumeConfig& vConfig) const
  {
    return (rMin >= vConfig.rMax);
  }

  /// Check if contained longitudinally
  ///
  /// @param [in] vConfig is the config against which is checked
  bool
  containsInZ(const VolumeConfig& vConfig) const
  {
    return (vConfig.zMin > zMin && vConfig.zMax < zMax);
  }

  /// Method for output formatting
  std::string
  toString() const
  {
    /// for screen output
    std::stringstream sl;
    sl << rMin << ", " << rMax << " / " << zMin << ", " << zMax;
    return sl.str();
  }

  /// Conversion operator to bool
  operator bool() const { return present; }
};

/// @brief The WrappingSetup that is happening here
struct WrappingConfig 
{
  public:
    /// the new volumes 
    VolumeConfig nVolumeConfig;
    VolumeConfig cVolumeConfig;
    VolumeConfig pVolumeConfig;
      
    /// the combined volume 
    VolumeConfig containerVolumeConfig;

    /// existing volume config with potential gaps
    VolumeConfig existingVolumeConfig;
    VolumeConfig fGapVolumeConfig;
    VolumeConfig sGapVolumeConfig;
    
    /// externally provided config, this can only change the 
    /// the ncp volumes 
    VolumeConfig externalVolumeConfig;
    
    // WrappingCondition
    WrappingCondition wCondition = Undefined;
    
    /// constructor
    WrappingConfig() {}
    
    /// configure the new Volume 
    void configureContainerVolume() {
      containerVolumeConfig.present = true;
      // adapt the new volume config to the existing configs
      if (nVolumeConfig) containerVolumeConfig.adapt(nVolumeConfig);
      if (cVolumeConfig) containerVolumeConfig.adapt(cVolumeConfig);
      if (pVolumeConfig) containerVolumeConfig.adapt(pVolumeConfig);
      // adapt the external one
      if (externalVolumeConfig)
        containerVolumeConfig.adapt(externalVolumeConfig);
      // attach the volume configs 
      if (nVolumeConfig && cVolumeConfig) 
        nVolumeConfig.midPointAttachZ(cVolumeConfig);
      if (cVolumeConfig && pVolumeConfig) 
        cVolumeConfig.midPointAttachZ(pVolumeConfig);
      // adapt r afterwards
      nVolumeConfig.adaptR(containerVolumeConfig);
      cVolumeConfig.adaptR(containerVolumeConfig);
      pVolumeConfig.adaptR(containerVolumeConfig);      
    }
    
    /// wrap, insert, attach
    void wrapInsertAttach(){
      // action is only needed if an existing volume
      // is present
      if (existingVolumeConfig){
        // 0 - simple atachment case
        if (!cVolumeConfig){
          // check if it can be easily attached 
          if (nVolumeConfig 
              && nVolumeConfig.zMax < existingVolumeConfig.zMin){
            nVolumeConfig.attachZ(existingVolumeConfig);         
            // will attach the new volume(s)
            wCondition = Attaching;
          }  
          if (pVolumeConfig 
              && pVolumeConfig.zMin < existingVolumeConfig.zMax){
            pVolumeConfig.attachZ(existingVolumeConfig);  
            // will attach the new volume(s)
            wCondition = Attaching;
          }
          // see if inner glue volumes are needed 
          if (containerVolumeConfig.rMin > existingVolumeConfig.rMin){
              nVolumeConfig.rMin = existingVolumeConfig.rMin;
              pVolumeConfig.rMin = existingVolumeConfig.rMin;
          } else {
              fGapVolumeConfig.present = true;
              // get the zMin/zMax boundaries
              fGapVolumeConfig.adaptZ(existingVolumeConfig);
              fGapVolumeConfig.rMin = containerVolumeConfig.rMin;
              fGapVolumeConfig.rMax = existingVolumeConfig.rMin;
          }
          // see if outer glue volumes are needed
          if (containerVolumeConfig.rMax < existingVolumeConfig.rMax){
            nVolumeConfig.rMax = existingVolumeConfig.rMax;
            pVolumeConfig.rMax = existingVolumeConfig.rMax;
          } else {
            sGapVolumeConfig.present = true;
            // get the zMin/zMax boundaries
            sGapVolumeConfig.adaptZ(existingVolumeConfig);
            sGapVolumeConfig.rMin = existingVolumeConfig.rMax;
            sGapVolumeConfig.rMax = containerVolumeConfig.rMax;
          }
        } else {
          // wrapping or insertion case 
          if (existingVolumeConfig.rMax < containerVolumeConfig.rMin){
            // wrapping case
            nVolumeConfig.rMin = existingVolumeConfig.rMax;
            cVolumeConfig.rMin = existingVolumeConfig.rMax;
            pVolumeConfig.rMin = existingVolumeConfig.rMax;
            // will wrap the new volume(s) around existing
            wCondition = Wrapping;
          } else if (existingVolumeConfig.rMin > containerVolumeConfig.rMax){
            // insertion case
            nVolumeConfig.rMax = existingVolumeConfig.rMin;
            cVolumeConfig.rMax = existingVolumeConfig.rMin;
            pVolumeConfig.rMax = existingVolumeConfig.rMin;
            // will insert the new volume(s) into existing
            wCondition = Inserting;
          }
          // check if gaps are needed - at negative sector
          if (existingVolumeConfig.zMin > containerVolumeConfig.zMin){
            fGapVolumeConfig.present = true;
            fGapVolumeConfig.adaptR(existingVolumeConfig);
            fGapVolumeConfig.zMin = containerVolumeConfig.zMin;
            fGapVolumeConfig.zMax = existingVolumeConfig.zMin;
          } else {
            // adapt lower z boundary
            if (nVolumeConfig) nVolumeConfig.zMin = existingVolumeConfig.zMin;
            else if (cVolumeConfig) cVolumeConfig.zMin = existingVolumeConfig.zMin;
          }
          // check if gaps are needed - at positive sector
          if (existingVolumeConfig.zMax < containerVolumeConfig.zMax){
            sGapVolumeConfig.present = true;
            sGapVolumeConfig.adaptR(existingVolumeConfig);
            sGapVolumeConfig.zMin = existingVolumeConfig.zMax;
            sGapVolumeConfig.zMax = containerVolumeConfig.zMax;
          } else {
            // adapt higher z boundary
            if (pVolumeConfig) pVolumeConfig.zMax = existingVolumeConfig.zMax;
            else if (cVolumeConfig) cVolumeConfig.zMax = existingVolumeConfig.zMax;
          }
        }
      }
      return;
    }
    
    /// Method for output formatting
    std::string
    toString() const
    {
      // for screen output
      std::stringstream sl;
      if (containerVolumeConfig)
         sl << "New contaienr built with       configuration: " 
           << containerVolumeConfig.toString() << '\n';
      // go throug the new new ones first
      if (nVolumeConfig) 
         sl << " - n: Negative Endcap, current configuration: "  
            << nVolumeConfig.toString() << '\n';
      if (cVolumeConfig) 
         sl << " - c: Barrel, current          configuration: " 
            << cVolumeConfig.toString() << '\n';
      if (pVolumeConfig) 
         sl << " - p: Negative Endcap, current configuration: " 
            << pVolumeConfig.toString() << '\n';
      if (existingVolumeConfig){
         sl << "Existing volume with           configuration: " 
           << existingVolumeConfig.toString() << '\n';
        if (fGapVolumeConfig)
          sl << " - g1: First gap volume,       configuration : " 
            << fGapVolumeConfig.toString() << '\n';
        if (sGapVolumeConfig)
          sl << " - g2: Second gap volume,      configuration : " 
            << sGapVolumeConfig.toString() << '\n';
        if (wCondition != Undefined){
          sl << "WrappingCondition = " << wCondition << '\n';
          
        }
      }      
      return sl.str();
    }
      
};


/// @brief The sub volume struct gives the possibility to hand over the
/// dimensions of the sub volumes (i.e. barrel,endcaps)
///
/// For the sub volume config the minimum and maximum radius of all the sub
/// volumes (i.e. barrel,endcaps) need to be set. Furthermore two (if only
/// barrel layers are given) or four boundaries (if barrel and endcap layers are
/// given) can be set.
struct SubVolumeConfig
{
  /// the minimum radius of the sub volume config in the central barrel region
  double centralRmin;
  /// the minimum radius of the sub volume config in the outer endcap region
  /// which can differ from the barrel rmin and be smaller
  double outerRmin;
  /// the maximum radius of the sub volume config - for both regions use the
  /// bigger one
  double rMax;
  /// possible cases for zBounadries:
  ///| Negative Endcap | Barrel | Positive Endcap |	=> four z boundaries needed
  ///                  | Barrel |				            => two z boundaries needed 
  /// @note the zBoundaries need to be handed over sorted in ascending order
  std::vector<double> zBoundaries;

  /// Default constructor
  /// Per default the SubVolumeConfig is not set
  SubVolumeConfig() : centralRmin(-1.), outerRmin(-1.), rMax(-1.), zBoundaries()
  {
  }

  /// Conversion operator to bool needed for checks if the sub volume config is
  /// set by the user
  operator bool() const
  {
    return (zBoundaries.size() > 1 && centralRmin >= 0. && rMax > 0.
            && centralRmin < rMax
            && outerRmin < rMax);
  }
};

/// @class CylinderVolumeBuilder
///
/// A volume builder to be used for building a concentrical cylindrical volumes
///  - a) configured volume
///  - b) wrapping around a cylindrical/disk layer config
///
/// All are optionally wrapped around a given volume which has to by a cylinder
/// volume and which has to be center at z == 0
///
///  To receive the tracking volume it is possible to also hand over a triple of
/// layers, which is a C++ tuple of three pointers to layer vectors (defined in
/// the ITrackingVolumeBuilder). This functionality is needed for a possible
/// translation of an geometry existing in another format. The first entry
/// represents the layers of the negative endcap, the second the layers of the
/// barrel and the third the layers of the positive endcap. If the one of these
/// pointers is a nullptr no layers will be created for this volume

class CylinderVolumeBuilder : public ITrackingVolumeBuilder
{
public:
  /// @struct Config
  /// Nested configuration struct for this CylinderVolumeBuilder
  struct Config
  {
    /// the trackign volume helper for construction
    std::shared_ptr<const ITrackingVolumeHelper> trackingVolumeHelper = nullptr;
    /// the string based indenfication
    std::string volumeName = "";
    /// The dimensions of the manually created world
    std::vector<double> volumeDimension = {};
    /// the world material
    std::shared_ptr<const Material> volumeMaterial = nullptr;
    /// build the volume to the beam line
    bool buildToRadiusZero = false;
    /// needed to build layers within the volume if no SubVolumeConfig is given
    std::shared_ptr<const ILayerBuilder> layerBuilder = nullptr;
    /// the envelope covering the potential layers rMin, rMax if no
    /// SubVolumeConfig is given
    std::pair<double, double> layerEnvelopeR
        = {5. * Acts::units::_mm, 5. * Acts::units::_mm};
    /// the envelope covering the potential layers inner/outer
    double layerEnvelopeZ = 10. * Acts::units::_mm;
    /// possible SubVoumeConfig - it can provide the boundaries for the sub
    /// volumes and replaces the layerEnvelope configuration in that case
    SubVolumeConfig subVolumeConfig;
    /// the volume signature
    int volumeSignature = -1;
  };

  /// Constructor
  ///
  /// @param [in] cvbConfig is the configuraiton struct to steer the builder
  /// @param [in] logger logging instance
  CylinderVolumeBuilder(const Config&                 cvbConfig,
                        std::unique_ptr<const Logger> logger
                        = getDefaultLogger("CylinderVolumeBuilder",
                                           Logging::INFO));

  /// Destructor
  virtual ~CylinderVolumeBuilder();

  /// CylinderVolumeBuilder main call method
  ///
  /// @param [in] existingVolume is an (optional) volume to be included
  /// @param [in] externalBounds are (optional) external confinement
  ///             constraints
  /// @return a mutable pointer to a new TrackingVolume which includes the 
  ///         optionally provided exisitingVolume consistently for further
  ///         processing
  MutableTrackingVolumePtr
  trackingVolume(TrackingVolumePtr existingVolume  = nullptr,
                 VolumeBoundsPtr   externalBounds = nullptr) const override;

  /// Set configuration method
  ///
  /// @param [in] cvbConfig is the new configuration to be set
  void
  setConfiguration(const Config& cvbConfig);

  /// Get configuration method
  /// 
  /// @return a copy of the config object
  Config
  getConfiguration() const;

  /// set logging instance
  ///
  /// @param [in] logger is the logging istance to be set
  void
  setLogger(std::unique_ptr<const Logger> logger);

private:
  /// Configuration struct
  Config m_cfg;

  /// Private access to the logger
  ///
  /// @return a const reference to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Analyze the layer config to gather needed dimension
  ///
  /// @param [in] lVector is the vector of layers that are parsed
  /// @return a VolumeConfig representing this layer
  VolumeConfig
  analyzeLayers(const LayerVector& lVector) const;

  /// Helper method check the layer containment,
  /// both for inside / outside.
  ///
  /// @param [in] layerConfig is the VolumeConfig to be tested
  ///        the wrapping flag may be set
  /// @param [in] insideConfig is the inside volume in order to
  ///        check the wrapping
  /// @param [in] volumeConfig is the volume to be tested
  /// @param [in] sign distinguishes inside/outside testing
  ///
  /// @return boolean that indicates the test result
  bool
  checkLayerContainment(VolumeConfig&       layerConfig,
                        const VolumeConfig& insideConfig,
                        const VolumeConfig& volumeConfig,
                        int                 sign) const;

};

/// Return the configuration object
inline CylinderVolumeBuilder::Config
CylinderVolumeBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H
