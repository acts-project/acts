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

#include <string>
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

/// VolumeConfig struct to understand the layer config
struct VolumeConfig
{
  bool        present;   ///< layers are present
  bool        wrapping;  ///< in what way they are binned
  double      rMin;      ///< min/max parameters
  double      rMax;      ///< min/max parameters
  double      zMin;      ///< min/max parameters
  double      zMax;      ///< min/max parameters
  LayerVector layers;    ///< the layers you have

  /// Default constructor
  VolumeConfig()
    : present(false)
    , wrapping(false)
    , rMin(10e10)
    , rMax(10e-10)
    , zMin(10e10)
    , zMax(-10e10)
    , layers()
  {
  }

  /// Adapt to the dimensions of another config
  ///
  /// @param lConfig is the config to which it should be adapded
  void
  adapt(const VolumeConfig& lConfig)
  {
    takeSmaller(rMin, lConfig.rMin);
    takeBigger(rMax, lConfig.rMax);
    takeSmaller(zMin, lConfig.zMin);
    takeBigger(zMax, lConfig.zMax);
  }

  /// Overlap check radially
  ///
  /// @param vConfig is the config against which is checked
  bool
  overlapsInR(const VolumeConfig& vConfig) const
  {
    if (!present) return false;
    return std::max(rMin, vConfig.rMin) <= std::min(rMax, vConfig.rMax);
  }

  /// Overlap check longitudinally
  ///
  /// @param vConfig is the config against which is checked
  bool
  overlapsInZ(const VolumeConfig& vConfig) const
  {
    if (!present) return false;
    return std::max(zMin, vConfig.zMin) <= std::min(zMax, vConfig.zMax);
  }

  /// Compatibility check full set
  ///
  /// @param vConfig is the config against which is checked
  bool
  wraps(const VolumeConfig& vConfig) const
  {
    if ((zMax <= vConfig.zMin) || (zMin >= vConfig.zMax)) return true;
    return containesInR(vConfig);
  }

  /// Check if contained full set
  ///
  /// @param vConfig is the config against which is checked
  bool
  containes(const VolumeConfig& vConfig) const
  {
    return (containesInR(vConfig) && containesInZ(vConfig));
  }

  /// Check if contained radially
  ///
  /// @param vConfig is the config against which is checked
  bool
  containesInR(const VolumeConfig& vConfig) const
  {
    return (rMin >= vConfig.rMax);
  }

  /// Check if contained longitudinally
  ///
  /// @param vConfig is the config against which is checked
  bool
  containesInZ(const VolumeConfig& vConfig) const
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
  ///| Negative Endcap | Barrel | Positive Endcap |	=> four boundaries needed in
  /// z
  ///                  | Barrel |				    => two boundaries needed in z
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
///  All are optionally wrapped around a given volume which has to by a cylinder
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

  /// @enum WrappingCondition
  enum WrappingCondition {

    SynchronizationError = 0,
    NoWrapping           = 1,
    BarrelWrapping       = 2,
    BarrelWrappingGaps   = 3,
    TripleWrapping       = 4,
    TripleWrappingGaps   = 5,
    Undefined            = 6

  };

  /// Constructor
  ///
  /// @param cvbConfig is the configuraiton struct to steer the builder
  /// @param logger logging instance
  CylinderVolumeBuilder(const Config&                 cvbConfig,
                        std::unique_ptr<const Logger> logger
                        = getDefaultLogger("CylinderVolumeBuilder",
                                           Logging::INFO));

  /// Destructor
  virtual ~CylinderVolumeBuilder();

  /// CylinderVolumeBuilder main call method
  ///
  /// @param insideVolume is an (optional) volume to be wrapped
  /// @param outsideBounds is an (optional) outside confinement
  MutableTrackingVolumePtr
  trackingVolume(TrackingVolumePtr insideVolume  = nullptr,
                 VolumeBoundsPtr   outsideBounds = nullptr) const override;

  /// Set configuration method
  ///
  /// @param cvbConfig is the new configuration to be set
  void
  setConfiguration(const Config& cvbConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// set logging instance
  ///
  /// @param logger is the logging istance to be set
  void
  setLogger(std::unique_ptr<const Logger> logger);

private:
  /// Configuration struct
  Config m_cfg;

  /// Private access to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Analyze the layer config to gather needed dimension
  ///
  /// @param lVector is the vector of layers that are parsed
  /// @return a VolumeConfig representing this layer
  VolumeConfig
  analyzeLayers(const LayerVector& lVector) const;

  /// Helper method check the layer containment,
  /// both for inside / outside.
  ///
  /// @param layerConfig is the VolumeConfig to be tested
  ///        the wrapping flag may be set
  /// @param insideConfig is the inside volume in order to
  ///        check the wrapping
  /// @param volumeConfig is the volume to be tested
  /// @param sign distinguishes inside/outside testing
  ///
  /// @return boolean that indicates the test result
  bool
  checkLayerContainment(VolumeConfig&       layerConfig,
                        const VolumeConfig& insideConfig,
                        const VolumeConfig& volumeConfig,
                        int                 sign) const;

  /// Synchronize the layer configs with given
  /// inside / outside constraints.
  ///
  /// This is the last method to be called in the building
  /// chain. It adapts the configs accordingly and sets the
  /// right boundaries.
  ///
  /// @param nConfig the config of negative EC layers (if present)
  /// @param cConfig the config of the barrel layers
  /// @param pConfig the config of positive EC layers (if present)
  /// @param insideConfig is the inside volume config/dimensions
  /// @param outsideBoundConfig is the outside and final volume
  /// config/dimensions
  ///
  /// @note non-const references may be changed
  ///
  /// @return a wrapping condition @todo check if needed
  WrappingCondition
  synchronizeVolumeConfigs(VolumeConfig&       nConfig,
                           VolumeConfig&       cConfig,
                           VolumeConfig&       pConfig,
                           const VolumeConfig& insideConfig,
                           VolumeConfig&       outsideBoundConfig) const;
};

/// Return the configuration object
inline CylinderVolumeBuilder::Config
CylinderVolumeBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H
