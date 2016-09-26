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

/// VolumeSetup struct to understand the layer setup
struct VolumeSetup
{
  bool        present;   ///< layers are present
  bool        wrapping;  ///< in what way they are binned
  double      rMin;      ///< min/max parameters
  double      rMax;      ///< min/max parameters
  double      zMin;      ///< min/max parameters
  double      zMax;      ///< min/max parameters
  LayerVector layers;    ///< the layers you have

  /// Default constructor
  VolumeSetup()
    : present(false)
    , wrapping(false)
    , rMin(10e10)
    , rMax(10e-10)
    , zMin(10e10)
    , zMax(-10e10)
    , layers()
  {
  }

  /// Adapt to the dimensions of another setup
  ///
  /// @param lSetup is the setup to which it should be adapded
  void
  adapt(const VolumeSetup& lSetup)
  {
    takeSmaller(rMin, lSetup.rMin);
    takeBigger(rMax, lSetup.rMax);
    takeSmaller(zMin, lSetup.zMin);
    takeBigger(zMax, lSetup.zMax);
  }

  /// Overlap check radially
  bool
  overlapsInR(const VolumeSetup& vSetup) const
  {
    if (!present) return false;
    if (rMin < vSetup.rMax && rMax > vSetup.rMax) return true;
    if (vSetup.rMin < rMax && vSetup.rMax > rMin) return true;
    return false;
  }

  /// Overlap check longitudinally
  bool
  overlapsInZ(const VolumeSetup& vSetup) const
  {
    if (!present) return false;
    if (zMin < vSetup.zMax && zMax > vSetup.zMax) return true;
    if (vSetup.zMin < zMax && vSetup.zMax > zMin) return true;
    return false;
  }

  /// Compatibility check radially
  bool
  wrapsInR(const VolumeSetup& vSetup) const
  {
    if (vSetup.rMax > rMin) return false;
    return true;
  }

  /// Compatibility check longitudinally
  bool
  wrapsInZ(const VolumeSetup& vSetup) const
  {
    if (vSetup.zMin < zMin || vSetup.zMax > zMax) return false;
    return true;
  }

  /// Compatibility check full set
  bool
  wraps(const VolumeSetup& vSetup) const
  {
    // it wraps
    return wrapsInR(vSetup) && wrapsInZ(vSetup);
  }

  /// toString
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

/// @class CylinderVolumeBuilder
///
/// A volume builder to be used for building a concentrical cylindrical volumes
///  - a) configured volume
///  - b) wrapping around a cylindrical/disk layer setup
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
    std::shared_ptr<ITrackingVolumeHelper> trackingVolumeHelper = nullptr;
    /// the string based indenfication
    std::string volumeName = "";
    /// The dimensions of the manually created world
    std::vector<double> volumeDimension = {};
    /// the world material
    std::shared_ptr<Material> volumeMaterial = nullptr;
    /// build the volume to the beam line
    bool buildToRadiusZero = false;
    /// needed to build layers within the volume
    std::shared_ptr<ILayerBuilder> layerBuilder = nullptr;
    /// the envelope covering the potential layers rMin, rMax
    std::pair<double, double> layerEnvelopeR
        = {5. * Acts::units::_mm, 5. * Acts::units::_mm};
    /// the envelope covering the potential layers inner/outer
    double layerEnvelopeZ = 10. * Acts::units::_mm;
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
  /// @param cvbConfig is the configuraiton struct to steer the builder
  /// @param logger logging instance
  CylinderVolumeBuilder(const Config&           cvbConfig,
                        std::unique_ptr<Logger> logger
                        = getDefaultLogger("CylinderVolumeBuilder",
                                           Logging::INFO));

  /// Destructor
  virtual ~CylinderVolumeBuilder();

  /// CylinderVolumeBuilder main call method
  ///
  /// @param insideVolume is an (optional) volume to be wrapped
  /// @param outsideBounds is an (optional) outside confinement
  /// @param layerTriple is an (optional) triplet of layers
  /// @param volumeTriple is an (optional) triplet of volumes
  TrackingVolumePtr
  trackingVolume(TrackingVolumePtr  insideVolume  = nullptr,
                 VolumeBoundsPtr    outsideBounds = nullptr,
                 const LayerTriple* layerTriple   = nullptr) const override;

  /// Set configuration method
  /// @param cvbConfig is the new configuration to be set
  void
  setConfiguration(const Config& cvbConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

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
  std::unique_ptr<Logger> m_logger;

  /// Analyze the layer setup to gather needed dimension
  ///
  /// @param lVector is the vector of layers that are parsed
  /// @return a VolumeSetup representing this layer
  VolumeSetup
  analyzeLayers(const LayerVector& lVector) const;

  /// Helper method check the layer containment,
  /// both for inside / outside.
  ///
  /// @param layerSetup is the VolumeSetup to be tested
  ///        the wrapping flag may be set
  /// @param insideSetup is the inside volume in order to
  ///        check the wrapping
  /// @param volumeSetup is the volume to be tested
  /// @param sign distinguishes inside/outside testing
  /// @return boolean that indicates the test result
  bool
  checkLayerContainment(VolumeSetup&       layerSetup,
                        const VolumeSetup& insideSetup,
                        const VolumeSetup& volumeSetup,
                        int                sign) const;

  /// Wynchronize the layer setups with given
  /// inside / outside constraints.
  ///
  /// This is the last method to be called in the building
  /// chain. It adapts the setups accordingly and sets the
  /// right boundaries.
  ///
  /// @param nSetup the setup of negative EC layers (if present)
  /// @param cSetup the setup of the barrel layers
  /// @param pSetup the setup of positive EC layers (if present)
  /// @param is the inside volume setup/dimensions
  /// @param is the outside and final volume setup/dimensions
  ///
  /// @note non-const references may be changed
  ///
  /// @return a wrapping condition @TODO check if needed
  WrappingCondition
  synchronizeVolumeSetups(VolumeSetup&       nSetup,
                          VolumeSetup&       cSetup,
                          VolumeSetup&       pSetup,
                          const VolumeSetup& insideSetup,
                          VolumeSetup&       outsideBoundSetup) const;
};

/// Return the configuration object
inline CylinderVolumeBuilder::Config
CylinderVolumeBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#endif  // ACTS_TOOLS_CYLINDERVOLUMEBUILDER_H
