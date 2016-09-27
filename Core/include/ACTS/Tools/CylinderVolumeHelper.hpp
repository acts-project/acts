// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_CYLINDERVOLUMEHELPER_H
#define ACTS_TOOLS_CYLINDERVOLUMEHELPER_H 1

#ifndef ACTS_TOOLS_TAKESMALLERBIGGER
#define ACTS_TOOLS_TAKESMALLERBIGGER
#define takeSmaller(current, test) current = current < test ? current : test
#define takeBigger(current, test) current  = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test)                           \
  takeSmaller(cSmallest, test);                                                \
  takeBigger(cBiggest, test)
#endif

#include <memory>
#include <string>
#include <vector>
#include "ACTS/Tools/ILayerArrayCreator.hpp"
#include "ACTS/Tools/ITrackingVolumeArrayCreator.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"

namespace Acts {

class Layer;
class TrackingVolume;
class VolumeBounds;
class CylinderVolumeBounds;
class Material;

/// @class CylinderVolumeHelper
///
/// The concrete implementation for cylindrical TrackingVolume
/// objects of the ITrackingVolumeCreator interface
///
class CylinderVolumeHelper : public ITrackingVolumeHelper
{
public:
  /// @struct Config
  /// Nested configuration struct for this CylinderVolumeHelper
  struct Config
  {
    /// a tool for coherent LayerArray creation
    std::shared_ptr<ILayerArrayCreator> layerArrayCreator = nullptr;
    /// Helper Tool to create TrackingVolume
    std::shared_ptr<ITrackingVolumeArrayCreator> trackingVolumeArrayCreator
        = nullptr;
    /// thickness of passive layers
    double passiveLayerThickness = 1;
    /// bins in phi for the passive layer
    int passiveLayerPhiBins = 1;
    /// bins in r/z for the passive layer
    int passiveLayerRzBins = 100;
  };

  /// Constructor
  /// @param cvhConfig is the configuration struct for this builder
  /// @param logger logging instance
  CylinderVolumeHelper(const Config&           cvhConfig,
                       std::unique_ptr<Logger> logger
                       = getDefaultLogger("CylinderVolumeHelper",
                                          Logging::INFO));

  /// Destructor
  virtual ~CylinderVolumeHelper() = default;

  // clang-format off
  /// @copydoc ITrackingVolumeHelper::createTrackingVolume(const LayerVector&,std::shared_ptr<Material>,VolumeBoundsPtr,std::shared_ptr<Transform3D>,const std::string&,BinningType) const
  // clang-format on
  TrackingVolumePtr
  createTrackingVolume(const LayerVector&           layers,
                       std::shared_ptr<Material>    matprop,
                       VolumeBoundsPtr              volBounds,
                       std::shared_ptr<Transform3D> transform = nullptr,
                       const std::string& volumeName = "UndefinedVolume",
                       BinningType        btype      = arbitrary) const;
  // clang-format off
  /// @copydoc ITrackingVolumeHelper::createTrackingVolume(const LayerVector&,std::shared_ptr<Material>,double,double,double,double,const std::string&,BinningType) const
  // clang-format on
  TrackingVolumePtr
  createTrackingVolume(const LayerVector&        layers,
                       std::shared_ptr<Material> matprop,
                       double                    loc1Min,
                       double                    loc1Max,
                       double                    loc2Min,
                       double                    loc2Max,
                       const std::string&        volumeName = "UndefinedVolume",
                       BinningType               btype      = arbitrary) const;
  // clang-format off
  /// @copydoc ITrackingVolumeHelper::createGapTrackingVolume(std::shared_ptr<Material>,double,double,double,double,unsigned int,bool,const std::string&) const
  // clang-format on
  /// @note for cylindrical implementation loc1Min = rMin, loc1Max = rMax,
  /// loc2Min = zMin, loc2Max = zMax
  TrackingVolumePtr
  createGapTrackingVolume(std::shared_ptr<Material> matprop,
                          double                    loc1Min,
                          double                    loc1Max,
                          double                    loc2Min,
                          double                    loc2Max,
                          unsigned int              materialLayers,
                          bool                      cylinder = true,
                          const std::string&        volumeName
                          = "UndefinedVolume") const;
  // clang-format off
  /// @copydoc ITrackingVolumeHelper::createGapTrackingVolume(std::shared_ptr<Material>  matprop,double,double,double,double,const std::vector<double>&,bool,const std::string&,BinningType) const
  // clang-format on
  /// @note for cylindrical implementation loc1Min = rMin, loc1Max = rMax,
  /// loc2Min = zMin, loc2Max = zMax
  TrackingVolumePtr
  createGapTrackingVolume(std::shared_ptr<Material>  matprop,
                          double                     loc1Min,
                          double                     loc1Max,
                          double                     loc2Min,
                          double                     loc2Max,
                          const std::vector<double>& layerPositions,
                          bool                       cylinder = true,
                          const std::string& volumeName = "UndefinedVolume",
                          BinningType        btype      = arbitrary) const;

  /// Create a container volumes from sub volumes, input volumes are ordered in
  /// R or Z by convention */
  TrackingVolumePtr
  createContainerTrackingVolume(const TrackingVolumeVector& volumes) const;

  /// Set configuration method
  void
  setConfiguration(const Config& cvbConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

protected:
  /// Configuration object
  Config m_cfg;

private:
  /// Private access method to the logging instance
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// the looging instance
  std::unique_ptr<Logger> m_logger;

  /// Private method - it estimates the CylinderBounds and Translation
  /// of layers, these are checked against the layer positions/dimensions.
  bool
  estimateAndCheckDimension(const LayerVector&                 layers,
                            const Acts::CylinderVolumeBounds*& cylBounds,
                            std::shared_ptr<Transform3D>&      transform,
                            double&                            rMinClean,
                            double&                            rMaxClean,
                            double&                            zMinClean,
                            double&                            zMaxClean,
                            BinningValue&                      bValue,
                            BinningType bType = arbitrary) const;

  /// Private method - interglue all volumes contained by a TrackingVolume
  /// and set the outside glue volumes in the descriptor
  bool
  interGlueTrackingVolume(TrackingVolumePtr tVolume,
                          bool              rBinned,
                          double            rMin,
                          double            rGlueMin,
                          double            rMax,
                          double            zMin,
                          double            zMax) const;

  /// Private method - glue volume to the other -- use trackingVolume helper
  void
  glueTrackingVolumes(TrackingVolumePtr   volumeOne,
                      BoundarySurfaceFace faceOne,
                      TrackingVolumePtr   volumeTwo,
                      BoundarySurfaceFace faceTwod,
                      double              rMin,
                      double              rGlueMin,
                      double              rMax,
                      double              zMin,
                      double              zMax) const;

  /// Private method - helper method not to duplicate code
  void
  addFaceVolumes(TrackingVolumePtr         tVolume,
                 Acts::BoundarySurfaceFace bsf,
                 TrackingVolumeVector&     vols) const;

  /// Private method - helper method to save some code
  LayerPtr
  createCylinderLayer(double z,
                      double r,
                      double halflength,
                      double thickness,
                      int    binsPhi,
                      int    binsZ) const;

  /// Private method - helper method to save some code
  LayerPtr
  createDiscLayer(double z,
                  double rMin,
                  double rMax,
                  double thickness,
                  int    binsPhi,
                  int    binsR) const;
};

inline CylinderVolumeHelper::Config
CylinderVolumeHelper::getConfiguration() const
{
  return m_cfg;
}
}

#endif  // ACTS_TOOLS_CYLINDERVOLUMEHELPER_H
