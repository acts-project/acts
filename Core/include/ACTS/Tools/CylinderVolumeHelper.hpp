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


  /// Create a TrackingVolume* from a set of layers and (optional) parameters
  ///
  /// @param layers vector of static layers confined by the TrackingVolume
  /// if no bounds or HepTransform is given, they define the size
  /// together with the volume enevlope parameters
  /// @param matprop dense material properties for this TrackingVolume
  /// @param volBounds (optional) bounds of this TrackingVolume - ownership
  /// given
  /// @param transform (optional) placement of this TrackingVolume - ownership
  /// given
  /// @param volumeName  volume name to be given
  /// @param btype (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  TrackingVolumePtr
  createTrackingVolume(const LayerVector&           layers,
                       std::shared_ptr<Material>    matprop,
                       VolumeBoundsPtr              volBounds,
                       std::shared_ptr<Transform3D> transform = nullptr,
                       const std::string& volumeName = "UndefinedVolume",
                       BinningType        btype      = arbitrary) const;
                       
  /// Create a TrackingVolume* from a set of layers and (optional) parameters
  ///
  /// @param layers vector of static layers confined by the TrackingVolume
  /// if no bounds or HepTransform is given, they define the size
  /// together with the volume enevlope parameters
  /// @param matprop dense material properties for this TrackingVolume
  /// @param loc1Min, loc1Max, loc2Min, loc2Max : local position in space,
  /// this TrackingVolume is restricted to Translation only
  /// @param volumeName  volume name to be given
  /// @param btype (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
  TrackingVolumePtr
  createTrackingVolume(const LayerVector&        layers,
                       std::shared_ptr<Material> matprop,
                       double                    loc1Min,
                       double                    loc1Max,
                       double                    loc2Min,
                       double                    loc2Max,
                       const std::string&        volumeName = "UndefinedVolume",
                       BinningType               btype      = arbitrary) const;
                       
  /// Create a gap volume from dimensions and
  ///
  /// @param matprop dense material properties for this TrackingVolume
  /// @param loc1Min, loc1Max, loc2Min, loc2Max : local position in space,
  /// this TrackingVolume is restricted to Translation only
  /// @param materialLayers number of material layers (aequidistant binning)
  /// @param cylinder type of layers
  /// @param volumeName  volume name to be given
  ///
  /// @return shared pointer to a new TrackingVolume                       
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
                          
  /// Create a gap volume from dimensions and
  ///
  /// @param matprop dense material properties for this TrackingVolume
  /// @param loc1Min, loc1Max, loc2Min, loc2Max local position in space,
  /// @param layerPositions custom layer positions
  /// @param cylinder type of layers
  /// @param volumeName  : volume name to be given
  /// @param btype (optional) BinningType - arbitrary(default) or equidistant
  ///
  /// @return shared pointer to a new TrackingVolume
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
  /// R or Z by convention 
  ///
  /// @param volumes the volumes to be contained
  ///
  ///
  /// @return shared pointer to a new TrackingVolume                          
  TrackingVolumePtr
  createContainerTrackingVolume(const TrackingVolumeVector& volumes) const;

  /// Set configuration method
  /// 
  /// @param cvbConfig is the configurtion struct assigned
  void
  setConfiguration(const Config& cvbConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// Set logging instance
  ///
  /// @param logger is the logger isntance to be set
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
  ///
  /// @param layers the layers for which the dimensions are checked
  /// @param cylBounds the cylinder volume bounds needed for wrapping
  /// @param transform a transformation of the layers, volume
  /// @param rMinClean the smallest radius given by layers
  /// @param rMaxClean the maximal radius given by layers
  /// @param zMinClean the smallest z extend given by layers
  /// @param zMaxClean the maximal z extend given by layers
  /// @param bValue the binning value in which the binning works
  /// @param bType is the type of binning: equidistant, arbitrary
  bool
  estimateAndCheckDimension(const LayerVector&             layers,
                            const CylinderVolumeBounds*&   cylBounds,
                            std::shared_ptr<Transform3D>&  transform,
                            double&                        rMinClean,
                            double&                        rMaxClean,
                            double&                        zMinClean,
                            double&                        zMaxClean,
                            BinningValue&                  bValue,
                            BinningType                    bType = arbitrary) const;

  /// Private method - interglue all volumes contained by a TrackingVolume
  /// and set the outside glue volumes in the descriptor
  /// 
  /// @param tVolume the tracking volume that is glued together
  /// @param rBinned a boolean indicating if it is binned in r
  /// @param rMin the minimum radius of the volume
  /// @param rGlueMin the minimum glue radius (@todo check and document)
  /// @param rMax the maximim radius of the volume
  /// @param zMin the minimum z extend of the volume
  /// @param zMax the maximum z extend of the volume
  bool
  interGlueTrackingVolume(TrackingVolumePtr tVolume,
                          bool              rBinned,
                          double            rMin,
                          double            rGlueMin,
                          double            rMax,
                          double            zMin,
                          double            zMax) const;

  /// Private method - glue volume to the other 
  ///
  /// @param volumeOne is the first volume in the glue process
  /// @param faceOne is the first boundary face of the glue process
  /// @param volumeTwo is the second volume in the glue process
  /// @param faceTwo is the second boundary face of the glue process
  /// @param rMin the minimum radius of the volume
  /// @param rGlueMin the minimum glue radius (@todo check and document)
  /// @param rMax the maximim radius of the volume
  /// @param zMin the minimum z extend of the volume
  /// @param zMax the maximum z extend of the volume
  void
  glueTrackingVolumes(TrackingVolumePtr   volumeOne,
                      BoundarySurfaceFace faceOne,
                      TrackingVolumePtr   volumeTwo,
                      BoundarySurfaceFace faceTwo,
                      double              rMin,
                      double              rGlueMin,
                      double              rMax,
                      double              zMin,
                      double              zMax) const;

  /// Private method - helper method not to duplicate code
  ///
  /// @param tVolume is the volume to which faces are added
  /// @param bsf is the boundary surface to which faces are added
  /// @param vols are the voluems which are added                                                             
  void
  addFaceVolumes(TrackingVolumePtr      tVolume,
                 BoundarySurfaceFace    bsf,
                 TrackingVolumeVector&  vols) const;

  /// Private method - helper method to save some code
  ///
  /// @param z is the z position of the layer (@todo use Transform)
  /// @param r is the radius of the layer
  /// @param halflength is the half lengthz in z of the cylinder
  /// @param thickness is the thickness of the cylinder
  /// @param binsPhi are the bins for the material in phi
  /// @param binsZ are the bins for the material in z
  ///
  /// @return shared pointer to newly created cylinder layer
  LayerPtr
  createCylinderLayer(double z,
                      double r,
                      double halflength,
                      double thickness,
                      int    binsPhi,
                      int    binsZ) const;

  /// Private method - helper method to save some code
  ///
  /// @param z is the z position of the layer (@todo use Transform)
  /// @param rMin is the minimum radius of the layer
  /// @param rMax is the maximal radius of the layer
  /// @param thickness is the thickness of the cylinder
  /// @param binsPhi are the bins for the material in phi
  /// @param binsR are the bins for the material in R
  ///
  /// @return shared pointer to newly created cylinder layer
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
