// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_SURFACERARRAYCREATOR_H
#define ACTS_TOOLS_SURFACERARRAYCREATOR_H 1

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Tools/ISurfaceArrayCreator.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class Surface;
class BinUtility;

typedef std::pair<const Surface*, Vector3D>  SurfacePosition;
typedef std::pair<SurfacePosition, Vector3D> SurfacePositionDirection;

/// @class SurfaceArrayCreator
///
/// It is designed create sub surface arrays to be ordered on Surfaces
///
class SurfaceArrayCreator : virtual public ISurfaceArrayCreator
{
public:
  /// @struct Config
  /// Nested configuration class for this surface array creator
  struct Config
  {
    std::shared_ptr<Logger> logger;  ///< logging instance
    
    Config() : logger(getDefaultLogger("SurfaceArrayCreator", Logging::INFO)) {}
  };

  /// Constructor 
  /// @param cfg is the configuration struct for this component
  SurfaceArrayCreator(const Config& cfg) : m_config(cfg) {}
  
  /// Destructor 
  virtual ~SurfaceArrayCreator() = default;

  /// SurfaceArrayCreator interface method 
  /// - create an array in a cylinder, binned in phi, z 
  /// @param surfaces are the sensitive surfaces to be ordered on the cylinder
  /// @param R is the radius of the cylinder
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param halfZ is the half length in z of the cylinder
  /// @param binsPhi is the number of bins in phi for the surfaces
  /// @param binsX is the number of bin in Z for the surfaces   
  /// @param transform is the (optional) additional transform applied
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                         double                             R,
                         double                             minPhi,
                         double                             maxPhi,
                         double                             halfZ,
                         size_t                             binsPhi,
                         size_t                             binsZ,
                         std::shared_ptr<Transform3D>       transform
                         = nullptr) const final;

  /// SurfaceArrayCreator interface method 
  /// - create an array on a disc, binned in r, phi 
  /// @param surfaces are the sensitive surfaces to be 
  /// @param rMin is the minimimal radius of the disc
  /// @param rMax is the maximal radius of the disc
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param rBoundaries are the optional boundaris of the r rings 
  /// @param transform is the (optional) additional transform applied
  /// @return a unique pointer a new SurfaceArray    
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                     double                             rMin,
                     double                             rMax,
                     double                             minPhi,
                     double                             maxPhi,
                     size_t                             binsR,
                     size_t                             binsZ,
                     const std::vector<double>&         rBoundaries = {},
                     std::shared_ptr<Transform3D>       transform
                     = nullptr) const final;

  /// SurfaceArrayCreator interface method 
  /// - create an array on a plane
  /// @param surfaces are the sensitive surfaces to be 
  /// @param halflengthX is the half length in X
  /// @param halflengthY is the half length in Y
  /// @param binsX is the number of bins in X
  /// @param binsY is the number of bins in Y
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnPlane(const std::vector<const Surface*>& surfaces,
                      double                             halflengthX,
                      double                             halflengthY,
                      size_t                             binsX,
                      size_t                             binsY,
                      std::shared_ptr<Transform3D>       transform
                      = nullptr) const final;

  /// set the configuration class
  /// @param cfg is a Configuraiton struct                    
  void
  setConfiguration(const Config& cfg)
  {
    m_config = cfg;
  }

  /// Get configuration method 
  /// @return the current congifuration
  Config
  getConfiguration() const
  {
    return m_config;
  }

private:
  /// Configuration struct 
  Config m_config;
  
  /// Private access to logger 
  const Logger&
  logger() const
  {
    return *m_config.logger;
  }

  /// Private helper method to complete the binning
  /// This is being called when you chose to use more bins thans surfaces
  /// I.e. to put a finer granularity binning onto your surface
  /// Neighbour bins are then filled to contain pointers as well
  /// @param surfaces are the contained surfaces of the layer
  /// @param binUtility is the unitility class that describes the binning
  /// sVector is the filled vector of Surface and binning position 
  /// binSystem is the full system of bins 
  void
  completeBinning(
      const std::vector<const Surface*>&                  surfaces,
      const BinUtility&                                   binUtility,
      std::vector<SurfacePosition>&                       sVector,
      std::vector<std::vector<SurfacePositionDirection>>& binSystem) const;

  /// Register the neighbours on a Grid - needs to be a BinnedArray1D or
  //  BinnedArray2D type binning 
  void
  registerNeighboursGrid(
      const std::vector<std::vector<const Surface*>>& surfaceArrayObjects,
      bool                                            open0,
      bool                                            open1) const;
};

}  // end of namespace

#endif  // ACTS_TOOLS_SURFACERARRAYCREATOR_H
