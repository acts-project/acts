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

#include "ACTS/Tools/ISurfaceArrayCreator.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class Surface;
class BinUtility;

typedef std::vector<const Surface*> SurfaceVector;
typedef std::vector<SurfaceVector>  SurfaceMatrix;
typedef std::vector<SurfaceMatrix>  SurfaceGrid;

typedef std::vector<Vector3D> V3Vector;
typedef std::vector<V3Vector> V3Matrix;

/// @class SurfaceArrayCreator
///
/// It is designed create sub surface arrays to be ordered on Surfaces
///
class SurfaceArrayCreator : virtual public ISurfaceArrayCreator
{
public:
  /// Constructor
  /// @param logger logging instance
  SurfaceArrayCreator(std::unique_ptr<Logger> logger
                      = getDefaultLogger("SurfaceArrayCreator", Logging::INFO))
    : m_logger(std::move(logger))
  {
  }

  /// Destructor
  virtual ~SurfaceArrayCreator() = default;

  /// SurfaceArrayCreator interface method
  /// - create an array in a cylinder, binned in phi, z when extremas and
  /// bin numbers are known
  /// @param surfaces are the sensitive surfaces to be ordered on the cylinder
  /// @param R is the radius of the cylinder
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param halfZ is the half length in z of the cylinder
  /// @param binsPhi is the number of bins in phi for the surfaces
  /// @param binsX is the number of bin in Z for the surfaces
  /// @param transform is the (optional) additional transform applied
  /// @return a unique pointer to a new SurfaceArray
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
  ///
  /// - create an array in a cylinder, binned in phi, z when extremas and bin
  /// numbers are unknown - this method goes through the surfaces and finds out
  /// the needed information
  /// @param surfaces are the sensitive surfaces to be ordered on the cylinder
  /// @param bTypePhi the binning type in phi direction (equidistant/aribtrary)
  /// @param bTypeZ the binning type in z direction (equidistant/aribtrary)
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<Acts::SurfaceArray>
  surfaceArrayOnCylinder(const std::vector<const Acts::Surface*>& surfaces,
                         Acts::BinningType bTypePhi = equidistant,
                         Acts::BinningType bTypeZ   = equidistant,
                         std::shared_ptr<Acts::Transform3D> transform
                         = nullptr) const final;

  /// SurfaceArrayCreator interface method
  /// - create an array on a disc, binned in r, phi when extremas and
  /// bin numbers are known
  /// @param surfaces are the sensitive surfaces to be
  /// @param rMin is the minimimal radius of the disc
  /// @param rMax is the maximal radius of the disc
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param transform is the (optional) additional transform applied
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                     double                             rMin,
                     double                             rMax,
                     double                             minPhi,
                     double                             maxPhi,
                     size_t                             binsR,
                     size_t                             binsPhi,
                     std::shared_ptr<Transform3D>       transform
                     = nullptr) const final;

  /// SurfaceArrayCreator interface method
  ///
  /// - create an array in a cylinder, binned in phi, r when extremas and bin
  /// numbers are unknown - this method goes through the surfaces and finds out
  /// the needed information
  /// @param surfaces are the sensitive surfaces to be ordered on the cylinder
  /// @param bTypeR the binning type in r direction (equidistant/aribtrary)
  /// @param bTypePhi the binning type in phi direction (equidistant/aribtrary)
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<Acts::SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Acts::Surface*>& surfaces,
                     Acts::BinningType                        bTypeR,
                     Acts::BinningType                        bTypePhi,
                     std::shared_ptr<Acts::Transform3D> transform) const final;

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

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger)
  {
    m_logger = std::move(logger);
  }

private:
  /// Private access to logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }
  /// SurfaceArrayCreator internal method
  /// - create an arbitrary BinUtility when the extremas and the bin numer are
  /// not known
  /// - it will loop through the surfaces and find out the needed information
  /// currently implemented for phi, r and z bining
  /// @TODD implement for x,y binning
  /// @param surfaces are the sensitive surfaces to be
  /// @param bValue the BinningValue in which direction should be binned
  /// (currently possible: binPhi, binR, binZ)
  /// @param bType the BinningType (equidistant or arbitrary binning)
  /// @return a unique pointer a one dimensional BinUtility
  Acts::BinUtility
  createArbitraryBinUtility(const std::vector<const Acts::Surface*>& surfaces,
                            Acts::BinningValue                       bValue,
                            std::shared_ptr<Acts::Transform3D>       transform
                            = nullptr) const;
  /// SurfaceArrayCreator internal method
  /// - create an equidistant BinUtility when the extremas and the bin numer are
  /// - it will loop through the surfaces and find out the needed information
  /// currently implemented for phi, r and z bining
  /// @TODD implement for x,y binning
  /// @param surfaces are the sensitive surfaces to be
  /// @param bValue the BinningValue in which direction should be binned
  /// (currently possible: binPhi, binR, binZ)
  /// @param bType the BinningType (equidistant or arbitrary binning)
  /// @return a unique pointer a one dimensional BinUtility
  Acts::BinUtility
  createEquidistantBinUtility(const std::vector<const Acts::Surface*>& surfaces,
                              Acts::BinningValue                       bValue,
                              std::shared_ptr<Acts::Transform3D>       transform
                              = nullptr) const;
  /// SurfaceArrayCreator internal method
  /// - create an equidistant BinUtility with all parameters given
  /// - if parameters are known this function is preferred
  /// currently implemented for phi, r and z bining
  /// @TODD implement for x,y binning
  /// @param surfaces are the sensitive surfaces to be
  /// @param bValue the BinningValue in which direction should be binned
  /// (currently possible: binPhi, binR, binZ)
  /// @param bType the BinningType (equidistant or arbitrary binning)
  /// @param bins the number of bins
  /// @param min is the minimal position of the surfaces in the binning
  /// direction
  /// @param max is the maximal position of the surfaces in the binning
  /// direction
  /// @return a unique pointer a one dimensional BinUtility
  Acts::BinUtility
  createBinUtility(const std::vector<const Acts::Surface*>& surfaces,
                   Acts::BinningValue                       bValue,
                   Acts::BinningType                        bType,
                   size_t                                   bins,
                   double                                   min,
                   double                                   max,
                   std::shared_ptr<Acts::Transform3D>       transform
                   = nullptr) const;

  /// logging instance
  std::unique_ptr<Logger> m_logger;

  /// Private helper method to complete the binning
  ///
  /// @TODO implement closest neighbour search
  ///
  ///  given a grid point o
  ///    |  0  |  1 |  2  |  3 |  4  |
  ///    ------------------------------
  ///  0 |  x  |    |     |    |  x  |
  ///  1 |     |    |  o  |    |     |
  ///  2 |  x  |    |     |    |  x  |
  ///
  /// This is being called when you chose to use more bins thans surfaces
  /// I.e. to put a finer granularity binning onto your surface
  /// Neighbour bins are then filled to contain pointers as well
  /// @param surfaces are the contained surfaces of the layer
  /// @param binUtility is the unitility class that describes the binning
  /// sVector is the filled vector of Surface and binning position
  /// binSystem is the full system of bins
  void
  completeBinning(const BinUtility& binUtility,
                  const V3Matrix&,
                  const SurfaceVector& sVector,
                  SurfaceGrid&         sGrid) const;

  /// Private helper method to complete the binning
  ///
  ///
  /// @TODO write documentation
  ///
  void
  registerNeighbourHood(const SurfaceArray& sArray) const;
  /// Private helper method to create the bin boundaries out of a vector of
  /// float pairs which represent the boundaries of the surfaces, which do not
  /// need to be attached to each other or can possibly be overlapping
  std::vector<float>
  createBinValues(std::vector<std::pair<float, float>> old) const;
};

}  // end of namespace

#endif  // ACTS_TOOLS_SURFACERARRAYCREATOR_H
