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
  /// - create an array in a cylinder, binned in phi, z
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
                         Acts::BinningType                  bType,
                         double                             R,
                         double                             minPhi  = 10e-10,
                         double                             maxPhi  = 10e-10,
                         double                             halfZ   = 10e-10,
                         size_t                             binsPhi = 0,
                         size_t                             binsZ   = 0,
                         std::shared_ptr<Transform3D>       transform
                         = nullptr) const final;

  /// SurfaceArrayCreator interface method
  /// - create an array on a disc, binned in r, phi
  /// @param surfaces are the sensitive surfaces to be
  /// @param rMin is the minimimal radius of the disc
  /// @param rMax is the maximal radius of the disc
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param transform is the (optional) additional transform applied
  /// @return a unique pointer a new SurfaceArray
  std::unique_ptr<SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                     Acts::BinningType                  bType,
                     double                             rMin    = 10e-10,
                     double                             rMax    = 10e-10,
                     double                             minPhi  = 10e-10,
                     double                             maxPhi  = 10e-10,
                     size_t                             binsR   = 0,
                     size_t                             binsPhi = 0,
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

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger)
  {
    m_logger = std::move(logger);
  }

  /// SurfaceArrayCreator internal method
  /// - create a BinUtility
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
  std::unique_ptr<Acts::BinUtility>
  createBinUtility(const std::vector<const Acts::Surface*>& surfaces,
                   Acts::BinningValue                       bValue,
                   Acts::BinningType                        bType,
                   size_t                                   bins = 0,
                   double                                   min  = 10e-10,
                   double                                   max  = 10e-10,
                   std::shared_ptr<Acts::Transform3D>       transform
                   = nullptr) const;

private:
  /// Private access to logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

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
