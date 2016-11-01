// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ISurfaceArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_ISURFACEARRAYCREATOR_H
#define ACTS_TOOLS_ISURFACEARRAYCREATOR_H 1

#include <vector>
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinningType.hpp"

namespace Acts {

/// forward declarations & typedef
class Surface;
typedef BinnedArray<const Surface*> SurfaceArray;

/// @class ISurfaceArrayCreator
///
/// Interface class ISurfaceArrayCreators, it inherits from IAlgTool.
///
class ISurfaceArrayCreator
{
public:
  /// Virtual destructor
  virtual ~ISurfaceArrayCreator() {}
  /// SurfaceArrayCreator interface method
  ///
  /// - create an array in a cylinder, binned in phi, z when extremas and
  /// bin numbers are known
  /// @param surfaces are the sensitive surfaces to be ordered on the cylinder
  /// @param R is the radius of the cylinder
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param halfZ is the half length in z of the cylinder
  /// @param binsPhi is the number of bins in phi for the surfaces
  /// @param binsZ is the number of bin in Z for the surfaces
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  virtual std::unique_ptr<SurfaceArray>
  surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                         double                             R,
                         double                             minPhi,
                         double                             maxPhi,
                         double                             halfZ,
                         size_t                             binsPhi,
                         size_t                             binsZ,
                         std::shared_ptr<Transform3D>       transform
                         = nullptr) const = 0;

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
  virtual std::unique_ptr<Acts::SurfaceArray>
  surfaceArrayOnCylinder(const std::vector<const Acts::Surface*>& surfaces,
                         Acts::BinningType bTypePhi = equidistant,
                         Acts::BinningType bTypeZ   = equidistant,
                         std::shared_ptr<Acts::Transform3D> transform
                         = nullptr) const = 0;

  /// SurfaceArrayCreator interface method
  ///
  /// - create an array on a disc, binned in r, phi when extremas and
  /// bin numbers are known
  /// @param surfaces are the sensitive surfaces to be
  /// @param rMin is the minimimal radius of the disc
  /// @param rMax is the maximal radius of the disc
  /// @param minPhi is the minimal phi position of the surfaces
  /// @param maxPhi is the maximal phi position of the surfaces
  /// @param binsR is the number of bins in r for the surfaces
  /// @param binsPhi is the number of bins in phi for the surfaces
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  virtual std::unique_ptr<SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                     double                             rMin,
                     double                             rMax,
                     double                             minPhi,
                     double                             maxPhi,
                     size_t                             binsR,
                     size_t                             binsPhi,
                     std::shared_ptr<Transform3D>       transform
                     = nullptr) const = 0;

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
  virtual std::unique_ptr<Acts::SurfaceArray>
  surfaceArrayOnDisc(const std::vector<const Acts::Surface*>& surfaces,
                     Acts::BinningType                        bTypeR,
                     Acts::BinningType                        bTypePhi,
                     std::shared_ptr<Acts::Transform3D>       transform
                     = nullptr) const = 0;

  /// SurfaceArrayCreator interface method
  /// - create an array on a plane
  /// @param surfaces are the sensitive surfaces to be
  /// @param halflengthX is the half length in X
  /// @param halflengthY is the half length in Y
  /// @param binsX is the number of bins in X
  /// @param binsY is the number of bins in Y
  /// @param transform is the (optional) additional transform applied
  ///
  /// @return a unique pointer a new SurfaceArray
  virtual std::unique_ptr<SurfaceArray>
  surfaceArrayOnPlane(const std::vector<const Surface*>& surfaces,
                      double                             halflengthX,
                      double                             halflengthY,
                      size_t                             binsX,
                      size_t                             binsY,
                      std::shared_ptr<Transform3D>       transform
                      = nullptr) const = 0;
};

}  // end of namespace

#endif  // ACTS_TOOLS_ISURFACEARRAYCREATOR_H
