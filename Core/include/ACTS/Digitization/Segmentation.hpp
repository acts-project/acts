// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DIGITIZATION_SEGMENTATION_H
#define ACTS_DIGITIZATION_SEGMENTATION_H 1

#include <memory>
#include <vector>
#include "ACTS/Digitization/DigitizationCell.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

class SurfaceBounds;
class Surface;
typedef std::shared_ptr<const Surface> SurfacePtr;
typedef std::vector<SurfacePtr>        SurfacePtrVector;

/// @class Segmentation Base class
///
/// This helper class allows to define an arbitrary readout
/// segmentation for the geoemtric digitization, provided a shape of the module,
/// it creates the segmentation surfaces and hands them to the digitization
/// module
///
/// Since the segmentation description might be identical for many elements
/// while
/// the lorentz angle may change, lorentzAngle and readoutDirection are provided
/// and th the segmenation class only creates the surfaces for the module, but
/// hosts the binning information.
///
class Segmentation
{
public:
  /// Virtual Destructor
  virtual ~Segmentation() {}
  
  /// Create the segmentation surfaces in X
  ///
  /// @param boundarySurfaces vector to be filled
  /// @param segmentationSurfacesX are the segmetation boundaries in X
  /// @param segmentationSurfacesY are the segmetation boundaries in Y
  /// @param halfThickness is the half thickness in z of the module
  virtual void
  createSegmenationSurfaces(SurfacePtrVector& boundarySurfaces,
                            SurfacePtrVector& segmentationSurfacesX,
                            SurfacePtrVector& segmentationSurfacesY,
                            double            halfThickness,
                            int               readoutDirection,
                            double            lorentzAngle) const = 0;

  /// Get the digitization cell fropm a 3D position
  /// - ignores the shift, i.e. assumenes in to be in cell frame
  ///
  /// @param position is the position for which the cell is requested
  ///                          
  /// @return is a cell with cell ids
  virtual const DigitizationCell
  cell(const Vector3D& position) const = 0;

  /// Get the digitization cell fropm a 2D position
  /// - ignores the shift, i.e. assumenes in to be in cell frame
  ///
  /// @param position is the position for which the cell is requested
  ///
  /// @return is a cell with cell ids
  virtual const DigitizationCell
  cell(const Vector2D& position) const = 0;

  /// Calculate the cell Position from the Id
  ///
  /// @param dCell is the digitization cell 
  ///
  /// @return the center position of the associated cell
  virtual const Vector2D
  cellPosition(const DigitizationCell& dCell) const = 0;

  /// Fill the associated digitsation cell from the start and end position in 3D
  /// correct for lorentz effect if needed
  ///
  /// @param start is the start position of the step
  /// @param end is the end position of the step
  /// @param halfThickness is the half thickness in z
  /// @param readoutDirection is the readout direction with respect to local z
  /// @param lorentzAngle is the lorentz angle measured from local z towards x
  ///
  /// @return is a fully calculated digitzation step
  virtual const DigitizationStep
  digitizationStep(const Vector3D& start,
                   const Vector3D& end,
                   double          halfThickness,
                   int             readoutDirection,
                   double          lorentzAngle) const = 0;

  // return the surface bounds by reference
  virtual const SurfaceBounds&
  moduleBounds() const = 0;
};
}  // end of namespace Acts

#endif
