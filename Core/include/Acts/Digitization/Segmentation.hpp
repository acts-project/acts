// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Digitization/DigitizationCell.hpp"

#include <memory>
#include <vector>

namespace Acts {

class SurfaceBounds;
class Surface;
class BinUtility;
using SurfacePtr = std::shared_ptr<const Surface>;
using SurfacePtrVector = std::vector<SurfacePtr>;

/// @brief Segmentation Base class
///
/// This helper class allows to define an arbitrary readout
/// segmentation for the geoemtric digitization, provided a shape of the module,
/// it creates the segmentation surfaces and hands them to the digitization
/// module
///
/// Since the segmentation description might be identical for many elements
/// while the lorentz angle may change, lorentzAngle and readoutDirection
/// are provided and th the segmenation class only creates the surfaces for the
/// module,
/// but hosts the binning information.
///
class Segmentation {
 public:
  /// Virtual Destructor
  virtual ~Segmentation() = default;

  /// Create the segmentation surfaces in X
  ///
  /// This method is only used if the full 3D digitization is done
  ///
  /// @param boundarySurfaces vector to be filled
  /// @param segmentationSurfacesX are the segmetation boundaries in X
  /// @param segmentationSurfacesY are the segmetation boundaries in Y
  /// @param halfThickness is the half thickness in z of the module
  /// @param readoutDirection is the direction w.r.t normal vector
  /// where the readout is given : -1, 0, 1 possible
  /// @param lorentzAngle is the lorentz angle measured from the local z
  /// towards x axis
  virtual void createSegmentationSurfaces(
      SurfacePtrVector& boundarySurfaces,
      SurfacePtrVector& segmentationSurfacesX,
      SurfacePtrVector& segmentationSurfacesY, double halfThickness,
      int readoutDirection, double lorentzAngle) const = 0;

  /// Get the digitization cell from a 3D position
  /// - ignores the shift, i.e. assumenes in to be in cell frame
  ///
  /// @param position is the position for which the cell is requested
  ///
  /// @return is a cell with cell ids
  virtual DigitizationCell cell(const Vector3& position) const = 0;

  /// Get the digitization cell from a 2D position
  /// - ignores the shift, i.e. assumenes in to be in cell frame
  ///
  /// @param position is the position for which the cell is requested
  ///
  /// @return is a cell with cell ids
  virtual DigitizationCell cell(const Vector2& position) const = 0;

  /// Calculate the cell Position from the Id
  ///
  /// @param dCell the digitization cell
  ///
  /// @return the center position of the associated cell
  virtual Vector2 cellPosition(const DigitizationCell& dCell) const = 0;

  /// Fill the associated digitization cell from the start and end position in
  /// 3D correct for lorentz effect if needed
  ///
  /// @param start is the start position of the step
  /// @param end is the end position of the step
  /// @param halfThickness is the half thickness in z
  /// @param readoutDirection is the readout direction with respect to local z
  /// @param lorentzAngle is the lorentz angle measured from local z towards x
  ///
  /// @return is a fully calculated digitzation step
  virtual DigitizationStep digitizationStep(const Vector3& start,
                                            const Vector3& end,
                                            double halfThickness,
                                            int readoutDirection,
                                            double lorentzAngle) const = 0;

  /// return the surface bounds by reference
  virtual const SurfaceBounds& moduleBounds() const = 0;

  /// return the bin utility that defines the readout
  virtual const BinUtility& binUtility() const = 0;
};
}  // namespace Acts
