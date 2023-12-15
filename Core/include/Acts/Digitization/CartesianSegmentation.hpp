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
#include "Acts/Digitization/Segmentation.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {

/// @brief Segmentation Base class
///
/// Segmentation class for generic pixel, strixels and strip segmentations
/// in a cartesian frame, this uses a cartesian X/Y local surface definition
///
/// The calculation can be done in full 3D, i.e. the segments of the path
/// through
/// the planar module are calculated in a 3D model - or in 2D, when the entire
/// calculation is done on the projective surface. When the 2D option is used,
/// segmentation surfaces are not created. The 2D calculation is faster and uses
/// less memory, however, effects within the sensor volume can not be easily
/// integrated
///
/// Conventions:
///   - 3D positions are within the 3D frame of the module
///   - 2D positions are corrected to the readout surface
///     they need to be corrected by the lorentzShift
///     for the parameter surface in the center of the surface)
///
class CartesianSegmentation : public Segmentation {
 public:
  /// Constructor for all same-size pixels or strips
  /// (in cas numCellsY is set to 1)
  ///
  /// @param mBounds are the rectangle bounds of the sensitive volume
  /// @param numCellsX is the number of cells in X
  /// @param numCellsY is the number of cells in Y
  CartesianSegmentation(const std::shared_ptr<const PlanarBounds>& mBounds,
                        std::size_t numCellsX, std::size_t numCellsY = 1);

  /// @todo constructor from BinUtilities for more complex readouts
  ///
  /// @param bUtility is the bin Utility,
  //  it will define the RectangleBounds if none are provided
  /// @param mBounds are the rectangle bounds if provided for memory
  /// optimisation
  ///
  /// @note if both RectangleBounds and BinUtility are provided, no check is
  /// done for consistency
  CartesianSegmentation(std::shared_ptr<const BinUtility> bUtility,
                        std::shared_ptr<const PlanarBounds> mBounds = nullptr);

  /// Virtual Destructor
  ~CartesianSegmentation() override;

  /// @copydoc Acts::Segmentation::createSegmentationSurfaces
  ///
  /// Create the segmentation surfaces in X and Y for rectangular shapes
  /// These are needed for a full three dimensional module
  void createSegmentationSurfaces(SurfacePtrVector& boundarySurfaces,
                                  SurfacePtrVector& segmentationSurfacesX,
                                  SurfacePtrVector& segmentationSurfacesY,
                                  double halfThickness,
                                  int readoutDirection = 1.,
                                  double lorentzAngle = 0.) const final;

  /// @copydoc Segmentation::cell
  DigitizationCell cell(const Vector3& position) const final;

  /// @copydoc Segmentation::cell
  DigitizationCell cell(const Vector2& position) const final;

  /// @copydoc Segmentation::cellPosition
  Vector2 cellPosition(const DigitizationCell& dCell) const final;

  /// Fill the associated digitization cell from this start and end position
  /// correct for lorentz effect if needed
  ///
  /// @copydoc Segmentation::digitizationStep
  DigitizationStep digitizationStep(const Vector3& start, const Vector3& end,
                                    double halfThickness,
                                    int readoutDirection = 1,
                                    double lorentzAngle = 0.) const final;

  /// return the surface bounds by reference
  /// specialization for Rectangle Bounds
  const PlanarBounds& moduleBounds() const final;

  /// return the bin utility that defines the
  /// readout segmentation
  const BinUtility& binUtility() const final;

  /// return the pitch sizes as a pair
  std::pair<double, double> pitch() const;

 private:
  template <class T>
  DigitizationCell cellT(const T& position) const;

  std::shared_ptr<const PlanarBounds> m_activeBounds;  /// active area size
  std::shared_ptr<const BinUtility> m_binUtility;      /// bin Utility
};

inline const PlanarBounds& CartesianSegmentation::moduleBounds() const {
  return (*(m_activeBounds.get()));
}

inline const BinUtility& CartesianSegmentation::binUtility() const {
  return (*(m_binUtility.get()));
}

template <class T>
DigitizationCell CartesianSegmentation::cellT(const T& position) const {
  return DigitizationCell(m_binUtility->bin(position, 0),
                          m_binUtility->bin(position, 1));
}

inline DigitizationCell CartesianSegmentation::cell(
    const Vector3& position) const {
  return cellT<Vector3>(position);
}

inline DigitizationCell CartesianSegmentation::cell(
    const Vector2& position) const {
  return cellT<Vector2>(position);
}

inline std::pair<double, double> CartesianSegmentation::pitch() const {
  auto boundingBox = m_activeBounds->boundingBox();
  auto values = boundingBox.values();
  double pitchX = 2. * values[0] / m_binUtility->bins(0);
  double pitchY = 2. * values[1] / m_binUtility->bins(1);
  return std::pair<double, double>(pitchX, pitchY);
}

}  // namespace Acts
