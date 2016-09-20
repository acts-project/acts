// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DIGITIZATION_RECTANGULARSEGMENTATION_H
#define ACTS_DIGITIZATION_RECTANGULARSEGMENTATION_H

#include <memory>
#include "ACTS/Digitization/DigitizationCell.hpp"
#include "ACTS/Digitization/Segmentation.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/// @class Segmentation Base class
///
/// Segementation class for generic pixel, strixels and strip segmentations on a
/// rectangle shape
///
/// Conventions:
///   - 3D positions are within the 3D frame of the module
///   - 2D positions are corrected to the readout surface
///     they need to be corrected by the lorentzShift
///     for the parameter surface in the center of the surface)
///
class RectangularSegmentation : public Segmentation
{
public:
  /// Constructor for all same-size pixels or strips (in cas numCellsY is set to
  /// 1)
  RectangularSegmentation(std::shared_ptr<const RectangleBounds>,
                          size_t numCellsX,
                          size_t numCellsY = 1);

  /// @TODO contructor from BinUtilities for more complex readouts

  /// Virtual Destructor */
  virtual ~RectangularSegmentation();

  /// @copydoc Segmentation::createSegmenationSurfaces
  ///
  /// Create the segmentation surfaces in X and Y for rectangular shapes
  void
  createSegmenationSurfaces(SurfacePtrVector& boundarySurfaces,
                            SurfacePtrVector& segmentationSurfacesX,
                            SurfacePtrVector& segmentationSurfacesY,
                            double            halfThickness,
                            int               readoutDirection = 1.,
                            double            lorentzAngle = 0.) const override;

  /// @copydoc Segmentation::cell
  const DigitizationCell
  cell(const Vector3D& position) const override;

  /// @copydoc Segmentation::cell
  const DigitizationCell
  cell(const Vector2D& position) const override;

  /// @copydoc Segmentation::cellPosition
  const Vector2D
  cellPosition(const DigitizationCell& cId) const override;

  /// Fill the associated digitsation cell from this start and end positio
  /// correct for lorentz effect if needed
  ///
  /// @copydoc Segmentation::digitizationStep
  const DigitizationStep
  digitizationStep(const Vector3D& start,
                   const Vector3D& end,
                   double          halfThickness,
                   int             readoutDirection = 1,
                   double          lorentzAngle     = 0.) const override;

  /// return the surface bounds by reference
  /// specialization for Rectangle Bounds
  const RectangleBounds&
  moduleBounds() const override;

  /// Return the simple binning parameters
  size_t
  numCellsX() const;

  /// Return the simple binning parameters
  size_t
  numCellsY() const;

private:
  template <class T>
  const DigitizationCell
  cellT(const T& position) const;

  std::shared_ptr<const RectangleBounds> m_activeBounds;  /// active area size
  std::unique_ptr<BinUtility>            m_binUtility;    /// bin Utility
  size_t                                 m_binsX;  /// number of bins in X
  size_t                                 m_binsY;  /// number of bins in Y
};

inline const RectangleBounds&
RectangularSegmentation::moduleBounds() const
{
  return (*(m_activeBounds.get()));
}

inline size_t
RectangularSegmentation::numCellsX() const
{
  return m_binsX;
}

inline size_t
RectangularSegmentation::numCellsY() const
{
  return m_binsY;
}

template <class T>
const DigitizationCell
RectangularSegmentation::cellT(const T& position) const
{
  if (m_binsX == 1)
    return DigitizationCell(0, m_binUtility->bin(position, 0));
  else if (m_binsY == 1)
    return DigitizationCell(m_binUtility->bin(position, 0), 0);
  return DigitizationCell(m_binUtility->bin(position, 0),
                          m_binUtility->bin(position, 1));
}

inline const DigitizationCell
RectangularSegmentation::cell(const Vector3D& position) const
{
  return cellT<Vector3D>(position);
}

inline const DigitizationCell
RectangularSegmentation::cell(const Vector2D& position) const
{
  return cellT<Vector2D>(position);
}
}

#endif
