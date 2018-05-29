// This file is part of the Acts project.
//
// Copyright (C) 2016-2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// #ifndef ACTS_DIGITIZATION_POLARSEGMENTATION_H
// #define ACTS_DIGITIZATION_POLARSEGMENTATION_H
//
// #include <memory>
// #include "Acts/Digitization/DigitizationCell.hpp"
// #include "Acts/Digitization/Segmentation.hpp"
// #include "Acts/Surfaces/TrapezoidBounds.hpp"
// #include "Acts/Utilities/BinUtility.hpp"
// #include "Acts/Utilities/Definitions.hpp"

// namespace Acts {
//
// /// @class TrapezoidSegmentation class
// ///
// ///  Segementation class for generic pixel, strixels and strip segmentations
// on
// ///  a trapezoidal shape
// ///
// ///  Conventions:
// ///    - 3D positions are within the 3D frame of the module
// ///    - 2D positions are corrected to the readout surface
// ///      (the lorentzShift is not taken into account fro trapezoidal elements
// ///      yet @TODO)
// ///
// class TrapezoidSegmentation : public Segmentation
// {
// public:
//   /// Constructor for all same-size pixels or strips
//   // (in case numCellsY is set to 1)/
//   TrapezoidSegmentation(std::shared_ptr<const TrapezoidBounds>,
//                         size_t numCellsX,
//                         size_t numCellsY = 1);
//
//   /// TODO contructor from BinUtilities for more complex readouts
//
//   /// Virtual Destructor
//   virtual ~TrapezoidSegmentation();
//
//   /// @copydoc Acts::Segmentation::createSegmentationSurfaces
//   void
//   createSegmentationSurfaces(SurfacePtrVector& boundarySurfaces,
//                              SurfacePtrVector& segmentationSurfacesX,
//                              SurfacePtrVector& segmentationSurfacesY,
//                              double            halfThickness,
//                              int               readoutDirection = 1.,
//                              double            lorentzAngle = 0.) const
//                              override;
//
//   /// @copydoc Segmentation::cell
//   const DigitizationCell
//   cell(const Vector3D& position) const override;
//
//   /// @copydoc Segmentation::cell
//   const DigitizationCell
//   cell(const Vector2D& position) const override;
//
//   /// @copydoc Segmentation::cellPosition
//   const Vector2D
//   cellPosition(const DigitizationCell& cId) const override;
//
//   /// @copydoc Segmentation::digitizationStep
//   const DigitizationStep
//   digitizationStep(const Vector3D& start,
//                    const Vector3D& end,
//                    double          halfThickness,
//                    int             readoutDirection = 1,
//                    double          lorentzAngle     = 0.) const override;
//
//   /// return the surface bounds by reference
//   /// specification for TrapezoidBounds
//   const TrapezoidBounds&
//   moduleBounds() const override;
//
//   /// Return the simple binning parameters
//   size_t
//   numCellsX() const;
//
//   /// Return the simple binning parameters
//   size_t
//   numCellsY() const;
//
// private:
//   /// private helper method - templated
//   template <class T>
//   const DigitizationCell
//   cellT(const T& position) const;
//
//   /// Return the local pitch X at Y
//   double
//   PitchX(const Vector2D& localPos) const;
//
//   /// Return the local sinStereo
//   double
//   sinStereoLocal(const Vector2D& localPos) const;
//
//   /// Return the projected x value on the y=0
//   double
//   projectLocX(const Vector2D& localPos) const;
//
//   /// Return the radius correponding to the given module
//   double
//   radius() const;
//
//   std::shared_ptr<const TrapezoidBounds> m_activeBounds;
//   std::unique_ptr<const BinUtility>      m_binUtility;
//   size_t                                 m_binsX;
//   size_t                                 m_binsY;
// };
//
// inline const TrapezoidBounds&
// TrapezoidSegmentation::moduleBounds() const
// {
//   return (*(m_activeBounds.get()));
// }
//
// inline size_t
// TrapezoidSegmentation::numCellsX() const
// {
//   return m_binsX;
// }
//
// inline size_t
// TrapezoidSegmentation::numCellsY() const
// {
//   return m_binsY;
// }
//
// template <class T>
// const DigitizationCell
// TrapezoidSegmentation::cellT(const T& position) const
// {
//   if (m_binsX == 1)
//     return DigitizationCell(0, m_binUtility->bin(position, 0));
//   else if (m_binsY == 1)
//     return DigitizationCell(m_binUtility->bin(position, 0), 0);
//   return DigitizationCell(m_binUtility->bin(position, 0),
//                           m_binUtility->bin(position, 1));
// }
//
// inline const DigitizationCell
// TrapezoidSegmentation::cell(const Vector3D& position) const
// {
//   Vector3D CorrPosition = position;
//   CorrPosition.x() = projectLocX(Vector2D(CorrPosition.x(),
//   CorrPosition.y()));
//   return cellT<Vector3D>(CorrPosition);
// }
//
// inline const DigitizationCell
// TrapezoidSegmentation::cell(const Vector2D& position) const
// {
//   Vector2D corrPosition = position;
//   corrPosition[0] = projectLocX(Vector2D(corrPosition[0], corrPosition[1]));
//   return cellT<Vector2D>(corrPosition);
// }
//
// double
// TrapezoidSegmentation::radius() const
// {
//   return m_activeBounds->halflengthY()
//       / (m_activeBounds->maxHalflengthX() - m_activeBounds->minHalflengthX())
//       * (m_activeBounds->maxHalflengthX() +
//       m_activeBounds->minHalflengthX());
// }
//
// double
// TrapezoidSegmentation::sinStereoLocal(const Vector2D& localPos) const
// {
//   double oneOverRadius = 1. / radius();
//   double x             = localPos.x();
//   double y             = localPos.y();
//   return -x * oneOverRadius
//       / sqrt((1 + y * oneOverRadius) * (1 + y * oneOverRadius)
//              + x * oneOverRadius * x * oneOverRadius);
// }
//
// double
// TrapezoidSegmentation::PitchX(const Vector2D& localPos) const
// {
//   double tanPhi
//       = (m_activeBounds->maxHalflengthX() - m_activeBounds->minHalflengthX())
//       / m_activeBounds->halflengthY();
//   double lengthXatHit = (radius() + localPos[1]) * tanPhi;
//   return lengthXatHit / float(m_binsX);
// }
//
// double
// TrapezoidSegmentation::projectLocX(const Vector2D& localPos) const
// {
//   return -radius() * tan(asin(sinStereoLocal(localPos)));
// }
// }
//
// #endif
//
