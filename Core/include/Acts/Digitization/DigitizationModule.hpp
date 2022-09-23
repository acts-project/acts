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

#include <memory>

namespace Acts {

class Surface;
using SurfacePtr = std::shared_ptr<const Surface>;
using SurfacePtrVector = std::vector<SurfacePtr>;

/// @class DigitizationModule
///
/// Class that holds the surfaces for a planar digitization detector module.
///
/// It needs a descriptor to design different pixel/strixels/strip setups
/// (with a segmentation class) in order to define the readout segmentation
///
/// The digitizaiton is done in the local frame of the surface.
///
/// The lorentz angle is assumed to be only in x-direction and constant for the
/// module, it is measured from the local z-direction towards the local
/// x-direction.
///
/// The readout direction defines the charge drift either:
/// a) towards the surface at -halfThickness if readout is defined at -1
/// b) towards the surface at +halfThickness if readout is defined at +1
///
/// Conventions:
///   - 3D positions are within the 3D frame of the module
///   - 2D positions are corrected to parameter surface  at the center of the
///   module (and not the readout surface)
///
/// The lorenzShift is the correction from the readout surface to the parameter
/// surface
///
class DigitizationModule {
 public:
  /// Constructor from a Segmentation descriptor
  ///
  /// @param moduleSegmentation is the segmentation descriptions
  /// @param halfThickness is the half thickness of the module
  /// @param readoutDirection is the readout drift direction
  /// @param lorentzAngle is the lorentz drift angle
  /// @param energyThreshold Optional energy threshold for digitization
  /// @param analogue Run analogue digitization (defaults to false)
  DigitizationModule(std::shared_ptr<const Segmentation> moduleSegmentation,
                     double halfThickness, int readoutDirection,
                     double lorentzAngle, double energyThreshold = 0.,
                     bool analogue = false);

  /// Virtual Destructor
  virtual ~DigitizationModule() = default;

  /// Return the internal test segmentation surfaces to test between entry
  /// and exit given by their cell id's - the boundaries are not given
  ///
  /// @param entryCids are the entry digitisation cell ids
  /// @param exitCids are the exit digitisation cell ids
  ///
  /// @return object is a vector of shared surfaces
  const SurfacePtrVector segmentationSurfaces(
      const DigitizationCell& entryCids,
      const DigitizationCell& exitCids) const;

  /// Get the digitization cell from a position
  /// @param position The position to query
  /// @return
  const DigitizationCell cell(const Vector2& position) const;

  /// Return the module thickness
  double halfThickness() const;

  /// Return the readout direction
  int readoutDirection() const;

  /// Return the lorentz Angle
  double lorentzAngle() const;

  /// Return the energy threshold per cell of the module
  double energyThreshold() const;

  /// Indicates if the readout of the module is analogue, default is digital
  bool analogue() const;

  /// return the segmenation
  const Segmentation& segmentation() const;

  /// Return the test surfaces between these points
  ///
  /// @param start is the start position of the step
  /// @param end is the end position of the step
  ///
  /// @return stepSurfaces are the surfaces to test
  SurfacePtrVector stepSurfaces(const Vector3& start, const Vector3& end) const;

  /// Fill the associated digitization cell from this start and end position,
  /// correct for lorentz effect if needed
  ///
  /// @param start is the start position of the step
  /// @param end is the end position of the step
  /// @return the digitization step
  DigitizationStep digitizationStep(const Vector3& start,
                                    const Vector3& end) const;

  /// Return the bounding surfaces inlcuding top and bottom
  const SurfacePtrVector& boundarySurfaces() const;

  /// Return all surfaces in X - excluding the boundaries
  const SurfacePtrVector& segmentationSurfacesX() const;

  /// Return all surfaces in Y - excluding the boundaries
  const SurfacePtrVector& segmentationSurfacesY() const;

 private:
  /// half thickness of the module
  double m_halfThickness;
  /// readout is along (+1) / (-1) wrt local z axis
  int m_readoutDirection;
  /// the lorentz angle
  double m_lorentzAngle;
  /// and the tangent of it
  double m_tanLorentzAngle;
  /// energy threshold per cell
  double m_energyThreshold;
  /// flag indicating if module is read out analogue
  bool m_analogue;
  /// segmentation descriptor
  std::shared_ptr<const Segmentation> m_segmentation;
  /// boundary surfaces z, x, y
  SurfacePtrVector m_boundarySurfaces;
  /// segmentation surfaces in X - without boundaries
  SurfacePtrVector m_segmentationSurfacesX;
  /// segmentation surfaces in Y - without boundaries
  SurfacePtrVector m_segmentationSurfacesY;
};

inline double DigitizationModule::halfThickness() const {
  return m_halfThickness;
}

inline int DigitizationModule::readoutDirection() const {
  return m_readoutDirection;
}

inline double DigitizationModule::lorentzAngle() const {
  return m_lorentzAngle;
}

inline double DigitizationModule::energyThreshold() const {
  return m_energyThreshold;
}

inline bool DigitizationModule::analogue() const {
  return m_analogue;
}

inline const Segmentation& DigitizationModule::segmentation() const {
  return (*(m_segmentation.get()));
}

inline const SurfacePtrVector& DigitizationModule::boundarySurfaces() const {
  return m_boundarySurfaces;
}

inline const SurfacePtrVector& DigitizationModule::segmentationSurfacesX()
    const {
  return m_segmentationSurfacesX;
}

inline const SurfacePtrVector& DigitizationModule::segmentationSurfacesY()
    const {
  return m_segmentationSurfacesY;
}

inline DigitizationStep DigitizationModule::digitizationStep(
    const Vector3& start, const Vector3& end) const {
  return m_segmentation->digitizationStep(start, end, m_halfThickness,
                                          m_readoutDirection, m_lorentzAngle);
}
}  // namespace Acts
