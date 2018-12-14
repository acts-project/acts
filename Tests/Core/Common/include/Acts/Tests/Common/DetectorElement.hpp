// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  ///
  /// @brief Stub implementation of the detector element
  ///
  class DetectorElement : public DetectorElementBase
  {
  public:
    /// @brief Constructor
    ///
    /// @param [in] transformation Transformation of the detector element
    /// @param [in] pBounds Planar boundaries of the plane surface
    /// @param [in] thickness Thickness of the detector element
    /// @param [in] moduleMaterial The material associated to teh module
    DetectorElement(std::shared_ptr<const Transform3D>     transformation,
                    std::shared_ptr<const PlanarBounds>    pBounds,
                    double                                 thickness,
                    std::shared_ptr<const SurfaceMaterial> moduleMaterial
                    = nullptr)
      : DetectorElementBase()
      , m_transformation(transformation)
      , m_surface(Surface::makeShared<PlaneSurface>(pBounds, *this))
      , m_thickness(thickness)
    {
      // Assign material if you have
      if (moduleMaterial) {
        m_surface->setAssociatedMaterial(moduleMaterial);
      }
    }

    /// @brief Getter of the transformation
    virtual const Transform3D&
    transform() const
    {
      return *m_transformation;
    }

    /// @brief Getter of the surface
    virtual const Surface&
    surface() const
    {
      return *(m_surface);
    }

    /// @brief Getter of the thickness
    virtual double
    thickness() const
    {
      return m_thickness;
    }

    // Pointer to the transformation
    std::shared_ptr<const Transform3D> m_transformation = nullptr;

    // Surface related to the detector element
    std::shared_ptr<Surface> m_surface = nullptr;

    // Thickness of the detector element
    double m_thickness = 0.;
  };

}  // namespace Test
}  // namespace Acts
