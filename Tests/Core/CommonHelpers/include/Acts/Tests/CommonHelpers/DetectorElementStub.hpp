// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetectorElementStub.h, Acts project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Utilities/Definitions.hpp"
// Geometry module
#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

// required because LineSurface is abstract
#include "Acts/Tests/CommonHelpers/LineSurfaceStub.hpp"

namespace Acts {

class PlanarBounds;
class DiscBounds;
class SurfaceMaterial;
class LineBounds;

namespace Test {

  /// @class DetectorElementStub
  ///
  /// This is a lightweight type of detector element,
  /// it simply implements the base class.
  class DetectorElementStub : public DetectorElementBase
  {
  public:
    DetectorElementStub() : DetectorElementBase() {}

    DetectorElementStub(std::shared_ptr<const Transform3D> transform)
      : DetectorElementBase(), m_elementTransform(std::move(transform))
    {
    }

    /// Constructor for single sided detector element
    /// - bound to a Plane Surface
    ///
    /// @param transform is the transform that element the layer in 3D frame
    /// @param pBounds is the planar bounds for the planar detector element
    /// @param thickness is the module thickness
    /// @param material is the (optional) Surface material associated to it
    DetectorElementStub(std::shared_ptr<const Transform3D>     transform,
                        std::shared_ptr<const PlanarBounds>    pBounds,
                        double                                 thickness,
                        std::shared_ptr<const SurfaceMaterial> material
                        = nullptr)
      : DetectorElementBase()
      , m_elementTransform(std::move(transform))
      , m_elementThickness(thickness)
    {
      auto mutableSurface = Surface::makeShared<PlaneSurface>(pBounds, *this);
      mutableSurface->setAssociatedMaterial(material);
      m_elementSurface = mutableSurface;
    }

    /// Constructor for single sided detector element
    /// - bound to a Line Surface
    ///
    /// @param transform is the transform that element the layer in 3D frame
    /// @param dBounds is the line bounds for the line like detector element
    /// @param thickness is the module thickness
    /// @param material is the (optional) Surface material associated to it
    DetectorElementStub(std::shared_ptr<const Transform3D>     transform,
                        std::shared_ptr<const LineBounds>      lBounds,
                        double                                 thickness,
                        std::shared_ptr<const SurfaceMaterial> material
                        = nullptr)
      : DetectorElementBase()
      , m_elementTransform(std::move(transform))
      , m_elementThickness(thickness)
    {
      auto mutableSurface
          = Surface::makeShared<LineSurfaceStub>(lBounds, *this);
      mutableSurface->setAssociatedMaterial(material);
      m_elementSurface = mutableSurface;
    }

    ///  Destructor
    ~DetectorElementStub() override { /*nop */}

    /// Return local to global transform associated with this identifier
    ///
    /// @note this is called from the surface().transform() in the PROXY mode
    const Transform3D&
    transform() const override;

    /// Return surface associated with this detector element
    const Surface&
    surface() const override;

    /// The maximal thickness of the detector element wrt normal axis
    double
    thickness() const override;

  private:
    /// the transform for positioning in 3D space
    std::shared_ptr<const Transform3D> m_elementTransform;
    /// the surface represented by it
    std::shared_ptr<const Surface> m_elementSurface{nullptr};
    /// the element thickness
    double m_elementThickness{0.};
  };

  inline const Transform3D&
  DetectorElementStub::transform() const
  {
    return *m_elementTransform;
  }

  inline const Surface&
  DetectorElementStub::surface() const
  {
    return *m_elementSurface;
  }

  inline double
  DetectorElementStub::thickness() const
  {
    return m_elementThickness;
  }
}
}  // end of ns
