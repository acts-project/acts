// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericDetectorElement.cxx, ACTS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#include "ACTS/Examples/GenericDetectorElement.hpp"
#include "ACTS/Surfaces/DiscBounds.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"

Acts::GenericDetectorElement::GenericDetectorElement(
    const Identifier                             identifier,
    std::shared_ptr<Acts::Transform3D>           transform,
    std::shared_ptr<const Acts::PlanarBounds>    pBounds,
    double                                       thickness,
    std::shared_ptr<const Acts::SurfaceMaterial> material)
  : DetectorElementBase()
  , m_elementIdentifier(std::move(identifier))
  , m_elementTransform(std::move(transform))
  , m_elementSurface(new PlaneSurface(pBounds, *this))
  , m_elementThickness(thickness)
  , m_elementSurfaces({m_elementSurface})
  , m_elementPlanarBounds(std::move(pBounds))
  , m_elementDiscBounds(nullptr)
{
  m_elementSurface->setAssociatedMaterial(material);
}

Acts::GenericDetectorElement::GenericDetectorElement(
    const Identifier                             identifier,
    std::shared_ptr<Acts::Transform3D>           transform,
    std::shared_ptr<const Acts::DiscBounds>      dBounds,
    double                                       thickness,
    std::shared_ptr<const Acts::SurfaceMaterial> material)
  : DetectorElementBase()
  , m_elementIdentifier(std::move(identifier))
  , m_elementTransform(std::move(transform))
  , m_elementSurface(new DiscSurface(dBounds, *this))
  , m_elementThickness(thickness)
  , m_elementSurfaces({m_elementSurface})
  , m_elementPlanarBounds(nullptr)
  , m_elementDiscBounds(std::move(dBounds))
{
  m_elementSurface->setAssociatedMaterial(material);
}

Acts::GenericDetectorElement::~GenericDetectorElement()
{
}
