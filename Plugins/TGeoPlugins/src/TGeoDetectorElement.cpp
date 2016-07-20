// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
//#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
//#include "ACTS/Material/Material.hpp"
//#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "TGeoBBox.h"
#include "TGeoTrd2.h"
#include <iostream>

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier&                        identifier,
    TGeoNode*                                tGeoDetElement,
    std::shared_ptr<const Acts::Transform3D> motherTransform)
  : Acts::DetectorElementBase()
  , m_detElement(tGeoDetElement)
  , m_identifier(identifier)
  , m_thickness(0.)
{
  // get the placement and orientation in respect to its mother
  const Double_t* rotation = (m_detElement->GetMatrix()->GetRotationMatrix());
  const Double_t* translation = (m_detElement->GetMatrix()->GetTranslation());
  // currently only surfaces with rectangular or trapezoidal shape implemented
  TGeoBBox* box
      = dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
  TGeoTrd2* trapezoid
      = dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
//  TGeoMaterial* mat = m_detElement->GetVolume()->GetMaterial();
/*  Material      moduleMaterial(mat->GetRadLen(),
                          mat->GetIntLen(),
                          mat->GetA(),
                          mat->GetZ(),
                          mat->GetDensity());*/
  if (trapezoid) {
    // if the shape is TGeoTrd2 y and z axes needs to be exchanged, since in
    // TGei the description is different
    m_transform = std::make_shared<Acts::Transform3D>(
        Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
        Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
        Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
        Acts::Vector3D(
            translation[0] * cm, translation[1] * cm, translation[2] * cm));
    // now calculate the global transformation
    if (motherTransform)
      m_transform = std::make_shared<const Acts::Transform3D>(
          (*motherTransform) * (*m_transform)
          * AngleAxis3D(0.5 * M_PI, Vector3D::UnitX()));
    // extract the surface bounds
    auto trapezoidBounds = std::make_shared<const Acts::TrapezoidBounds>(
        trapezoid->GetDx1() * cm,
        trapezoid->GetDx2() * cm,
        trapezoid->GetDz() * cm);
    m_bounds  = trapezoidBounds;
    m_surface = std::make_shared<const Acts::PlaneSurface>(trapezoidBounds, 
                                                           *this, 
                                                           m_identifier);
      // ignore module material for the moment @TODO handle module material
/*
    MaterialProperties moduleMaterialProperties(
        moduleMaterial,
        0.5 * (trapezoid->GetDy1() * cm + trapezoid->GetDy1() * cm));
    m_surface->setAssociatedMaterial(std::shared_ptr<const SurfaceMaterial>(
        new HomogeneousSurfaceMaterial(moduleMaterialProperties)));*/
  } else {
    m_transform = std::make_shared<Acts::Transform3D>(
        Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
        Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
        Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
        Acts::Vector3D(
            translation[0] * cm, translation[1] * cm, translation[2] * cm));
    // now calculate the global transformation
    if (motherTransform) {
      m_transform = std::make_shared<const Acts::Transform3D>((*motherTransform)
                                                               *(*m_transform));
    }
    // extract the surface bounds
    auto rectangleBounds = std::make_shared<const Acts::RectangleBounds>(
        box->GetDX() * cm, box->GetDY() * cm);
    m_bounds  = rectangleBounds;
    m_surface = std::make_shared<const Acts::PlaneSurface>(rectangleBounds,
                                                           *this,
                                                          m_identifier);
      // ignore module material for the moment @TODO handle module material
/*    MaterialProperties moduleMaterialProperties(moduleMaterial,
                                                box->GetDZ() * cm);
    m_surface->setAssociatedMaterial(std::shared_ptr<const SurfaceMaterial>(
        new HomogeneousSurfaceMaterial(moduleMaterialProperties)));*/
  }
}

Acts::TGeoDetectorElement::~TGeoDetectorElement()
{
}
