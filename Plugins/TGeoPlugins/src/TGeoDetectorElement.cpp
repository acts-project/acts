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
    const Identifier&            identifier,
    TGeoNode*                    tGeoDetElement,
    const TGeoMatrix*            mGlobal,
    const std::string&           axes,
    double                       scalor)
  : Acts::DetectorElementBase()
  , m_detElement(tGeoDetElement)
  , m_identifier(identifier)
  , m_thickness(0.)
{
  // get the placement and orientation in respect to its mother
  const TGeoMatrix* nodeTransform = (m_detElement->GetMatrix());
  const Double_t* rotation        = nullptr;
  const Double_t* translation     = nullptr;
  
  if (mGlobal){
    // the new combined translation
    TGeoHMatrix nTransform = (*nodeTransform)*(*mGlobal);
    std::string nName  = tGeoDetElement->GetName();
    std::string suffix = "_transform";
    nTransform.SetName((nName+suffix).c_str());
    translation    = nTransform.GetTranslation();
    rotation       = nTransform.GetRotationMatrix();
  } else {
    translation    = (nodeTransform->GetTranslation());  
    rotation       = (nodeTransform->GetRotationMatrix());
  }
  // create the translation
  Vector3D colT(scalor*translation[0], scalor*translation[1], scalor*translation[2]);
  Vector3D colX(rotation[0], rotation[3], rotation[6]);
  Vector3D colY(rotation[1], rotation[4], rotation[7]);
  Vector3D colZ(rotation[2], rotation[5], rotation[8]);
  
  // check if it's a box - always true ... 
  TGeoBBox* box = dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
  // check if it's a trapezoid - unfortunately box is the base of everything
  TGeoTrd2* trapezoid = dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
  // 
  if (axes == "xyz"){
    // identical definition
    m_transform = std::make_shared<Transform3D>(colX,colY,colZ,colT);
    if (trapezoid){
      //
    } else {
      // bounds with x/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor*box->GetDX(), scalor*box->GetDY());
      // thickness 
      m_thickness = scalor*box->GetDZ();
      // assign them 
      m_bounds  = rectangleBounds;
      // create the surface 
      m_surface = std::make_shared<const PlaneSurface>(rectangleBounds,
                                                       *this,
                                                       m_identifier);
    
    }
  } else if (axes == "yzx"){
    // next possibility
    m_transform = std::make_shared<Transform3D>(colY,colZ,colX,colT);
    if (trapezoid){
      //
    } else {
      // bounds with y/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor*box->GetDY(), scalor*box->GetDZ());
      // thickness 
      m_thickness = scalor*box->GetDX();
      // assign them 
      m_bounds  = rectangleBounds;
      // create the surface 
      m_surface = std::make_shared<const PlaneSurface>(rectangleBounds,
                                                       *this,
                                                       m_identifier);
    
    }
  } else {
    // default is "zxy"
    // next possibility
    m_transform = std::make_shared<Transform3D>(colZ,colX,colY,colT);
    if (trapezoid){
      //
    } else {
      // bounds with z/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor*box->GetDZ(), scalor*box->GetDX());
      // thickness 
      m_thickness = scalor*box->GetDY();
      // assign them 
      m_bounds  = rectangleBounds;
      // create the surface 
      m_surface = std::make_shared<const PlaneSurface>(rectangleBounds,
                                                       *this,
                                                       m_identifier);
    }
  }


  //// currently only surfaces with rectangular or trapezoidal shape implemented
  //// check if it's a box
  //TGeoBBox* box = dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
  //// check if it's a trapezoid (only when it's not a box)
  //TGeoTrd2* trapezoid = box ? nullptr :
  //  dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
  //
  //
  //  
  //  
  //  
  //// get the material
  //TGeoMaterial* mat = m_detElement->GetVolume()->GetMaterial();
  //
  //  
  //  
  //  
  // Material      moduleMaterial(mat->GetRadLen(),
  //                        mat->GetIntLen(),
  //                        mat->GetA(),
  //                        mat->GetZ(),
  //                        mat->GetDensity());
  //  
  //  
  //  
  //  
  //  
  /**
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
    if (transformToGlobal)
      m_transform = std::make_shared<const Acts::Transform3D>(
          (*transformToGlobal) * (*m_transform)
          * AngleAxis3D(0.5 * M_PI, Vector3D::UnitX()));
    // extract the surface bounds
    auto trapezoidBounds = std::make_shared<const Acts::TrapezoidBounds>(
        scalor*trapezoid->GetDx1(), 
        scalor*trapezoid->GetDx2(), 
        scalor*trapezoid->GetDz());
    m_bounds  = trapezoidBounds;
    m_surface = std::make_shared<const Acts::PlaneSurface>(trapezoidBounds, 
                                                           *this, 
                                                           m_identifier);
    // ignore module material for the moment @TODO handle module material
    
        MaterialProperties moduleMaterialProperties(
            moduleMaterial,
            0.5 * scalor * (trapezoid->GetDy1() + trapezoid->GetDy1()));
        m_surface->setAssociatedMaterial(std::shared_ptr<const SurfaceMaterial>(
            new HomogeneousSurfaceMaterial(moduleMaterialProperties)));
    // set the thickness
    m_thickness = scalor*2*trapezoid->GetDy1();
  } else {
    m_transform = std::make_shared<Acts::Transform3D>(
        Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
        Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
        Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
        Acts::Vector3D(
            translation[0] * cm, translation[1] * cm, translation[2] * cm));
    // now calculate the global transformation
    if (transformToGlobal) {
      m_transform = std::make_shared<const Acts::Transform3D>((*transformToGlobal)
                                                               *(*m_transform));
    }
    // extract the surface bounds
    auto rectangleBounds = std::make_shared<const Acts::RectangleBounds>(scalor*box->GetDX(), scalor*box->GetDY());
    m_bounds  = rectangleBounds;
    m_surface = std::make_shared<const Acts::PlaneSurface>(rectangleBounds,
                                                           *this,
                                                          m_identifier);
    // set the thickness
    m_thickness = scalor*2*box->GetDZ();
  }

  **/
}

Acts::TGeoDetectorElement::~TGeoDetectorElement()
{
}
