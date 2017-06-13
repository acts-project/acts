// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "TGeoBBox.h"
#include "TGeoTrd2.h"

Acts::TGeoDetectorElement::TGeoDetectorElement(const Identifier& identifier,
                                               TGeoNode*         tGeoDetElement,
                                               const TGeoMatrix* mGlobal,
                                               const std::string& axes,
                                               double             scalor)
  : Acts::DetectorElementBase()
  , m_detElement(tGeoDetElement)
  , m_identifier(identifier)
  , m_thickness(0.)
{
  // get the placement and orientation in respect to its mother
  const TGeoMatrix* nodeTransform = (m_detElement->GetMatrix());
  const Double_t*   rotation      = nullptr;
  const Double_t*   translation   = nullptr;

  if (mGlobal) {
    // the new combined translation
    TGeoHMatrix nTransform = (*mGlobal) * (*nodeTransform);
    std::string nName      = tGeoDetElement->GetName();
    std::string suffix     = "_transform";
    nTransform.SetName((nName + suffix).c_str());
    translation = nTransform.GetTranslation();
    rotation    = nTransform.GetRotationMatrix();
  } else {
    translation = (nodeTransform->GetTranslation());
    rotation    = (nodeTransform->GetRotationMatrix());
  }
  // create the translation
  Vector3D colT(scalor * translation[0],
                scalor * translation[1],
                scalor * translation[2]);
  Vector3D colX(rotation[0], rotation[3], rotation[6]);
  Vector3D colY(rotation[1], rotation[4], rotation[7]);
  Vector3D colZ(rotation[2], rotation[5], rotation[8]);

  // check if it's a box - always true ...
  TGeoBBox* box
      = dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
  // check if it's a trapezoid - unfortunately box is the base of everything
  TGeoTrd2* trapezoid
      = dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
  //
  if (boost::iequals(axes, "XYZ")) {
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signX = -1;
    if (islower(axes.at(1))) signY = -1;
    if (islower(axes.at(2))) signZ = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colX, colY, colZ, colT);
    if (trapezoid) {
      // bounds with x/y
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDx1(),
          scalor * trapezoid->GetDx2(),
          scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2()));
      // thickness
      m_thickness = scalor * trapezoid->GetDz();
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with x/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDX(), scalor * box->GetDY());
      // thickness
      m_thickness = scalor * box->GetDZ();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else if (boost::iequals(axes, "XZY")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signX = -1;
    if (islower(axes.at(1))) signZ = -1;
    if (islower(axes.at(2))) signY = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colX, colZ, colY, colT);
    if (trapezoid) {
      // bounds with x/z
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDx1(),
          scalor * trapezoid->GetDx2(),
          scalor * trapezoid->GetDz());
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with x/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDX(), scalor * box->GetDZ());
      // thickness
      m_thickness = scalor * box->GetDY();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }

  } else if (boost::iequals(axes, "YZX")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signY = -1;
    if (islower(axes.at(1))) signZ = -1;
    if (islower(axes.at(2))) signX = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colY, colZ, colX, colT);
    if (trapezoid) {
      // bounds with y/z
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDy1(),
          scalor * trapezoid->GetDy2(),
          scalor * trapezoid->GetDz());
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with y/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDY(), scalor * box->GetDZ());
      // thickness
      m_thickness = scalor * box->GetDX();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else if (boost::iequals(axes, "YXZ")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signY = -1;
    if (islower(axes.at(1))) signX = -1;
    if (islower(axes.at(2))) signZ = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colY, colX, colZ, colT);
    if (trapezoid) {
      // bounds with y/x
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDy1(),
          scalor * trapezoid->GetDy2(),
          scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2()));
      // thickness
      m_thickness = scalor * trapezoid->GetDz();
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with y/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDY(), scalor * box->GetDX());
      // thickness
      m_thickness = scalor * box->GetDZ();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else if (boost::iequals(axes, "ZYX")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signZ = -1;
    if (islower(axes.at(1))) signY = -1;
    if (islower(axes.at(2))) signX = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colZ, colY, colX, colT);
    if (trapezoid) {
      // bounds with z/y
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDz(),
          scalor * trapezoid->GetDz(),
          scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2()));
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with z/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDZ(), scalor * box->GetDY());
      // thickness
      m_thickness = scalor * box->GetDX();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else {
    // default is "ZXY"
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signZ = -1;
    if (islower(axes.at(1))) signX = -1;
    if (islower(axes.at(2))) signY = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colZ, colX, colY, colT);
    if (trapezoid) {
      // bounds with z/x
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDz(),
          scalor * trapezoid->GetDz(),
          scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2()));
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with z/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDZ(), scalor * box->GetDX());
      // thickness
      m_thickness = scalor * box->GetDY();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  }
}

Acts::TGeoDetectorElement::TGeoDetectorElement(const Identifier& identifier,
                                               const TGeoMatrix& transform,
                                               TGeoNode*         tGeoDetElement,
                                               const std::string& axes,
                                               double             scalor)
  : Acts::DetectorElementBase()
  , m_detElement(tGeoDetElement)
  , m_identifier(identifier)
  , m_thickness(0.)
{
  // get the placement and orientation in respect to its mother
  const Double_t* rotation    = nullptr;
  const Double_t* translation = nullptr;

  // the new combined translation
  translation = transform.GetTranslation();
  rotation    = transform.GetRotationMatrix();

  // create the translation
  Vector3D colT(scalor * translation[0],
                scalor * translation[1],
                scalor * translation[2]);
  Vector3D colX(rotation[0], rotation[3], rotation[6]);
  Vector3D colY(rotation[1], rotation[4], rotation[7]);
  Vector3D colZ(rotation[2], rotation[5], rotation[8]);

  // check if it's a box - always true ...
  TGeoBBox* box
      = dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
  // check if it's a trapezoid - unfortunately box is the base of everything
  TGeoTrd2* trapezoid
      = dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
  //
  if (boost::iequals(axes, "XYZ")) {
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signX = -1;
    if (islower(axes.at(1))) signY = -1;
    if (islower(axes.at(2))) signZ = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colX, colY, colZ, colT);
    if (trapezoid) {
      // bounds with x/y
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDx1(),
          scalor * trapezoid->GetDx2(),
          scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2()));
      // thickness
      m_thickness = scalor * trapezoid->GetDz();
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with x/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDX(), scalor * box->GetDY());
      // thickness
      m_thickness = scalor * box->GetDZ();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else if (boost::iequals(axes, "XZY")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signX = -1;
    if (islower(axes.at(1))) signZ = -1;
    if (islower(axes.at(2))) signY = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colX, colZ, colY, colT);
    if (trapezoid) {
      // bounds with x/z
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDx1(),
          scalor * trapezoid->GetDx2(),
          scalor * trapezoid->GetDz());
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with x/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDX(), scalor * box->GetDZ());
      // thickness
      m_thickness = scalor * box->GetDY();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }

  } else if (boost::iequals(axes, "YZX")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signY = -1;
    if (islower(axes.at(1))) signZ = -1;
    if (islower(axes.at(2))) signX = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colY, colZ, colX, colT);
    if (trapezoid) {
      // bounds with y/z
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDy1(),
          scalor * trapezoid->GetDy2(),
          scalor * trapezoid->GetDz());
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with y/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDY(), scalor * box->GetDZ());
      // thickness
      m_thickness = scalor * box->GetDX();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else if (boost::iequals(axes, "YXZ")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signY = -1;
    if (islower(axes.at(1))) signX = -1;
    if (islower(axes.at(2))) signZ = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colY, colX, colZ, colT);
    if (trapezoid) {
      // bounds with y/x
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDy1(),
          scalor * trapezoid->GetDy2(),
          scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2()));
      // thickness
      m_thickness = scalor * trapezoid->GetDz();
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with y/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDY(), scalor * box->GetDX());
      // thickness
      m_thickness = scalor * box->GetDZ();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else if (boost::iequals(axes, "ZYX")) {
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signZ = -1;
    if (islower(axes.at(1))) signY = -1;
    if (islower(axes.at(2))) signX = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colZ, colY, colX, colT);
    if (trapezoid) {
      // bounds with z/y
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDz(),
          scalor * trapezoid->GetDz(),
          scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2()));
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with z/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDZ(), scalor * box->GetDY());
      // thickness
      m_thickness = scalor * box->GetDX();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  } else {
    // default is "ZXY"
    // next possibility
    // get the sign of the axes
    int signX                      = 1;
    int signY                      = 1;
    int signZ                      = 1;
    if (islower(axes.at(0))) signZ = -1;
    if (islower(axes.at(1))) signX = -1;
    if (islower(axes.at(2))) signY = -1;
    // the transformation matrix
    colX *= signX;
    colY *= signY;
    colZ *= signZ;
    m_transform = std::make_shared<const Transform3D>(colZ, colX, colY, colT);
    if (trapezoid) {
      // bounds with z/x
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
          scalor * trapezoid->GetDz(),
          scalor * trapezoid->GetDz(),
          scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2()));
      // thickness
      m_thickness = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      // assign them
      m_bounds = trapezoidBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          trapezoidBounds, *this, m_identifier);
    } else {
      // bounds with z/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(
          scalor * box->GetDZ(), scalor * box->GetDX());
      // thickness
      m_thickness = scalor * box->GetDY();
      // assign them
      m_bounds = rectangleBounds;
      // create the surface
      m_surface = std::make_shared<const PlaneSurface>(
          rectangleBounds, *this, m_identifier);
    }
  }
}

Acts::TGeoDetectorElement::~TGeoDetectorElement()
{
}

std::shared_ptr<const Acts::DigitizationModule>
Acts::TGeoDetectorElement::digitizationModule() const
{
  return DetectorElementBase::digitizationModule();
}
