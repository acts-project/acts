// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <utility>

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelpers.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "TGeoBBox.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, TGeoNode* tGeoDetElement,
    const TGeoMatrix* mGlobal, const std::string& axes, double scalor,
    bool isDisc, std::shared_ptr<const Acts::ISurfaceMaterial> material,
    std::shared_ptr<const Acts::DigitizationModule> digitizationModule)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(tGeoDetElement),
      m_identifier(identifier),
      m_digitizationModule(std::move(digitizationModule)) {
  using namespace TGeoPrimitivesHelpers;

  // create temporary local non const surface (to allow setting the material)
  std::shared_ptr<Surface> surface = nullptr;
  // get the placement and orientation in respect to its mother
  const TGeoMatrix* nodeTransform = (m_detElement->GetMatrix());
  const Double_t* rotation = nullptr;
  const Double_t* translation = nullptr;

  if (mGlobal != nullptr) {
    // the new combined translation
    // TGeoHMatrix nTransform = (*mGlobal) * (*nodeTransform);
    TGeoHMatrix nTransform =
        TGeoCombiTrans(*mGlobal) * TGeoCombiTrans(*nodeTransform);
    // TGeoHMatrix nTransform = mGlobal->operator*(*nodeTransform);
    std::string nName = tGeoDetElement->GetName();
    std::string suffix = "_transform";
    nTransform.SetName((nName + suffix).c_str());
    translation = nTransform.GetTranslation();
    rotation = nTransform.GetRotationMatrix();
  } else {
    translation = (nodeTransform->GetTranslation());
    rotation = (nodeTransform->GetRotationMatrix());
  }
  // Simply call the construct method
  construct(rotation, translation, axes, scalor, isDisc, material);
}

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoMatrix& transform,
    TGeoNode* tGeoDetElement, const std::string& axes, double scalor,
    bool isDisc, std::shared_ptr<const Acts::ISurfaceMaterial> material,
    std::shared_ptr<const Acts::DigitizationModule> digitizationModule)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(tGeoDetElement),
      m_identifier(identifier),
      m_digitizationModule(std::move(digitizationModule)) {
  // get the placement and orientation in respect to its mother
  const Double_t* rotation = transform.GetRotationMatrix();
  const Double_t* translation = transform.GetTranslation();
  // Simply call the construct method
  construct(rotation, translation, axes, scalor, isDisc, material);
}

void Acts::TGeoDetectorElement::construct(
    const Double_t* rotation, const Double_t* translation,
    const std::string& axes, double scalor, bool isDisc,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  using namespace TGeoPrimitivesHelpers;

  // create temporary local non const surface (to allow setting the material)
  std::shared_ptr<Surface> surface = nullptr;

  // create the translation
  Vector3D colT(scalor * translation[0], scalor * translation[1],
                scalor * translation[2]);
  Vector3D colX(rotation[0], rotation[3], rotation[6]);
  Vector3D colY(rotation[1], rotation[4], rotation[7]);
  Vector3D colZ(rotation[2], rotation[5], rotation[8]);

  // check if it's a box - always true ...
  TGeoBBox* box =
      dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
  // check if it's a trapezoid - unfortunately box is the base of everything
  TGeoTrd2* trapezoid =
      dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
  // check if it's a tube segment
  TGeoTubeSeg* tube =
      dynamic_cast<TGeoTubeSeg*>(m_detElement->GetVolume()->GetShape());
  if (tube != nullptr) {
    m_transform = std::make_shared<const Transform3D>(
        makeTransform(colX, colY, colZ, colT));
    double rMin = tube->GetRmin() * scalor;
    double rMax = tube->GetRmax() * scalor;
    double halfZ = tube->GetDz() * scalor;

    if (isDisc) {
      // create disc surface
      m_thickness = halfZ;
      auto radialBounds = std::make_shared<const RadialBounds>(rMin, rMax);
      m_bounds = radialBounds;
      surface = Surface::makeShared<DiscSurface>(radialBounds, *this);
    } else {
      // create a cylinder surface
      m_thickness = std::fabs(rMax - rMin);
      double radius = (rMin + rMax) * 0.5;
      auto cylinderBounds =
          std::make_shared<const CylinderBounds>(radius, halfZ);
      m_bounds = cylinderBounds;
      surface = Surface::makeShared<CylinderSurface>(cylinderBounds, *this);
    }
  } else {
    if (boost::iequals(axes, "XYZ")) {
      // get the sign of the axes
      int signX = 1;
      int signY = 1;
      int signZ = 1;
      if (islower(axes.at(0)) != 0) {
        signX = -1;
      }
      if (islower(axes.at(1)) != 0) {
        signY = -1;
      }
      if (islower(axes.at(2)) != 0) {
        signZ = -1;
      }
      // the transformation matrix
      colX *= signX;
      colY *= signY;
      colZ *= signZ;
      m_transform = std::make_shared<const Transform3D>(
          makeTransform(colX, colY, colZ, colT));
      if (trapezoid != nullptr) {
        // bounds with x/y
        auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
            scalor * trapezoid->GetDx1(), scalor * trapezoid->GetDx2(),
            scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2()));
        // thickness
        m_thickness = scalor * trapezoid->GetDz();
        // assign them
        m_bounds = trapezoidBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(trapezoidBounds, *this);
      } else {
        // bounds with x/y
        auto rectangleBounds = std::make_shared<const RectangleBounds>(
            scalor * box->GetDX(), scalor * box->GetDY());
        // thickness
        m_thickness = scalor * box->GetDZ();
        // assign them
        m_bounds = rectangleBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(rectangleBounds, *this);
      }
    } else if (boost::iequals(axes, "XZY")) {
      // next possibility
      // get the sign of the axes
      int signX = 1;
      int signY = 1;
      int signZ = 1;
      if (islower(axes.at(0)) != 0) {
        signX = -1;
      }
      if (islower(axes.at(1)) != 0) {
        signZ = -1;
      }
      if (islower(axes.at(2)) != 0) {
        signY = -1;
      }
      // the transformation matrix
      colX *= signX;
      colY *= signY;
      colZ *= signZ;
      m_transform = std::make_shared<const Transform3D>(
          makeTransform(colX, colZ, colY, colT));
      if (trapezoid != nullptr) {
        // bounds with x/z
        auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
            scalor * trapezoid->GetDx1(), scalor * trapezoid->GetDx2(),
            scalor * trapezoid->GetDz());
        // thickness
        m_thickness =
            scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
        // assign them
        m_bounds = trapezoidBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(trapezoidBounds, *this);
      } else {
        // bounds with x/z
        auto rectangleBounds = std::make_shared<const RectangleBounds>(
            scalor * box->GetDX(), scalor * box->GetDZ());
        // thickness
        m_thickness = scalor * box->GetDY();
        // assign them
        m_bounds = rectangleBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(rectangleBounds, *this);
      }

    } else if (boost::iequals(axes, "YZX")) {
      // next possibility
      // get the sign of the axes
      int signX = 1;
      int signY = 1;
      int signZ = 1;
      if (islower(axes.at(0)) != 0) {
        signY = -1;
      }
      if (islower(axes.at(1)) != 0) {
        signZ = -1;
      }
      if (islower(axes.at(2)) != 0) {
        signX = -1;
      }
      // the transformation matrix
      colX *= signX;
      colY *= signY;
      colZ *= signZ;
      m_transform = std::make_shared<const Transform3D>(
          makeTransform(colY, colZ, colX, colT));
      if (trapezoid != nullptr) {
        // bounds with y/z
        auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
            scalor * trapezoid->GetDy1(), scalor * trapezoid->GetDy2(),
            scalor * trapezoid->GetDz());
        // thickness
        m_thickness =
            scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
        // assign them
        m_bounds = trapezoidBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(trapezoidBounds, *this);
      } else {
        // bounds with y/z
        auto rectangleBounds = std::make_shared<const RectangleBounds>(
            scalor * box->GetDY(), scalor * box->GetDZ());
        // thickness
        m_thickness = scalor * box->GetDX();
        // assign them
        m_bounds = rectangleBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(rectangleBounds, *this);
      }
    } else if (boost::iequals(axes, "YXZ")) {
      // next possibility
      // get the sign of the axes
      int signX = 1;
      int signY = 1;
      int signZ = 1;
      if (islower(axes.at(0)) != 0) {
        signY = -1;
      }
      if (islower(axes.at(1)) != 0) {
        signX = -1;
      }
      if (islower(axes.at(2)) != 0) {
        signZ = -1;
      }
      // the transformation matrix
      colX *= signX;
      colY *= signY;
      colZ *= signZ;
      m_transform = std::make_shared<const Transform3D>(
          makeTransform(colY, colX, colZ, colT));
      if (trapezoid != nullptr) {
        // bounds with y/x
        auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
            scalor * trapezoid->GetDy1(), scalor * trapezoid->GetDy2(),
            scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2()));
        // thickness
        m_thickness = scalor * trapezoid->GetDz();
        // assign them
        m_bounds = trapezoidBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(trapezoidBounds, *this);
      } else {
        // bounds with y/x
        auto rectangleBounds = std::make_shared<const RectangleBounds>(
            scalor * box->GetDY(), scalor * box->GetDX());
        // thickness
        m_thickness = scalor * box->GetDZ();
        // assign them
        m_bounds = rectangleBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(rectangleBounds, *this);
      }
    } else if (boost::iequals(axes, "ZYX")) {
      // next possibility
      // get the sign of the axes
      int signX = 1;
      int signY = 1;
      int signZ = 1;
      if (islower(axes.at(0)) != 0) {
        signZ = -1;
      }
      if (islower(axes.at(1)) != 0) {
        signY = -1;
      }
      if (islower(axes.at(2)) != 0) {
        signX = -1;
      }
      // the transformation matrix
      colX *= signX;
      colY *= signY;
      colZ *= signZ;
      m_transform = std::make_shared<const Transform3D>(
          makeTransform(colZ, colY, colX, colT));
      if (trapezoid != nullptr) {
        // bounds with z/y
        auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
            scalor * trapezoid->GetDz(), scalor * trapezoid->GetDz(),
            scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2()));
        // thickness
        m_thickness =
            scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
        // assign them
        m_bounds = trapezoidBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(trapezoidBounds, *this);
      } else {
        // bounds with z/y
        auto rectangleBounds = std::make_shared<const RectangleBounds>(
            scalor * box->GetDZ(), scalor * box->GetDY());
        // thickness
        m_thickness = scalor * box->GetDX();
        // assign them
        m_bounds = rectangleBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(rectangleBounds, *this);
      }
    } else {
      // default is "ZXY"
      // next possibility
      // get the sign of the axes
      int signX = 1;
      int signY = 1;
      int signZ = 1;
      if (islower(axes.at(0)) != 0) {
        signZ = -1;
      }
      if (islower(axes.at(1)) != 0) {
        signX = -1;
      }
      if (islower(axes.at(2)) != 0) {
        signY = -1;
      }
      // the transformation matrix
      colX *= signX;
      colY *= signY;
      colZ *= signZ;
      m_transform = std::make_shared<const Transform3D>(
          makeTransform(colZ, colX, colY, colT));
      if (trapezoid != nullptr) {
        // bounds with z/x
        auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(
            scalor * trapezoid->GetDz(), scalor * trapezoid->GetDz(),
            scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2()));
        // thickness
        m_thickness =
            scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
        // assign them
        m_bounds = trapezoidBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(trapezoidBounds, *this);
      } else {
        // bounds with z/x
        auto rectangleBounds = std::make_shared<const RectangleBounds>(
            scalor * box->GetDZ(), scalor * box->GetDX());
        // thickness
        m_thickness = scalor * box->GetDY();
        // assign them
        m_bounds = rectangleBounds;
        // create the surface
        surface = Surface::makeShared<PlaneSurface>(rectangleBounds, *this);
      }
    }
  }
  // set the asscoiated material (non const method)
  if (surface) {
    surface->assignSurfaceMaterial(std::move(material));
  }
  // set the const member surface
  m_surface = surface;
}

Acts::TGeoDetectorElement::~TGeoDetectorElement() = default;
