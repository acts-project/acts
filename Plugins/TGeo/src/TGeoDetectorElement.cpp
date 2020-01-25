// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <utility>

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelpers.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
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
    TGeoHMatrix nTransform =
        TGeoCombiTrans(*mGlobal) * TGeoCombiTrans(*nodeTransform);
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
  construct(rotation, translation, axes, scalor, isDisc, std::move(material));
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
  construct(rotation, translation, axes, scalor, isDisc, std::move(material));
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

  auto sensor = m_detElement->GetVolume();

  // Special test for composite shape of silicon
  auto tgShape = sensor->GetShape();
  auto compShape = dynamic_cast<TGeoCompositeShape*>(tgShape);
  if (compShape != nullptr) {
    auto interNode = dynamic_cast<TGeoIntersection*>(compShape->GetBoolNode());
    if (interNode) {
      auto baseTube = dynamic_cast<TGeoTubeSeg*>(interNode->GetLeftShape());
      if (baseTube) {
        m_transform = std::make_shared<const Transform3D>(
            makeTransform(colX, colY, colZ, colT));
        double rMin = baseTube->GetRmin() * scalor;
        double rMax = baseTube->GetRmax() * scalor;
        auto maskShape = dynamic_cast<TGeoArb8*>(interNode->GetRightShape());
        if (maskShape) {
          auto rMatrix = interNode->GetRightMatrix();
          // get the placement and orientation in respect to its mother
          auto rRotation = rMatrix->GetRotationMatrix();
          auto rTranslation = rMatrix->GetTranslation();

          // create the translation
          Vector3D rColT(scalor * rTranslation[0], scalor * rTranslation[1],
                         scalor * rTranslation[2]);
          Vector3D rColX(rRotation[0], rRotation[3], rRotation[6]);
          Vector3D rColY(rRotation[1], rRotation[4], rRotation[7]);
          Vector3D rColZ(rRotation[2], rRotation[5], rRotation[8]);

          auto rTransform = makeTransform(rColX, rColY, rColZ, rColT);

          // Helper method for line line intersection
          auto det = [](double a, double b, double c, double d) -> double {
            return a * d - b * c;
          };

          auto lineLineIntersect = [&](double x1, double y1,  // Line 1 start
                                       double x2, double y2,  // Line 1 end
                                       double x3, double y3,  // Line 2 start
                                       double x4, double y4,  // Line 2 end
                                       double& ixOut, double& iyOut) -> bool {
            // http://mathworld.wolfram.com/Line-LineIntersection.html

            double detL1 = det(x1, y1, x2, y2);
            double detL2 = det(x3, y3, x4, y4);
            double x1mx2 = x1 - x2;
            double x3mx4 = x3 - x4;
            double y1my2 = y1 - y2;
            double y3my4 = y3 - y4;

            double xnom = det(detL1, x1mx2, detL2, x3mx4);
            double ynom = det(detL1, y1my2, detL2, y3my4);
            double denom = det(x1mx2, y1my2, x3mx4, y3my4);
            if (denom == 0.0) {
              ixOut = NAN;
              iyOut = NAN;
              return false;
            }

            ixOut = xnom / denom;
            iyOut = ynom / denom;
            return true;
          };

          auto vertices = maskShape->GetVertices();
          std::array<Vector3D, 4> vxyz;
          for (unsigned int iv = 0; iv < 4; ++iv) {
            Vector3D v = (*m_transform) * rTransform *
                         Vector3D(scalor * vertices[2 * iv],
                                  scalor * vertices[2 * iv + 1], 0.);
            vxyz[iv] = v;
          }

          double ix, iy;
          if (lineLineIntersect(vxyz[0].x(), vxyz[0].y(), vxyz[1].x(),
                                vxyz[1].y(), vxyz[2].x(), vxyz[2].y(),
                                vxyz[3].x(), vxyz[3].y(), ix, iy)) {

            // rMin / rMax are in module system
            // phiMin / rMax are in strip sytem                      
            Vector3D ssOffset(ix,iy,0.);
            double phiMin = VectorHelpers::phi(vxyz[3]+ssOffset);
            double phiMax = VectorHelpers::phi(vxyz[0]+ssOffset);
            // Create the bounds
            auto annulusBounds 
              =  std::make_shared<const AnnulusBounds>(rMin, rMax,
                                                       phiMin, phiMax,
                                                       Vector2D(ix,iy));   

            m_thickness = maskShape->GetDZ() * scalor;
            m_bounds = annulusBounds;
            surface = Surface::makeShared<DiscSurface>(annulusBounds, *this);
          }
        }
      }
    }
  }

  if (surface == nullptr) {
    // check if it's a box - always true ...
    TGeoBBox* box = dynamic_cast<TGeoBBox*>(tgShape);
    // check if it's a trapezoid - unfortunately box is the base of everything
    TGeoTrd2* trapezoid = dynamic_cast<TGeoTrd2*>(tgShape);
    // check if it's a tube segment
    TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>(tgShape);
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
  }
  // set the asscoiated material (non const method)
  if (surface != nullptr) {
    surface->assignSurfaceMaterial(std::move(material));
  }
  // set the const member surface
  m_surface = surface;
}

Acts::TGeoDetectorElement::~TGeoDetectorElement() = default;
