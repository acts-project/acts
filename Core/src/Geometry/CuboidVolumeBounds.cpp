// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include <cmath>
#include <iostream>
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

Acts::CuboidVolumeBounds::CuboidVolumeBounds()
    : VolumeBounds(), m_valueStore(bv_length, 0.) {}

Acts::CuboidVolumeBounds::CuboidVolumeBounds(double halex, double haley,
                                             double halez)
    : VolumeBounds(),
      m_valueStore(bv_length, 0.),
      m_xyBounds(nullptr),
      m_yzBounds(nullptr),
      m_zxBounds(nullptr) {
  m_valueStore.at(bv_halfX) = halex;
  m_valueStore.at(bv_halfY) = haley;
  m_valueStore.at(bv_halfZ) = halez;

  m_xyBounds = faceXYRectangleBounds();
  m_yzBounds = faceYZRectangleBounds();
  m_zxBounds = faceZXRectangleBounds();
}

Acts::CuboidVolumeBounds::CuboidVolumeBounds(const CuboidVolumeBounds& bobo)
    : VolumeBounds(),
      m_valueStore(bobo.m_valueStore),
      m_xyBounds(bobo.m_xyBounds),
      m_yzBounds(bobo.m_yzBounds),
      m_zxBounds(bobo.m_zxBounds) {}

Acts::CuboidVolumeBounds::~CuboidVolumeBounds() = default;

Acts::CuboidVolumeBounds& Acts::CuboidVolumeBounds::operator=(
    const CuboidVolumeBounds& bobo) {
  if (this != &bobo) {
    m_valueStore = bobo.m_valueStore;
    m_xyBounds = bobo.m_xyBounds;
    m_yzBounds = bobo.m_yzBounds;
    m_zxBounds = bobo.m_zxBounds;
  }
  return *this;
}

Acts::SurfacePtrVector Acts::CuboidVolumeBounds::decomposeToSurfaces(
    const Transform3D* transformPtr) const {
  // the transform - apply when given
  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);
  const Transform3D* tTransform = nullptr;

  SurfacePtrVector rSurfaces;
  rSurfaces.reserve(6);
  // face surfaces xy -------------------------------------
  //   (1) - at negative local z
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(0., 1., 0.)) *
                      Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_xyBounds));
  //   (2) - at positive local z
  tTransform = new Transform3D(transform *
                               Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_xyBounds));
  // face surfaces yz -------------------------------------
  // transmute cyclical
  //   (3) - at negative local x
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(0., 0., 1.)) *
                      Translation3D(Vector3D(halflengthX(), 0., 0)) *
                      AngleAxis3D(0.5 * M_PI, Vector3D(0., 1., 0)) *
                      AngleAxis3D(0.5 * M_PI, Vector3D(0., 0., 1.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_yzBounds));
  //   (4) - at positive local x
  tTransform = new Transform3D(transform *
                               Translation3D(Vector3D(halflengthX(), 0., 0.)) *
                               AngleAxis3D(0.5 * M_PI, Vector3D(0., 1., 0.)) *
                               AngleAxis3D(0.5 * M_PI, Vector3D(0., 0., 1.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_yzBounds));
  // face surfaces zx -------------------------------------
  //   (5) - at negative local y
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.)) *
                      Translation3D(Vector3D(0., halflengthY(), 0.)) *
                      AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
                      AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_zxBounds));
  //   (6) - at positive local y
  tTransform = new Transform3D(transform *
                               Translation3D(Vector3D(0., halflengthY(), 0.)) *
                               AngleAxis3D(-0.5 * M_PI, Vector3D(0., 1., 0.)) *
                               AngleAxis3D(-0.5 * M_PI, Vector3D(1., 0., 0.)));
  rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
      std::shared_ptr<const Transform3D>(tTransform), m_zxBounds));
  // return the surfaces
  return rSurfaces;
}

std::shared_ptr<const Acts::RectangleBounds>
Acts::CuboidVolumeBounds::faceXYRectangleBounds() const {
  return std::make_shared<const RectangleBounds>(m_valueStore.at(bv_halfX),
                                                 m_valueStore.at(bv_halfY));
}

std::shared_ptr<const Acts::RectangleBounds>
Acts::CuboidVolumeBounds::faceYZRectangleBounds() const {
  return std::make_shared<const RectangleBounds>(m_valueStore.at(bv_halfY),
                                                 m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::RectangleBounds>
Acts::CuboidVolumeBounds::faceZXRectangleBounds() const {
  return std::make_shared<const RectangleBounds>(m_valueStore.at(bv_halfZ),
                                                 m_valueStore.at(bv_halfX));
}

// ostream operator overload
std::ostream& Acts::CuboidVolumeBounds::toStream(std::ostream& sl) const {
  return dumpT(sl);
}

Acts::Volume::BoundingBox Acts::CuboidVolumeBounds::boundingBox(
    const Acts::Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  Vector3D vmin(-halflengthX(), -halflengthY(), -halflengthZ());
  Vector3D vmax(halflengthX(), halflengthY(), halflengthZ());

  Volume::BoundingBox box(entity, vmin - envelope, vmax + envelope);
  return trf == nullptr ? box : box.transformed(*trf);
}
