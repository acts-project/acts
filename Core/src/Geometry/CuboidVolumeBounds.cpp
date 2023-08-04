// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CuboidVolumeBounds.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <type_traits>
#include <utility>

Acts::CuboidVolumeBounds::CuboidVolumeBounds(double halex, double haley,
                                             double halez)
    : VolumeBounds(), m_values({halex, haley, halez}) {
  checkConsistency();
  buildSurfaceBounds();
}

Acts::CuboidVolumeBounds::CuboidVolumeBounds(const CuboidVolumeBounds& bobo)
    : VolumeBounds(),
      m_values(bobo.m_values),
      m_xyBounds(bobo.m_xyBounds),
      m_yzBounds(bobo.m_yzBounds),
      m_zxBounds(bobo.m_zxBounds) {}

Acts::CuboidVolumeBounds& Acts::CuboidVolumeBounds::operator=(
    const CuboidVolumeBounds& bobo) {
  if (this != &bobo) {
    m_values = bobo.m_values;
    m_xyBounds = bobo.m_xyBounds;
    m_yzBounds = bobo.m_yzBounds;
    m_zxBounds = bobo.m_zxBounds;
  }
  return *this;
}

Acts::OrientedSurfaces Acts::CuboidVolumeBounds::orientedSurfaces(
    const Transform3& transform) const {
  OrientedSurfaces oSurfaces;
  oSurfaces.reserve(6);
  // Face surfaces xy -------------------------------------
  //   (1) - at negative local z
  auto sf = Surface::makeShared<PlaneSurface>(
      transform * Translation3(0., 0., -get(eHalfLengthZ)), m_xyBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), Direction::Positive));
  //   (2) - at positive local z
  sf = Surface::makeShared<PlaneSurface>(
      transform * Translation3(0., 0., get(eHalfLengthZ)), m_xyBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), Direction::Negative));
  // Face surfaces yz -------------------------------------
  //   (3) - at negative local x
  sf = Surface::makeShared<PlaneSurface>(
      transform * Translation3(-get(eHalfLengthX), 0., 0.) * s_planeYZ,
      m_yzBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), Direction::Positive));
  //   (4) - at positive local x
  sf = Surface::makeShared<PlaneSurface>(
      transform * Translation3(get(eHalfLengthX), 0., 0.) * s_planeYZ,
      m_yzBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), Direction::Negative));
  // Face surfaces zx -------------------------------------
  //   (5) - at negative local y
  sf = Surface::makeShared<PlaneSurface>(
      transform * Translation3(0., -get(eHalfLengthY), 0.) * s_planeZX,
      m_zxBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), Direction::Positive));
  //   (6) - at positive local y
  sf = Surface::makeShared<PlaneSurface>(
      transform * Translation3(0., get(eHalfLengthY), 0.) * s_planeZX,
      m_zxBounds);
  oSurfaces.push_back(OrientedSurface(std::move(sf), Direction::Negative));

  return oSurfaces;
}

std::ostream& Acts::CuboidVolumeBounds::toStream(std::ostream& sl) const {
  return dumpT(sl);
}

Acts::Volume::BoundingBox Acts::CuboidVolumeBounds::boundingBox(
    const Acts::Transform3* trf, const Vector3& envelope,
    const Volume* entity) const {
  Vector3 vmin(-get(eHalfLengthX), -get(eHalfLengthY), -get(eHalfLengthZ));
  Vector3 vmax(get(eHalfLengthX), get(eHalfLengthY), get(eHalfLengthZ));

  Volume::BoundingBox box(entity, vmin - envelope, vmax + envelope);
  return trf == nullptr ? box : box.transformed(*trf);
}

void Acts::CuboidVolumeBounds::buildSurfaceBounds() {
  m_xyBounds = std::make_shared<const RectangleBounds>(get(eHalfLengthX),
                                                       get(eHalfLengthY));
  m_yzBounds = std::make_shared<const RectangleBounds>(get(eHalfLengthY),
                                                       get(eHalfLengthZ));
  m_zxBounds = std::make_shared<const RectangleBounds>(get(eHalfLengthZ),
                                                       get(eHalfLengthX));
}
