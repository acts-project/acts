// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DiamondVolumeBounds.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

DiamondVolumeBounds::DiamondVolumeBounds(double x1, double x2, double x3,
                                         double y1, double y2, double halez)

    : VolumeBounds() {
  m_values[eHalfLengthX1] = x1;
  m_values[eHalfLengthX2] = x2;
  m_values[eHalfLengthX3] = x3;
  m_values[eLengthY1] = y1;
  m_values[eLengthY2] = y2;
  m_values[eHalfLengthZ] = halez;
  m_values[eAlphaAngle] = std::numbers::pi - std::atan2(y1, (x2 - x1));
  m_values[eBetaAngle] = std::atan2(y2, (x2 - x3));

  checkConsistency();
  buildSurfaceBounds();
}

std::vector<double> DiamondVolumeBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

std::vector<OrientedSurface> DiamondVolumeBounds::orientedSurfaces(
    const Transform3& transform) const {
  std::vector<OrientedSurface> surfaces;
  surfaces.reserve(8);

  // (0) - At negative local z
  auto negZTransform = transform * Translation3(0., 0., -get(eHalfLengthZ));
  auto sf = Surface::makeShared<PlaneSurface>(negZTransform, m_FaceXYBounds);
  surfaces.push_back(OrientedSurface{std::move(sf), Direction::AlongNormal()});

  // (1) - At positive local z
  auto posZTransform = transform * Translation3(0., 0., get(eHalfLengthZ));
  sf = Surface::makeShared<PlaneSurface>(posZTransform, m_FaceXYBounds);
  surfaces.push_back(
      OrientedSurface{std::move(sf), Direction::OppositeNormal()});

  double posXOffset23 = 0.5 * get(eLengthY2) / std::tan(get(eBetaAngle));
  double posXOffset12 =
      0.5 * get(eLengthY1) / std::tan(std::numbers::pi - get(eAlphaAngle));

  // (2) - At negative x face yz12
  Vector3 nyz12Position(-get(eHalfLengthX1) - posXOffset12,
                        -0.5 * get(eLengthY1), 0.);
  auto nyz12Transform =
      transform * Translation3(nyz12Position) *
      AngleAxis3(-std::numbers::pi / 2. + get(eAlphaAngle), Vector3::UnitZ()) *
      s_planeYZ;
  sf = Surface::makeShared<PlaneSurface>(nyz12Transform, m_FaceYZ12Bounds);
  surfaces.push_back(OrientedSurface{std::move(sf), Direction::AlongNormal()});

  // (3) - At positive x face yz12
  Vector3 pyz12Position(get(eHalfLengthX1) + posXOffset12,
                        -0.5 * get(eLengthY1), 0.);
  auto pyz12Transform =
      transform * Translation3(pyz12Position) *
      AngleAxis3(-std::numbers::pi / 2. + get(eAlphaAngle), -Vector3::UnitZ()) *
      s_planeYZ;
  sf = Surface::makeShared<PlaneSurface>(pyz12Transform, m_FaceYZ12Bounds);
  surfaces.push_back(
      OrientedSurface{std::move(sf), Direction::OppositeNormal()});

  // (4) - At negative x face yz23
  Vector3 nyz23Position(-get(eHalfLengthX3) - posXOffset23,
                        0.5 * get(eLengthY2), 0.);
  auto nyz23Transform = transform * Translation3(nyz23Position) *
                        AngleAxis3(std::numbers::pi / 2. - get(eBetaAngle),
                                   Vector3(0., 0., -1.)) *
                        s_planeYZ;
  sf = Surface::makeShared<PlaneSurface>(nyz23Transform, m_FaceYZ23Bounds);
  surfaces.push_back(OrientedSurface{std::move(sf), Direction::AlongNormal()});

  // (5) - At positive x face yz23
  Vector3 pyz23Position(get(eHalfLengthX3) + posXOffset23, 0.5 * get(eLengthY2),
                        0.);
  auto pyz23Transform =
      transform * Translation3(pyz23Position) *
      AngleAxis3(std::numbers::pi / 2. - get(eBetaAngle), Vector3::UnitZ()) *
      s_planeYZ;

  sf = Surface::makeShared<PlaneSurface>(pyz23Transform, m_FaceYZ23Bounds);
  surfaces.push_back(
      OrientedSurface{std::move(sf), Direction::OppositeNormal()});

  // (6) - At negative y face zx
  auto nyTransform =
      transform * Translation3(0., -get(eLengthY1), 0.) * s_planeZX;
  sf = Surface::makeShared<PlaneSurface>(nyTransform, m_negYFaceZXBounds);
  surfaces.push_back(OrientedSurface{std::move(sf), Direction::AlongNormal()});

  // (7) - At positive y face zx
  auto pyTransform =
      transform * Translation3(0., get(eLengthY2), 0.) * s_planeZX;
  sf = Surface::makeShared<PlaneSurface>(pyTransform, m_posYFaceZXBounds);
  surfaces.push_back(
      OrientedSurface{std::move(sf), Direction::OppositeNormal()});

  return surfaces;
};

void DiamondVolumeBounds::buildSurfaceBounds() {
  m_FaceXYBounds = std::make_shared<Acts::DiamondBounds>(
      get(eHalfLengthX1), get(eHalfLengthX2), get(eHalfLengthX3),
      get(eLengthY1), get(eLengthY2));

  double dx23 =
      fastHypot((get(eHalfLengthX2) - get(eHalfLengthX3)), get(eLengthY2));
  double dx12 =
      fastHypot((get(eHalfLengthX2) - get(eHalfLengthX1)), get(eLengthY1));

  m_FaceYZ12Bounds =
      std::make_shared<Acts::RectangleBounds>(0.5 * dx12, get(eHalfLengthZ));
  m_FaceYZ23Bounds =
      std::make_shared<Acts::RectangleBounds>(0.5 * dx23, get(eHalfLengthZ));

  m_posYFaceZXBounds = std::make_shared<Acts::RectangleBounds>(
      get(eHalfLengthZ), get(eHalfLengthX3));
  m_negYFaceZXBounds = std::make_shared<Acts::RectangleBounds>(
      get(eHalfLengthZ), get(eHalfLengthX1));
}

Volume::BoundingBox DiamondVolumeBounds::boundingBox(
    const Transform3* trf, const Vector3& envelope,
    const Volume* entity) const {
  double maxX =
      std::max({get(eHalfLengthX1), get(eHalfLengthX2), get(eHalfLengthX3)});
  double maxY = std::max(get(eLengthY1), get(eLengthY2));
  double halez = get(eHalfLengthZ);

  std::array<Vector3, 8> vertices = {{{-maxX, -maxY, -halez},
                                      {+maxX, -maxY, -halez},
                                      {-maxX, +maxY, -halez},
                                      {+maxX, +maxY, -halez},
                                      {-maxX, -maxY, +halez},
                                      {+maxX, -maxY, +halez},
                                      {-maxX, +maxY, +halez},
                                      {+maxX, +maxY, +halez}}};

  Transform3 transform = Transform3::Identity();
  if (trf != nullptr) {
    transform = *trf;
  }

  Vector3 vmin = transform * vertices[0];
  Vector3 vmax = transform * vertices[0];

  for (std::size_t i = 1; i < 8; i++) {
    const Vector3 vtx = transform * vertices[i];
    vmin = vmin.cwiseMin(vtx);
    vmax = vmax.cwiseMax(vtx);
  }

  return {entity, vmin - envelope, vmax + envelope};
}
bool DiamondVolumeBounds::inside(const Vector3& pos, double tol) const {
  if (std::abs(pos.z()) > get(eHalfLengthZ) + tol) {
    return false;
  }

  if (std::abs(pos.y()) > std::max(get(eLengthY1), get(eLengthY2)) + tol) {
    return false;
  }

  Vector2 locPos(pos.x(), pos.y());
  return m_FaceXYBounds->inside(locPos,
                                BoundaryTolerance::AbsoluteEuclidean(tol));
}

void DiamondVolumeBounds::checkConsistency() noexcept(false) {
  if (get(eHalfLengthX1) <= 0. || get(eHalfLengthX2) <= 0. ||
      get(eHalfLengthX3) <= 0.) {
    throw std::invalid_argument(
        "DiamondVolumeBounds: invalid polygon parameters in x.");
  }
  if (get(eLengthY1) <= 0. || get(eLengthY2) <= 0.) {
    throw std::invalid_argument("DiamondVolumeBounds: invalid y extrusion.");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument("DiamondVolumeBounds: invalid y extrusion.");
  }

  // make sure this is a convex polygon - do not allow angles > 180 deg
  if (get(eHalfLengthX2) < get(eHalfLengthX1) ||
      get(eHalfLengthX2) < get(eHalfLengthX3)) {
    throw std::invalid_argument(
        "DiamondVolumeBounds: invalid polygon shape - not convex.");
  }
}

std::ostream& DiamondVolumeBounds::toStream(std::ostream& os) const {
  os << std::setiosflags(std::ios::fixed);
  os << std::setprecision(5);
  os << "DiamondVolumeBounds: (halfX1, halfX2, halfX3, halfY1, halfY2, "
        "halfZ) = ";
  os << "(" << get(eHalfLengthX1) << ", " << get(eHalfLengthX2) << ", "
     << get(eHalfLengthX3) << ", " << get(eLengthY1) << ", " << get(eLengthY2)
     << ", " << get(eHalfLengthZ) << ")";
  return os;
}

}  // namespace Acts
