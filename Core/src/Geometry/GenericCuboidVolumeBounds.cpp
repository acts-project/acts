// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <array>
#include <ostream>

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/Definitions.hpp"

Acts::GenericCuboidVolumeBounds::GenericCuboidVolumeBounds(
    const std::array<Acts::Vector3D, 8>& vertices) noexcept(false)
    : m_vertices(vertices) {
  construct();
}

Acts::GenericCuboidVolumeBounds::GenericCuboidVolumeBounds(
    const std::array<double, GenericCuboidVolumeBounds::eSize>&
        values) noexcept(false)
    : m_vertices() {
  for (size_t iv = 0; iv < 8; ++iv) {
    m_vertices[iv] =
        Vector3D(values[iv * 3], values[iv * 3 + 1], values[iv * 3 + 2]);
  }
  construct();
}

bool Acts::GenericCuboidVolumeBounds::inside(const Acts::Vector3D& gpos,
                                             double tol) const {
  constexpr std::array<size_t, 6> vtxs = {0, 4, 0, 1, 2, 1};
  // needs to be on same side, get ref
  bool ref = std::signbit((gpos - m_vertices[vtxs[0]]).dot(m_normals[0]));
  for (size_t i = 1; i < 6; i++) {
    double dot = (gpos - m_vertices[vtxs[i]]).dot(m_normals[i]);
    if (std::signbit(dot) != ref) {
      // technically outside, but how far?
      if (std::abs(dot) > tol) {
        // distance greater than tol
        return false;
      }
      // distance smaller than tol, ignore
    }
  }
  return true;
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::GenericCuboidVolumeBounds::decomposeToSurfaces(
    const Acts::Transform3D* transform) const {
  std::vector<std::shared_ptr<const Acts::Surface>> surfaces;

  // approximate cog of the volume
  Vector3D cog(0, 0, 0);

  for (size_t i = 0; i < 8; i++) {
    cog += m_vertices[i];
  }

  cog *= 0.125;  // 1/8.

  auto make_surface = [&](const auto& a, const auto& b, const auto& c,
                          const auto& d) {
    // calculate centroid of these points
    Vector3D ctrd = (a + b + c + d) / 4.;
    // create normal
    const Vector3D ab = b - a, ac = c - a;
    Vector3D normal = ab.cross(ac).normalized();

    if ((cog - d).dot(normal) > 0) {
      // normal points inwards, flip normal
      normal *= -1.;
    }
    // get rid of -0 values if present
    normal += Vector3D::Zero();

    // normal should point away from volume center now

    // build transform from z unit to normal
    // z is normal in local coordinates
    // Volume local to surface local
    Transform3D vol2srf;
    vol2srf = (Eigen::Quaternion<double>().setFromTwoVectors(
        normal, Vector3D::UnitZ()));

    vol2srf = vol2srf * Translation3D(-ctrd);

    // now calculate position of vertices in surface local frame
    Vector3D a_l, b_l, c_l, d_l;
    a_l = vol2srf * a;
    b_l = vol2srf * b;
    c_l = vol2srf * c;
    d_l = vol2srf * d;

    std::vector<Vector2D> vertices({{a_l.x(), a_l.y()},
                                    {b_l.x(), b_l.y()},
                                    {c_l.x(), c_l.y()},
                                    {d_l.x(), d_l.y()}});

    auto polyBounds = std::make_shared<const ConvexPolygonBounds<4>>(vertices);

    auto srfTrf = std::make_shared<Transform3D>(vol2srf.inverse());
    if (transform != nullptr) {
      *srfTrf = (*transform) * (*srfTrf);
    }

    auto srf = Surface::makeShared<PlaneSurface>(std::move(srfTrf), polyBounds);

    surfaces.push_back(std::move(srf));
  };

  make_surface(m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
  make_surface(m_vertices[4], m_vertices[5], m_vertices[6], m_vertices[7]);
  make_surface(m_vertices[0], m_vertices[3], m_vertices[7], m_vertices[4]);
  make_surface(m_vertices[1], m_vertices[2], m_vertices[6], m_vertices[5]);
  make_surface(m_vertices[2], m_vertices[3], m_vertices[7], m_vertices[6]);
  make_surface(m_vertices[1], m_vertices[0], m_vertices[4], m_vertices[5]);

  return surfaces;
}

std::ostream& Acts::GenericCuboidVolumeBounds::toStream(
    std::ostream& sl) const {
  sl << "Acts::GenericCuboidVolumeBounds: vertices (x, y, z) =\n";
  for (size_t i = 0; i < 8; i++) {
    if (i > 0) {
      sl << ",\n";
    }
    sl << "[" << m_vertices[i].transpose() << "]";
  }
  return sl;
}

void Acts::GenericCuboidVolumeBounds::construct() noexcept(false) {
  // calculate approximate center of gravity first, so we can make sure
  // the normals point inwards
  Vector3D cog(0, 0, 0);

  for (size_t i = 0; i < 8; i++) {
    cog += m_vertices[i];
  }

  cog *= 0.125;  // 1/8.

  size_t idx = 0;

  auto handle_face = [&](const auto& a, const auto& b, const auto& c,
                         const auto& d) {
    // we assume a b c d to be counter clockwise
    const Vector3D ab = b - a, ac = c - a;
    Vector3D normal = ab.cross(ac).normalized();

    if ((cog - a).dot(normal) < 0) {
      // normal points outwards, flip normal
      normal *= -1.;
    }

    // get rid of -0 values if present
    normal += Vector3D::Zero();

    // check if d is on the surface
    if (std::abs((a - d).dot(normal)) > 1e-6) {
      throw(std::invalid_argument(
          "GenericCuboidBounds: Four points do not lie on the same plane!"));
    }

    m_normals[idx] = normal;
    idx++;
  };

  // handle faces
  handle_face(m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
  handle_face(m_vertices[4], m_vertices[5], m_vertices[6], m_vertices[7]);
  handle_face(m_vertices[0], m_vertices[3], m_vertices[7], m_vertices[4]);
  handle_face(m_vertices[1], m_vertices[2], m_vertices[6], m_vertices[5]);
  handle_face(m_vertices[2], m_vertices[3], m_vertices[7], m_vertices[6]);
  handle_face(m_vertices[1], m_vertices[0], m_vertices[4], m_vertices[5]);
}

std::vector<double> Acts::GenericCuboidVolumeBounds::values() const {
  std::vector<double> rvalues;
  rvalues.reserve(eSize);
  for (size_t iv = 0; iv < 8; ++iv) {
    for (size_t ic = 0; ic < 3; ++ic) {
      rvalues.push_back(m_vertices[iv][ic]);
    }
  }
  return rvalues;
}

Acts::Volume::BoundingBox Acts::GenericCuboidVolumeBounds::boundingBox(
    const Acts::Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  Vector3D vmin, vmax;

  Transform3D transform = Transform3D::Identity();
  if (trf != nullptr) {
    transform = *trf;
  }

  vmin = transform * m_vertices[0];
  vmax = transform * m_vertices[0];

  for (size_t i = 1; i < 8; i++) {
    Vector3D vtx = transform * m_vertices[i];
    vmin = vmin.cwiseMin(vtx);
    vmax = vmax.cwiseMax(vtx);
  }

  return {entity, vmin - envelope, vmax + envelope};
}

void Acts::GenericCuboidVolumeBounds::draw(IVisualization& helper,
                                           const Transform3D& transform) const {
  auto draw_face = [&](const auto& a, const auto& b, const auto& c,
                       const auto& d) {
    helper.face(std::vector<Vector3D>(
        {transform * a, transform * b, transform * c, transform * d}));
  };

  draw_face(m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
  draw_face(m_vertices[4], m_vertices[5], m_vertices[6], m_vertices[7]);
  draw_face(m_vertices[0], m_vertices[3], m_vertices[7], m_vertices[4]);
  draw_face(m_vertices[1], m_vertices[2], m_vertices[6], m_vertices[5]);
  draw_face(m_vertices[2], m_vertices[3], m_vertices[7], m_vertices[6]);
  draw_face(m_vertices[1], m_vertices[0], m_vertices[4], m_vertices[5]);
}
