// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Helpers.hpp"

template <typename value_t, size_t DIM, size_t SIDES>
template <size_t D, std::enable_if_t<D == 2, int>>
Acts::Frustum<value_t, DIM, SIDES>::Frustum(const vertex_type& origin,
                                            const vertex_type& dir,
                                            value_type opening_angle)
    : m_origin(origin) {
  using rotation_t = Eigen::Rotation2D<value_type>;

  static_assert(SIDES == 2, "2D frustum can only have 2 sides");
  assert(opening_angle < M_PI);

  translation_t translation(origin);
  value_type angle = VectorHelpers::phi(dir);
  Eigen::Rotation2D<value_type> rot(angle);

  value_type normal_angle = 0.5 * M_PI - 0.5 * opening_angle;
  vertex_type normal1 = rotation_t(normal_angle) * vertex_type::UnitX();
  vertex_type normal2 = rotation_t(-normal_angle) * vertex_type::UnitX();

  m_normals = {rot * vertex_type::UnitX(), rot * normal1, rot * normal2};
}

template <typename value_t, size_t DIM, size_t SIDES>
template <size_t D, std::enable_if_t<D == 3, int>>
Acts::Frustum<value_t, DIM, SIDES>::Frustum(const vertex_type& origin,
                                            const vertex_type& dir,
                                            value_type opening_angle)
    : m_origin(origin) {
  static_assert(SIDES > 2, "3D frustum must have 3 or more sides");
  assert(opening_angle < M_PI);
  using angle_axis_t = Eigen::AngleAxis<value_type>;

  const vertex_type ldir = vertex_type::UnitZ();
  const vertex_type lup = vertex_type::UnitX();

  transform_type transform;
  transform = (Eigen::Quaternion<value_type>().setFromTwoVectors(ldir, dir));

  m_normals[0] = ldir;

  const value_type phi_sep = 2 * M_PI / sides;
  transform_type rot;
  rot = angle_axis_t(phi_sep, ldir);

  value_type half_opening_angle = opening_angle / 2.;
  auto calculate_normal =
      [&ldir, &half_opening_angle](const vertex_type& out) -> vertex_type {
    const vertex_type tilt_axis = -1 * out.cross(ldir);
    return (-1 * (angle_axis_t(half_opening_angle, tilt_axis) * out))
        .normalized();
  };

  vertex_type current_outward = lup;
  m_normals[1] = calculate_normal(current_outward);

  for (size_t i = 1; i < sides; i++) {
    current_outward = rot * current_outward;
    m_normals[i + 1] = calculate_normal(current_outward);
  }

  for (auto& normal : m_normals) {
    normal = transform * normal;
  }
}

template <typename value_t, size_t DIM, size_t SIDES>
template <size_t D, std::enable_if_t<D == 3, int>>
void Acts::Frustum<value_t, DIM, SIDES>::draw(IVisualization& helper,
                                              value_type far_distance) const {
  static_assert(DIM == 3, "Drawing is only supported in 3D");

  // Iterate around normals, calculate cross with "far" plane
  // to get intersection lines.
  // Work in local reference frame of the frustum, and only convert to global
  // right before drawing.
  vertex_type far_normal = m_normals[0];  // far has same normal as pseudo-near
  vertex_type far_center = m_normals[0] * far_distance;
  std::array<std::pair<vertex_type, vertex_type>, SIDES> planeFarIXs;

  auto ixPlanePlane =
      [](const auto& n1, const auto& p1, const auto& n2,
         const auto& p2) -> std::pair<vertex_type, vertex_type> {
    const vertex_type m = n1.cross(n2).normalized();
    const double j = (n2.dot(p2 - p1)) / (n2.dot(n1.cross(m)));
    const vertex_type q = p1 + j * n1.cross(m);
    return {m, q};
  };

  auto ixLineLine = [](const auto& p1, const auto& d1, const auto& p2,
                       const auto& d2) -> vertex_type {
    return p1 + (((p2 - p1).cross(d2)).norm() / (d1.cross(d2)).norm()) * d1;
  };

  // skip i=0 <=> pseudo-near
  for (size_t i = 1; i < n_normals; i++) {
    const auto ixLine =
        ixPlanePlane(far_normal, far_center, m_normals[i], vertex_type::Zero());
    planeFarIXs.at(i - 1) = ixLine;
  }

  std::array<vertex_type, SIDES> points;

  for (size_t i = 0; i < std::size(planeFarIXs); i++) {
    size_t j = (i + 1) % std::size(planeFarIXs);
    const auto& l1 = planeFarIXs.at(i);
    const auto& l2 = planeFarIXs.at(j);
    const vertex_type ix =
        m_origin + ixLineLine(l1.second, l1.first, l2.second, l2.first);
    points.at(i) = ix;
  }

  for (size_t i = 0; i < std::size(points); i++) {
    size_t j = (i + 1) % std::size(points);
    helper.face(
        std::vector<vertex_type>({m_origin, points.at(i), points.at(j)}));
  }
}

template <typename value_t, size_t DIM, size_t SIDES>
template <size_t D, std::enable_if_t<D == 2, int>>
std::ostream& Acts::Frustum<value_t, DIM, SIDES>::svg(std::ostream& os,
                                                      value_type w,
                                                      value_type h,
                                                      value_type far_distance,
                                                      value_type unit) const {
  static_assert(DIM == 2, "SVG is only supported in 2D");

  vertex_type mid(w / 2., h / 2.);

  // set up transform for svg. +y is down, normally, and unit is pixels.
  // We flip the y axis, and scale up by `unit`.
  transform_type trf = transform_type::Identity();
  trf.translate(mid);
  trf = trf * Eigen::Scaling(vertex_type(1, -1));
  trf.scale(unit);

  std::array<std::string, 3> colors({"orange", "blue", "red"});

  auto draw_line = [&](const vertex_type& left_, const vertex_type& right_,
                       std::string color, size_t width) {
    vertex_type left = trf * left_;
    vertex_type right = trf * right_;
    os << "<line ";

    os << "x1=\"" << left.x() << "\" ";
    os << "y1=\"" << left.y() << "\" ";
    os << "x2=\"" << right.x() << "\" ";
    os << "y2=\"" << right.y() << "\" ";

    os << " stroke=\"" << color << "\" stroke-width=\"" << width << "\"/>\n";
  };

  auto draw_point = [&](const vertex_type& p_, std::string color, size_t r) {
    vertex_type p = trf * p_;
    os << "<circle ";
    os << "cx=\"" << p.x() << "\" cy=\"" << p.y() << "\" r=\"" << r << "\"";
    os << " fill=\"" << color << "\"";
    os << "/>\n";
  };

  using vec3 = ActsVector<value_type, 3>;
  auto ixLineLine = [](const vertex_type& p1_2, const vertex_type& d1_2,
                       const vertex_type& p2_2,
                       const vertex_type& d2_2) -> vertex_type {
    const vec3 p1(p1_2.x(), p1_2.y(), 0);
    const vec3 p2(p2_2.x(), p2_2.y(), 0);
    const vec3 d1(d1_2.x(), d1_2.y(), 0);
    const vec3 d2(d2_2.x(), d2_2.y(), 0);

    vec3 num = (p2 - p1).cross(d2);
    vec3 den = d1.cross(d2);

    value_type f = 1.;
    value_type dot = num.normalized().dot(den.normalized());
    if (std::abs(dot) - 1 < 1e-9 && dot < 0) {
      f = -1.;
    }

    const vec3 p = p1 + f * (num.norm() / den.norm()) * d1;
    assert(std::abs(p.z()) < 1e-9);
    return {p.x(), p.y()};
  };

  const vertex_type far_dir = {m_normals[0].y(), -m_normals[0].x()};
  const vertex_type far_point = m_normals[0] * far_distance;

  std::array<vertex_type, 2> points;

  for (size_t i = 1; i < n_normals; i++) {
    vertex_type plane_dir(m_normals[i].y(), -m_normals[i].x());

    const vertex_type ix = ixLineLine(far_point, far_dir, {0, 0}, plane_dir);
    draw_point(m_origin + ix, "black", 3);
    draw_line(m_origin, m_origin + ix, "black", 2);
    points.at(i - 1) = ix;
  }

  draw_line(m_origin + points.at(0), m_origin + points.at(1), "black", 2);

  draw_line(m_origin, m_origin + m_normals[0] * 2, colors[0], 3);

  draw_point({0, 0}, "yellow", 5);
  draw_point(m_origin, "green", 5);

  return os;
}

template <typename value_t, size_t DIM, size_t SIDES>
Acts::Frustum<value_t, DIM, SIDES>
Acts::Frustum<value_t, DIM, SIDES>::transformed(
    const transform_type& trf) const {
  const auto& rot = trf.rotation();

  std::array<vertex_type, n_normals> new_normals;
  for (size_t i = 0; i < n_normals; i++) {
    new_normals[i] = rot * m_normals[i];
  }

  return Frustum<value_t, DIM, SIDES>(trf * m_origin, std::move(new_normals));
}
