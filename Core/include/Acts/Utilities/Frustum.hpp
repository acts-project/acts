// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/IVisualization.hpp"

#include <ostream>

namespace Acts {

template <typename value_t, size_t DIM, size_t SIDES>
class Frustum {
  using translation_t = Eigen::Translation<value_t, DIM>;

  static constexpr size_t n_normals = SIDES + 1;

 public:
  using transform_type = Eigen::Transform<value_t, DIM, Eigen::Affine>;
  using value_type = value_t;
  using vertex_type = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;

  static constexpr size_t dim = DIM;
  static constexpr size_t sides = SIDES;

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  Frustum(const vertex_type& origin, const vertex_type& dir,
          value_type opening_angle);

  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  Frustum(const vertex_type& origin, const vertex_type& dir,
          value_type opening_angle);

  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  void draw(IVisualization& helper, value_type far_distance = 10) const;

  template <size_t D = DIM, std::enable_if_t<D == 2, int> = 0>
  std::ostream& svg(std::ostream& os, value_type w, value_type h,
                    value_type far_distance = 1, value_type unit = 20.) const;

  const vertex_type& origin() const { return m_origin; }

  const vertex_type& dir() const { return m_normals[0]; }

  const std::array<vertex_type, SIDES + 1>& normals() const {
    return m_normals;
  }

  Frustum<value_t, DIM, SIDES> transformed(const transform_type& trf) const;

 private:
  // private constructor with pre-calculated normals
  Frustum(const vertex_type& origin, std::array<vertex_type, SIDES + 1> normals)
      : m_origin(origin), m_normals(std::move(normals)) {}

  vertex_type m_origin;
  // need one more for direction we're facing
  std::array<vertex_type, SIDES + 1> m_normals;
};

}  // namespace Acts

#include "Acts/Utilities/Frustum.ipp"
