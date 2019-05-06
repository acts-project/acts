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

namespace Acts {

template <typename value_t, size_t DIM>
class Ray {
 public:
  using value_type = value_t;
  using vertex_type = ActsVector<value_t, DIM>;
  using vertex_array_type = Eigen::Array<value_t, DIM, 1>;
  using transform_type = Eigen::Transform<value_t, DIM, Eigen::Affine>;

  Ray(const vertex_type& origin, const vertex_type& dir);

  const vertex_type& origin() const { return m_origin; }

  const vertex_type& dir() const { return m_dir; }
  const vertex_array_type& idir() const { return m_idir; }

  Ray<value_t, DIM> transformed(const transform_type& trf) const;

  std::ostream& toStream(std::ostream& os) const;

  template <size_t D = DIM, std::enable_if_t<D == 3, int> = 0>
  void draw(IVisualization& helper, value_type far_distance = 10) const;

 private:
  vertex_type m_origin;
  vertex_type m_dir;
  vertex_array_type m_idir;
};

template <typename T, size_t D>
std::ostream& operator<<(std::ostream& os, const Ray<T, D>& ray) {
  ray.dump(os);
  return os;
}

using Ray3F = Ray<float, 3>;
using Ray3D = Ray<double, 3>;

}  // namespace Acts

#include "Acts/Utilities/Ray.ipp"
