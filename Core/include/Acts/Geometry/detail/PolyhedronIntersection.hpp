// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Geometry/Polyhedron.hpp>

namespace Acts::detail {

enum class TriIntersectionResult {
  eTriIntersectTrue = 0,
  eTriIntersectFalse,
  eTriCoplanar,
};

std::ostream &operator<<(std::ostream &os, TriIntersectionResult r);

// https://inria.hal.science/file/index/docid/72100/filename/RR-4488.pdf
TriIntersectionResult triangleIntersection3D(std::array<Vector3, 3> t1,
                                             std::array<Vector3, 3> t2,
                                             double epsilon = 1e-4);

bool triangleCoplanarIntersect(std::array<Vector3, 3> t1,
                               std::array<Vector3, 3> t2,
                               double epsilon = 1e-4);

bool generalTriangleIntersection(const std::array<Vector3, 3> &t1,
                                 const std::array<Vector3, 3> &t2,
                                 double epsilon = 1e-4);

bool polyhedronIntersection(const Acts::Polyhedron &p1,
                            const Acts::Polyhedron &p2, double epsilon = 1e-4);

}  // namespace Acts::detail
