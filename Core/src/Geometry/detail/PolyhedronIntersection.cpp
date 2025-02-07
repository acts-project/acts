// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Geometry/Polyhedron.hpp>
#include <Acts/Geometry/detail/PolyhedronIntersection.hpp>

#include <boost/geometry.hpp>

namespace Acts::detail {

int fourPointDeterminant(const Vector3 &a, const Vector3 &b, const Vector3 &c,
                         const Vector3 &d, double epsilon = 1.e-8) {
  SquareMatrix3 m;
  m.row(0) = a - d;
  m.row(1) = b - d;
  m.row(2) = c - d;
  auto det = m.determinant();
  if (det > epsilon) {
    return 1;
  } else if (det < -epsilon) {
    return -1;
  } else {
    return 0;
  }
}

std::ostream &operator<<(std::ostream &os, TriIntersectionResult r) {
  switch (r) {
    case TriIntersectionResult::eTriIntersectTrue:
      os << "eTriIntersectTrue";
      break;
    case TriIntersectionResult::eTriIntersectFalse:
      os << "eTriIntersectFalse";
      break;
    case TriIntersectionResult::eTriCoplanar:
      os << "eTriCoplanar";
      break;
  }
  return os;
}

void permuteTriangleVertexAlone(std::array<Vector3, 3> &t, int o1, int o2,
                                int o3) {
  auto permuteTriRight = [](std::array<Vector3, 3> &tt) {
    auto tmp = tt[0];
    tt[0] = tt[1];
    tt[1] = tt[2];
    tt[2] = tmp;
  };
  auto permuteTriLeft = [](std::array<Vector3, 3> &tt) {
    auto tmp = tt[2];
    tt[2] = tt[1];
    tt[1] = tt[0];
    tt[0] = tmp;
  };

  // Permute a, b, c so that a is alone on its side
  if (o1 == o2) {
    // c is alone, permute right so c becomes a
    permuteTriRight(t);
  } else if (o1 == o3) {
    // b is alone, permute so b becomes a
    permuteTriLeft(t);
  } else if (o2 != o3) {
    // In case a, b, c have different orientation, put a on positive side
    if (o2 > 0) {
      permuteTriLeft(t);
    } else if (o3 > 0) {
      permuteTriRight(t);
    }
  }
}

// https://inria.hal.science/file/index/docid/72100/filename/RR-4488.pdf
TriIntersectionResult triangleIntersection3D(std::array<Vector3, 3> t1,
                                             std::array<Vector3, 3> t2,
                                             double epsilon) {
  auto o11 = fourPointDeterminant(t2[0], t2[1], t2[2], t1[0], epsilon);
  auto o12 = fourPointDeterminant(t2[0], t2[1], t2[2], t1[1], epsilon);
  auto o13 = fourPointDeterminant(t2[0], t2[1], t2[2], t1[2], epsilon);

  // All three points of t2 lie within the open half space induced by the plane
  // of t1
  // -> no intersection
  using enum TriIntersectionResult;
  if ((o11 > 0 && o12 > 0 && o13 > 0) || (o11 < 0 && o12 < 0 && o13 < 0)) {
    return eTriIntersectFalse;
  } else if (o11 == 0 && o12 == 0 && o13 == 0) {
    return eTriCoplanar;
  }

  auto o21 = fourPointDeterminant(t1[0], t1[1], t1[2], t2[0], epsilon);
  auto o22 = fourPointDeterminant(t1[0], t1[1], t1[2], t2[1], epsilon);
  auto o23 = fourPointDeterminant(t1[0], t1[1], t1[2], t2[2], epsilon);

  if (o21 == 0 && o23 == 0 && o23 == 0) {
    return eTriCoplanar;
  } else if ((o21 > 0 && o22 > 0 && o23 > 0) ||
             (o21 < 0 && o22 < 0 && o23 < 0)) {
    return eTriIntersectFalse;
  }

  permuteTriangleVertexAlone(t1, o11, o12, o13);
  permuteTriangleVertexAlone(t2, o21, o22, o23);

  o11 = fourPointDeterminant(t2[0], t2[1], t2[2], t1[0], epsilon);
  if (o11 < 0) {
    std::swap(t2[1], t2[2]);
  }

  o21 = fourPointDeterminant(t1[0], t1[1], t1[2], t2[0], epsilon);
  if (o21 < 0) {
    std::swap(t1[1], t1[2]);
  }

  auto o1 = fourPointDeterminant(t1[0], t1[1], t2[0], t2[1], epsilon);
  auto o2 = fourPointDeterminant(t1[0], t1[2], t2[2], t2[0], epsilon);

  if (o1 <= 0 && o2 <= 0) {
    return eTriIntersectTrue;
  } else {
    return eTriIntersectFalse;
  }
}

bool triangleCoplanarIntersect(std::array<Vector3, 3> t1,
                               std::array<Vector3, 3> t2, double epsilon) {
  auto transformTriangle = [](std::array<Vector3, 3> &t) {
    Vector3 v1 = (t[1] - t[0]).normalized();
    Vector3 tmp = (t[2] - t[0]).normalized();
    Vector3 v3 = v1.cross(tmp).normalized();
    Vector3 v2 = v3.cross(v1).normalized();

    SquareMatrix3 m;
    m.col(0) = v1;
    m.col(1) = v2;
    m.col(2) = v3;

    for (auto i = 0ul; i < 3; ++i) {
      t.at(i) = m.inverse() * t.at(i);
    }
  };

  transformTriangle(t1);
  assert(std::abs(t1[0].z() - t1[1].z()) < epsilon);
  assert(std::abs(t1[0].z() - t1[2].z()) < epsilon);

  transformTriangle(t2);
  assert(std::abs(t2[0].z() - t2[1].z()) < epsilon);
  assert(std::abs(t2[0].z() - t2[2].z()) < epsilon);

  // Do the intersection with boost geometry
  namespace bg = boost::geometry;
  typedef bg::model::point<double, 2, bg::cs::cartesian> point_2d;
  typedef bg::model::polygon<point_2d> polygon_2d;

  polygon_2d poly1, poly2;
  for (int i = 0; i < 3; ++i) {
    bg::append(poly1.outer(), point_2d(t1[i].x(), t1[i].y()));
    bg::append(poly2.outer(), point_2d(t2[i].x(), t2[i].y()));
  }

  return bg::intersects(poly1, poly2);
}

bool generalTriangleIntersection(const std::array<Vector3, 3> &t1,
                                 const std::array<Vector3, 3> &t2,
                                 double epsilon) {
  auto res = triangleIntersection3D(t1, t2, epsilon);
  switch (res) {
    case TriIntersectionResult::eTriIntersectTrue:
      return true;
    case TriIntersectionResult::eTriIntersectFalse:
      return false;
    case TriIntersectionResult::eTriCoplanar:
      return triangleCoplanarIntersect(t1, t2, epsilon);
  }
}

bool polyhedronIntersection(const Acts::Polyhedron &p1,
                            const Acts::Polyhedron &p2, double epsilon) {
  if (p1.faces.size() > 1 || p2.faces.size() > 1) {
    throw std::runtime_error(
        "PolyhedronIntersection only works for polyhedra with one face");
  }
  for (const auto &t1 : p1.triangularMesh) {
    for (const auto &t2 : p2.triangularMesh) {
      std::array<Vector3, 3> t1v = {p1.vertices.at(t1.at(0)),
                                    p1.vertices.at(t1.at(1)),
                                    p1.vertices.at(t1.at(2))};
      std::array<Vector3, 3> t2v = {p2.vertices.at(t2.at(0)),
                                    p2.vertices.at(t2.at(1)),
                                    p2.vertices.at(t2.at(2))};

      if (generalTriangleIntersection(t1v, t2v, epsilon)) {
        return true;
      }
    }
  }
  return false;
}

}  // namespace Acts::detail
