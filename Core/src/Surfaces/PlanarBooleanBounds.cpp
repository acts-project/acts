// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain left at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PlanarBooleanBounds.hpp"
#include <iostream>

Acts::PlanarBooleanBounds::PlanarBooleanBounds(
    std::shared_ptr<const PlanarBounds> left, const Vector2D& leftShift,
    std::shared_ptr<const PlanarBounds> right, const Vector2D& rightShift,
    BooleanOperation bOperation) noexcept(false)
    : m_left(left),
      m_leftShift(leftShift),
      m_right(right),
      m_rightShift(rightShift),
      m_bOperation(bOperation),
      m_boundingBox(0., 0.) {
  auto leftBB = m_left->boundingBox();
  auto rightBB = m_right->boundingBox();

  auto lVertices = m_left->vertices(1);
  auto rVertices = m_right->vertices(1);

  m_overlap = false;
  for (const auto& rv : rVertices) {
    m_overlap = m_left->inside(rv + m_rightShift - m_leftShift, true);
    if (m_overlap) {
      break;
    }
  }

  if (m_bOperation != eUnion and not m_overlap) {
    throw std::invalid_argument(
        "PlanarBooleanBounds: "
        "shapes must overlap for boolean operation other than union.");
  }

  auto ivertices = vertices(1);
  double xmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::lowest();
  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::lowest();
  for (const auto v : ivertices) {
    xmin = std::min(xmin, v.x());
    xmax = std::max(xmax, v.x());
    ymin = std::min(ymin, v.y());
    ymax = std::max(ymax, v.y());
  }
  m_boundingBox = RectangleBounds({xmin, ymin}, {xmax, ymax});
}

std::vector<double> Acts::PlanarBooleanBounds::values() const {
  std::vector<double> rvalues = m_left->values();
  std::vector<double> ovalues = m_right->values();
  rvalues.insert(rvalues.end(), ovalues.begin(), ovalues.end());
  return rvalues;
}

bool Acts::PlanarBooleanBounds::inside(const Vector2D& lposition,
                                       const BoundaryCheck& bcheck) const {
  if (m_bOperation == eUnion) {
    return (m_left->inside(lposition - m_leftShift, bcheck) or
            m_right->inside(lposition - m_rightShift, bcheck));
  }
  if (m_bOperation == eIntersection) {
    return (m_left->inside(lposition - m_leftShift, bcheck) and
            m_right->inside(lposition - m_rightShift, bcheck));
  }
  if (m_bOperation == eNot) {
    return (m_left->inside(lposition - m_leftShift, bcheck) and
            not m_right->inside(lposition - m_rightShift, bcheck));
  }
  return (m_left->inside(lposition - m_leftShift, bcheck)) ^
         (m_right->inside(lposition - m_rightShift, bcheck));
}

double Acts::PlanarBooleanBounds::distanceToBoundary(
    const Vector2D& lposition) const {
  double dL = m_left->distanceToBoundary(lposition - m_leftShift);
  double dR = m_right->distanceToBoundary(lposition - m_rightShift);

  if (m_bOperation == eUnion) {
    if (dL * dR > 0) {
      return std::copysign(std::min(std::abs(dL), std::abs(dR)), dL);
    }
    return dL > 0 ? dL : dR;
  }

  return 0.;
}

std::vector<Acts::Vector2D> Acts::PlanarBooleanBounds::vertices(
    unsigned int lseg) const {
  auto lVertices = m_left->vertices(lseg);
  std::for_each(lVertices.begin(), lVertices.end(),
                [&](Vector2D& vtx) { vtx += m_leftShift; });

  auto rVertices = m_right->vertices(lseg);
  std::for_each(rVertices.begin(), rVertices.end(),
                [&](Vector2D& vtx) { vtx += m_rightShift; });

  // They don't overlap, take the two sets
  if (not m_overlap) {
    lVertices.insert(lVertices.begin(), rVertices.begin(), rVertices.end());
    return lVertices;
  }

  // Swtich through the cases
  std::vector<Vector2D> cVertices;
  switch (m_bOperation) {
    case eUnion: {
      cVertices = lVertices;
      for (const auto& rv : rVertices) {
        if (not m_left->inside(rv - m_rightShift + m_leftShift, true)) {
          lVertices.push_back(rv);
        }
      }
    } break;

    case eIntersection: {
      for (const auto& lv : lVertices) {
        if (m_right->inside(lv + m_rightShift - m_leftShift, true)) {
          cVertices.push_back(lv);
        }
      }

      for (const auto& rv : rVertices) {
        if (m_left->inside(rv - m_rightShift + m_leftShift, true)) {
          cVertices.push_back(rv);
        }
      }
    } break;

    case eNot: {
      for (const auto& lv : lVertices) {
        if (not m_right->inside(lv + m_rightShift - m_leftShift, true)) {
          cVertices.push_back(lv);
        }
      }

    } break;

    case eXor: {
      for (const auto& lv : lVertices) {
        if (not m_right->inside(lv + m_rightShift - m_leftShift, true)) {
          cVertices.push_back(lv);
        }
      }

      for (const auto& rv : rVertices) {
        if (not m_left->inside(rv - m_rightShift + m_leftShift, true)) {
          cVertices.push_back(rv);
        }
      }
    } break;
  }
  return cVertices;
}

std::ostream& Acts::PlanarBooleanBounds::toStream(std::ostream& os) const {
  os << "Acts::PlanarBooleanBounds - "
     << ((m_bOperation == eUnion) ? "uninon of " : "intersection of ") << '\n';
  os << "- left side  : ";
  m_left->toStream(os);
  os << "- right side : ";
  m_right->toStream(os);
  return os;
}
