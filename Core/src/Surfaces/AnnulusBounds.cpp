// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::AnnulusBounds::AnnulusBounds(double minR, double maxR, double minPhi,
                                   double maxPhi, const Vector2D& moduleOrigin,
                                   double avgPhi)
    : m_rMin(std::min(minR, maxR)),
      m_rMax(std::max(minR, maxR)),
      m_phiMin(std::min(minPhi, maxPhi)),
      m_phiMax(std::max(minPhi, maxPhi)),
      m_moduleOrigin(moduleOrigin),
      m_phiAvg(detail::radian_sym(avgPhi)) {
  m_rotationStripPC = Eigen::Translation<double, 2>(Vector2D(0, -m_phiAvg));
  m_translation = Eigen::Translation<double, 2>(m_moduleOrigin);

  m_shiftXY = m_moduleOrigin * -1;
  m_shiftPC =
      Vector2D(VectorHelpers::perp(m_shiftXY), VectorHelpers::phi(m_shiftXY));

  // we need the corner points of the module to do the inside
  // checking, calculate them here once, they don't change

  // find inner outer radius at edges in STRIP PC
  auto circIx = [](double O_x, double O_y, double r, double phi) -> Vector2D {
    //                      _____________________________________________
    //                     /      2  2                    2    2  2    2
    //     O_x + O_y*m - \/  - O_x *m  + 2*O_x*O_y*m - O_y  + m *r  + r
    // x = --------------------------------------------------------------
    //                                  2
    //                                 m  + 1
    //
    // y = m*x
    //
    double m = std::tan(phi);
    Vector2D dir(std::cos(phi), std::sin(phi));
    double x1 = (O_x + O_y * m -
                 std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) +
                           2 * O_x * O_y * m - std::pow(O_y, 2) +
                           std::pow(m, 2) * std::pow(r, 2) + std::pow(r, 2))) /
                (std::pow(m, 2) + 1);
    double x2 = (O_x + O_y * m +
                 std::sqrt(-std::pow(O_x, 2) * std::pow(m, 2) +
                           2 * O_x * O_y * m - std::pow(O_y, 2) +
                           std::pow(m, 2) * std::pow(r, 2) + std::pow(r, 2))) /
                (std::pow(m, 2) + 1);

    Vector2D v1(x1, m * x1);
    if (v1.dot(dir) > 0)
      return v1;
    return {x2, m * x2};
  };

  // calculate corners in STRIP XY, keep them we need them for minDistance()
  m_outLeftStripXY =
      circIx(m_moduleOrigin[eLOC_X], m_moduleOrigin[eLOC_Y], m_rMax, m_phiMax);
  m_inLeftStripXY =
      circIx(m_moduleOrigin[eLOC_X], m_moduleOrigin[eLOC_Y], m_rMin, m_phiMax);
  m_outRightStripXY =
      circIx(m_moduleOrigin[eLOC_X], m_moduleOrigin[eLOC_Y], m_rMax, m_phiMin);
  m_inRightStripXY =
      circIx(m_moduleOrigin[eLOC_X], m_moduleOrigin[eLOC_Y], m_rMin, m_phiMin);

  m_outLeftStripPC = {m_outLeftStripXY.norm(),
                      VectorHelpers::phi(m_outLeftStripXY)};
  m_inLeftStripPC = {m_inLeftStripXY.norm(),
                     VectorHelpers::phi(m_inLeftStripXY)};
  m_outRightStripPC = {m_outRightStripXY.norm(),
                       VectorHelpers::phi(m_outRightStripXY)};
  m_inRightStripPC = {m_inRightStripXY.norm(),
                      VectorHelpers::phi(m_inRightStripXY)};

  m_outLeftModulePC = stripXYToModulePC(m_outLeftStripXY);
  m_inLeftModulePC = stripXYToModulePC(m_inLeftStripXY);
  m_outRightModulePC = stripXYToModulePC(m_outRightStripXY);
  m_inRightModulePC = stripXYToModulePC(m_inRightStripXY);
}

std::vector<TDD_real_t> Acts::AnnulusBounds::valueStore() const {
  std::vector<TDD_real_t> values(AnnulusBounds::bv_length);
  values[AnnulusBounds::bv_minR] = rMin();
  values[AnnulusBounds::bv_maxR] = rMax();
  values[AnnulusBounds::bv_phiMin] = phiMin();
  values[AnnulusBounds::bv_phiMax] = phiMax();
  values[AnnulusBounds::bv_phiAvg] = 0.5 * (phiMin() + phiMax());
  values[AnnulusBounds::bv_originX] = m_moduleOrigin.x();
  values[AnnulusBounds::bv_originY] = m_moduleOrigin.y();
  return values;
}

std::vector<Acts::Vector2D> Acts::AnnulusBounds::corners() const {
  auto rot = m_rotationStripPC.inverse();

  return {rot * m_outRightStripPC, rot * m_outLeftStripPC,
          rot * m_inLeftStripPC, rot * m_inRightStripPC};
}

std::vector<Acts::Vector2D> Acts::AnnulusBounds::vertices(
    unsigned int lseg) const {
  // List of vertices counter-clockwise starting with left inner
  std::vector<Acts::Vector2D> rvertices;

  double phiMinInner = VectorHelpers::phi(m_inLeftStripXY);
  double phiMaxInner = VectorHelpers::phi(m_inRightStripXY);
  double phiMinOuter = VectorHelpers::phi(m_outRightStripXY);
  double phiMaxOuter = VectorHelpers::phi(m_outLeftStripXY);

  std::vector<double> phisInner =
      detail::VerticesHelper::phiSegments(phiMinInner, phiMaxInner);
  std::vector<double> phisOuter =
      detail::VerticesHelper::phiSegments(phiMinOuter, phiMaxOuter);

  // Inner bow from phi_min -> phi_max
  for (unsigned int iseg = 0; iseg < phisInner.size() - 1; ++iseg) {
    int addon = (iseg == phisInner.size() - 2) ? 1 : 0;
    detail::VerticesHelper::createSegment<Vector2D, Eigen::Affine2d>(
        rvertices, {rMin(), rMin()}, phisInner[iseg], phisInner[iseg + 1], lseg,
        addon);
  }
  // Upper bow from phi_min -> phi_max
  for (unsigned int iseg = 0; iseg < phisOuter.size() - 1; ++iseg) {
    int addon = (iseg == phisOuter.size() - 2) ? 1 : 0;
    detail::VerticesHelper::createSegment<Vector2D, Eigen::Affine2d>(
        rvertices, {rMax(), rMax()}, phisOuter[iseg], phisOuter[iseg + 1], lseg,
        addon);
  }

  return rvertices;
}

bool Acts::AnnulusBounds::inside(const Vector2D& lposition, double tolR,
                                 double tolPhi) const {
  // locpo is PC in STRIP SYSTEM
  // need to perform internal rotation induced by m_phiAvg
  Vector2D locpo_rotated = m_rotationStripPC * lposition;
  double phiLoc = locpo_rotated[eLOC_PHI];
  double rLoc = locpo_rotated[eLOC_R];

  if (phiLoc < (m_phiMin - tolPhi) || phiLoc > (m_phiMax + tolPhi)) {
    return false;
  }

  // calculate R in MODULE SYSTEM to evaluate R-bounds
  if (tolR == 0.) {
    // don't need R, can use R^2
    double r_mod2 =
        m_shiftPC[eLOC_R] * m_shiftPC[eLOC_R] + rLoc * rLoc +
        2 * m_shiftPC[eLOC_R] * rLoc * cos(phiLoc - m_shiftPC[eLOC_PHI]);

    if (r_mod2 < m_rMin * m_rMin || r_mod2 > m_rMax * m_rMax) {
      return false;
    }
  } else {
    // use R
    double r_mod =
        sqrt(m_shiftPC[eLOC_R] * m_shiftPC[eLOC_R] + rLoc * rLoc +
             2 * m_shiftPC[eLOC_R] * rLoc * cos(phiLoc - m_shiftPC[eLOC_PHI]));

    if (r_mod < (m_rMin - tolR) || r_mod > (m_rMax + tolR)) {
      return false;
    }
  }
  return true;
}

bool Acts::AnnulusBounds::inside(const Vector2D& lposition,
                                 const BoundaryCheck& bcheck) const {
  // locpo is PC in STRIP SYSTEM
  if (bcheck.type() == BoundaryCheck::Type::eAbsolute) {
    return inside(lposition, bcheck.tolerance()[eLOC_R],
                  bcheck.tolerance()[eLOC_PHI]);
  } else {
    // first check if inside. We don't need to look into the covariance if
    // inside
    if (inside(lposition, 0., 0.)) {
      return true;
    }

    // we need to rotated the locpo
    Vector2D locpo_rotated = m_rotationStripPC * lposition;

    // covariance is given in STRIP SYSTEM in PC
    // we need to convert the covariance to the MODULE SYSTEM in PC
    // via jacobian.
    // The following transforms into STRIP XY, does the shift into MODULE XY,
    // and then transforms into MODULE PC
    double dphi = m_phiAvg;
    double phi_strip = locpo_rotated[eLOC_PHI];
    double r_strip = locpo_rotated[eLOC_R];
    double O_x = m_shiftXY[eLOC_X];
    double O_y = m_shiftXY[eLOC_Y];

    // For a transformation from cartesian into polar coordinates
    //
    //              [         _________      ]
    //              [        /  2    2       ]
    //              [      \/  x  + y        ]
    //     [ r' ]   [                        ]
    // v = [    ] = [      /       y        \]
    //     [phi']   [2*atan|----------------|]
    //              [      |       _________|]
    //              [      |      /  2    2 |]
    //              [      \x + \/  x  + y  /]
    //
    // Where x, y are polar coordinates that can be rotated by dPhi
    //
    // [x]   [O_x + r*cos(dPhi - phi)]
    // [ ] = [                       ]
    // [y]   [O_y - r*sin(dPhi - phi)]
    //
    // The general jacobian is:
    //
    //        [d        d      ]
    //        [--(f_x)  --(f_x)]
    //        [dx       dy     ]
    // Jgen = [                ]
    //        [d        d      ]
    //        [--(f_y)  --(f_y)]
    //        [dx       dy     ]
    //
    // which means in this case:
    //
    //     [     d                   d           ]
    //     [ ----------(rMod)    ---------(rMod) ]
    //     [ dr_{strip}          dphiStrip       ]
    // J = [                                     ]
    //     [    d                   d            ]
    //     [----------(phiMod)  ---------(phiMod)]
    //     [dr_{strip}          dphiStrip        ]
    //
    // Performing the derivative one gets:
    //
    //     [B*O_x + C*O_y + rStrip  rStrip*(B*O_y + O_x*sin(dPhi - phiStrip))]
    //     [----------------------  -----------------------------------------]
    //     [          ___                               ___                  ]
    //     [        \/ A                              \/ A                   ]
    // J = [                                                                 ]
    //     [  -(B*O_y - C*O_x)           rStrip*(B*O_x + C*O_y + rStrip)     ]
    //     [  -----------------          -------------------------------     ]
    //     [          A                                 A                    ]
    //
    // where
    //        2                                          2 2
    // A = O_x  + 2*O_x*rStrip*cos(dPhi - phiStrip) + O_y  -
    // 2*O_y*rStrip*sin(dPhi - phiStrip) + rStrip B = cos(dPhi - phiStrip) C =
    // -sin(dPhi - phiStrip)

    double cosDPhiPhiStrip = std::cos(dphi - phi_strip);
    double sinDPhiPhiStrip = std::sin(dphi - phi_strip);

    double A = O_x * O_x + 2 * O_x * r_strip * cosDPhiPhiStrip + O_y * O_y -
               2 * O_y * r_strip * sinDPhiPhiStrip + r_strip * r_strip;
    double sqrtA = std::sqrt(A);

    double B = cosDPhiPhiStrip;
    double C = -sinDPhiPhiStrip;
    Eigen::Matrix<double, 2, 2> jacobianStripPCToModulePC;
    jacobianStripPCToModulePC(0, 0) = (B * O_x + C * O_y + r_strip) / sqrtA;
    jacobianStripPCToModulePC(0, 1) =
        r_strip * (B * O_y + O_x * sinDPhiPhiStrip) / sqrtA;
    jacobianStripPCToModulePC(1, 0) = -(B * O_y - C * O_x) / A;
    jacobianStripPCToModulePC(1, 1) =
        r_strip * (B * O_x + C * O_y + r_strip) / A;

    // covariance is given in STRIP PC
    auto covStripPC = bcheck.covariance();
    // calculate covariance in MODULE PC using jacobian from above
    auto covModulePC = jacobianStripPCToModulePC * covStripPC *
                       jacobianStripPCToModulePC.transpose();

    // Mahalanobis distance uses inverse covariance as weights
    auto weightStripPC = covStripPC.inverse();
    auto weightModulePC = covModulePC.inverse();

    double minDist = std::numeric_limits<double>::max();

    Vector2D currentClosest;
    double currentDist;

    // do projection in STRIP PC

    // first: STRIP system. locpo is in STRIP PC already
    currentClosest = closestOnSegment(m_inLeftStripPC, m_outLeftStripPC,
                                      locpo_rotated, weightStripPC);
    currentDist = squaredNorm(locpo_rotated - currentClosest, weightStripPC);
    minDist = currentDist;

    currentClosest = closestOnSegment(m_inRightStripPC, m_outRightStripPC,
                                      locpo_rotated, weightStripPC);
    currentDist = squaredNorm(locpo_rotated - currentClosest, weightStripPC);
    if (currentDist < minDist) {
      minDist = currentDist;
    }

    // now: MODULE system. Need to transform locpo to MODULE PC
    //  transform is STRIP PC -> STRIP XY -> MODULE XY -> MODULE PC
    Vector2D locpoStripXY(
        locpo_rotated[eLOC_R] * std::cos(locpo_rotated[eLOC_PHI]),
        locpo_rotated[eLOC_R] * std::sin(locpo_rotated[eLOC_PHI]));
    Vector2D locpoModulePC = stripXYToModulePC(locpoStripXY);

    // now check edges in MODULE PC (inner and outer circle)
    // assuming Mahalanobis distances are of same unit if covariance
    // is correctly transformed
    currentClosest = closestOnSegment(m_inLeftModulePC, m_inRightModulePC,
                                      locpoModulePC, weightModulePC);
    currentDist = squaredNorm(locpoModulePC - currentClosest, weightModulePC);
    if (currentDist < minDist) {
      minDist = currentDist;
    }

    currentClosest = closestOnSegment(m_outLeftModulePC, m_outRightModulePC,
                                      locpoModulePC, weightModulePC);
    currentDist = squaredNorm(locpoModulePC - currentClosest, weightModulePC);
    if (currentDist < minDist) {
      minDist = currentDist;
    }

    // compare resulting Mahalanobis distance to configured
    // "number of sigmas"
    // we square it b/c we never took the square root of the distance
    return minDist < bcheck.tolerance()[0] * bcheck.tolerance()[0];
  }
}

double Acts::AnnulusBounds::distanceToBoundary(
    const Vector2D& lposition) const {
  // find the closest point on all edges, calculate distance
  // return smallest one
  // closest distance is cartesian, we want the result in mm.

  Vector2D locpo_rotated = m_rotationStripPC * lposition;

  // locpo is given in STRIP PC, we need it in STRIP XY and possibly MODULE XY
  double rStrip = locpo_rotated[eLOC_R];
  double phiStrip = locpo_rotated[eLOC_PHI];
  Vector2D locpoStripXY(rStrip * std::cos(phiStrip),
                        rStrip * std::sin(phiStrip));
  Vector2D locpoModuleXY = locpoStripXY + m_shiftXY;
  double rMod = locpoModuleXY.norm();
  double phiMod = VectorHelpers::phi(locpoModuleXY);

  Vector2D closestStripPC;
  double minDist = std::numeric_limits<double>::max();
  ;
  double curDist;

  // for rmin
  if (m_inRightModulePC[eLOC_PHI] <= phiMod &&
      phiMod < m_inLeftModulePC[eLOC_PHI]) {
    // is inside phi bounds, to comparison to rmin and r max
    // r min
    curDist = std::abs(m_rMin - rMod);
    if (curDist < minDist) {
      minDist = curDist;
    }
  } else {
    // is outside phi bounds, closest can only be the edge points here

    // in left
    curDist = (m_inLeftStripXY - locpoStripXY).norm();
    if (curDist < minDist) {
      minDist = curDist;
    }

    // in right
    curDist = (m_inRightStripXY - locpoStripXY).norm();
    if (curDist < minDist) {
      minDist = curDist;
    }
  }

  if (m_phiMin <= phiStrip && phiStrip < m_phiMax) {
    // r max
    curDist = std::abs(m_rMax - rMod);
    if (curDist < minDist) {
      minDist = curDist;
    }
  } else {
    // out left
    curDist = (m_outLeftStripXY - locpoStripXY).norm();
    if (curDist < minDist) {
      minDist = curDist;
    }

    // out right
    curDist = (m_outRightStripXY - locpoStripXY).norm();
    if (curDist < minDist) {
      minDist = curDist;
    }
  }

  ActsMatrixD<2, 2> weight = ActsMatrixD<2, 2>::Identity();

  // phi left
  Vector2D closestLeft =
      closestOnSegment(m_inLeftStripXY, m_outLeftStripXY, locpoStripXY, weight);
  curDist = (closestLeft - locpoStripXY).norm();
  if (curDist < minDist) {
    minDist = curDist;
  }

  // phi right
  Vector2D closestRight = closestOnSegment(m_inRightStripXY, m_outRightStripXY,
                                           locpoStripXY, weight);
  curDist = (closestRight - locpoStripXY).norm();
  if (curDist < minDist) {
    minDist = curDist;
  }

  return minDist;
}

Acts::Vector2D Acts::AnnulusBounds::stripXYToModulePC(
    const Vector2D& vStripXY) const {
  Vector2D vecModuleXY = vStripXY + m_shiftXY;
  return {vecModuleXY.norm(), VectorHelpers::phi(vecModuleXY)};
}

Acts::Vector2D Acts::AnnulusBounds::closestOnSegment(
    const Vector2D& a, const Vector2D& b, const Vector2D& p,
    const Eigen::Matrix<double, 2, 2>& weight) const {
  // connecting vector
  auto n = b - a;
  // squared norm of line
  auto f = (n.transpose() * weight * n).value();
  // weighted scalar product of line to point and segment line
  auto u = ((p - a).transpose() * weight * n).value() / f;
  // clamp to [0, 1], convert to point
  return std::min(std::max(u, 0.0), 1.0) * n + a;
}

double Acts::AnnulusBounds::squaredNorm(
    const Vector2D& v, const Eigen::Matrix<double, 2, 2>& weight) const {
  return (v.transpose() * weight * v).value();
}

Acts::Vector2D Acts::AnnulusBounds::moduleOrigin() const {
  return Eigen::Rotation2D<double>(m_phiAvg) * m_moduleOrigin;
}

// Ostream operator overload
std::ostream& Acts::AnnulusBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::AnnulusBounds:  (innerRadius, outerRadius, minPhi, maxPhi) = ";
  sl << "(" << rMin() << ", " << rMax() << ", " << phiMin() << ", " << phiMax()
     << ")" << '\n';
  sl << " - shift xy = " << m_shiftXY.x() << ", " << m_shiftXY.y() << '\n';
  ;
  sl << " - shift pc = " << m_shiftPC.x() << ", " << m_shiftPC.y() << '\n';
  ;
  sl << std::setprecision(-1);
  return sl;
}
