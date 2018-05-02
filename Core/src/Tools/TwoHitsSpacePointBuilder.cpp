// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <limits>
#include <stdexcept>
#include "ACTS/Tools/TwoHitsSpacePointBuilder.hpp"

///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///

Acts::TwoHitsSpacePointBuilder::TwoHitsSpacePointBuilder(const Config& cfg)
  : m_cfg(cfg)
{
}

Acts::Vector2D
Acts::TwoHitsSpacePointBuilder::localCoords(
    const Acts::PlanarModuleCluster& hit) const
{
  // Local position information
  auto           par = hit.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

Acts::Vector3D
Acts::TwoHitsSpacePointBuilder::globalCoords(
    const Acts::PlanarModuleCluster& hit) const
{
  // Receive corresponding surface
  auto& clusterSurface = hit.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(localCoords(hit), mom, pos);

  return pos;
}

double
Acts::TwoHitsSpacePointBuilder::differenceOfHits(
    const Acts::PlanarModuleCluster& hit1,
    const Acts::PlanarModuleCluster& hit2) const
{
  // Calculate the global position of the hits
  Acts::Vector3D pos1 = globalCoords(hit1);
  Acts::Vector3D pos2 = globalCoords(hit2);

  // Check if measurements are close enough to each other
  if ((pos1 - pos2).norm() > m_cfg.diffDist) return -1.;

  // Calculate the angles of the hits
  double phi1, theta1, phi2, theta2;
  phi1   = (pos1 - m_cfg.vertex).phi();
  theta1 = (pos1 - m_cfg.vertex).theta();
  phi2   = (pos2 - m_cfg.vertex).phi();
  theta2 = (pos2 - m_cfg.vertex).theta();

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > m_cfg.diffTheta2) return -1.;

  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > m_cfg.diffPhi2) return -1.;

  // Return the squared distance between both hits
  return diffTheta2 + diffPhi2;
}

void
Acts::TwoHitsSpacePointBuilder::addHits(std::vector<std::vector<Acts::PlanarModuleCluster const*>>& hits)
{
  // Return if more/less surfaces given than 2
  if(hits.size() != 2) return;
  // Return if no hits are given
  if (hits[0].empty() && hits[1].empty()) return;
  
  // TODO: only the closest differences get selected -> some points are not
  // taken into account
  // Declare helper variables
  double                                       currentDiff;
  Acts::SpacePoint tmpSpacePoint;
  double                                       diffMin;
  unsigned int                                 hitMin;

  // Walk through all hits on both surfaces
  for (unsigned int iHits0 = 0; iHits0 < hits[0].size(); iHits0++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of hits
    hitMin = hits[1].size();
    for (unsigned int iHits1 = 0; iHits1 < hits[1].size(); iHits1++) {
      // Calculate the distances between the hits
      currentDiff = differenceOfHits(*(hits[0][iHits0]), *(hits[1][iHits1]));
      // Store the closest hits (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        hitMin  = iHits1;
      }
    }
    // Store the best (=closest) result
    if (hitMin < hits[1].size()) {
	  tmpSpacePoint.hitModule.resize(2);
	  tmpSpacePoint.hitModule.push_back(hits[0][iHits0]);
      tmpSpacePoint.hitModule.push_back(hits[1][hitMin]);
      m_allCombSpacePoints.push_back(tmpSpacePoint);
    }
  }
}

void
Acts::TwoHitsSpacePointBuilder::addSpacePoint(
    Acts::SpacePoint& sPoint)
{
  if (sPoint.hitModule.size() == 2 && sPoint.hitModule[0] && sPoint.hitModule[1])
    m_allCombSpacePoints.push_back(sPoint);
}

std::pair<Acts::Vector3D, Acts::Vector3D>
Acts::TwoHitsSpacePointBuilder::endsOfStrip(
    const Acts::PlanarModuleCluster& hit) const
{
  // Calculate the local coordinates of the hit
  const Acts::Vector2D local = localCoords(hit);

  // Receive the binning
  auto& sur     = hit.referenceSurface();
  auto  segment = dynamic_cast<const Acts::CartesianSegmentation*>(
      &(sur.associatedDetectorElement()->digitizationModule()->segmentation()));
  auto& binData     = segment->binUtility().binningData();
  auto& boundariesX = binData[0].boundaries();
  auto& boundariesY = binData[1].boundaries();

  // Search the x-/y-bin hit
  size_t binX = binData[0].searchLocal(local);
  size_t binY = binData[1].searchLocal(local);

  Acts::Vector2D topLocal, bottomLocal;

  if (boundariesX[binX + 1] - boundariesX[binX]
      < boundariesY[binY + 1] - boundariesY[binY]) {
    // Set the top and bottom end of the strip in local coordinates
    topLocal = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                boundariesY[binY + 1]};
    bottomLocal
        = {(boundariesX[binX] + boundariesX[binX + 1]) / 2, boundariesY[binY]};
  } else {
    // Set the top and bottom end of the strip in local coordinates
    topLocal
        = {boundariesX[binX], (boundariesY[binY] + boundariesY[binY + 1]) / 2};
    bottomLocal = {boundariesX[binX + 1],
                   (boundariesY[binY] + boundariesY[binY + 1]) / 2};
  }

  // Calculate the global coordinates of the top and bottom end of the strip
  Acts::Vector3D topGlobal, bottomGlobal, mom;  // mom is a dummy variable
  sur.localToGlobal(topLocal, mom, topGlobal);
  sur.localToGlobal(bottomLocal, mom, bottomGlobal);

  // Return the top and bottom end of the strip in global coordinates
  return std::make_pair(topGlobal, bottomGlobal);
}

double
Acts::TwoHitsSpacePointBuilder::calcPerpProj(const Acts::Vector3D& a,
                                             const Acts::Vector3D& c,
                                             const Acts::Vector3D& q,
                                             const Acts::Vector3D& r) const
{
  /// This approach assumes that no vertex is available. This option aims to
  /// approximate the space points from cosmic data.
  /// The underlying assumption is that the best point is given by the closest
  /// distance between both lines describing the SDEs.
  /// The point x on the first SDE is parametrized as a + lambda0 * q with the
  /// top end a of the strip and the vector q = a - b(ottom end of the strip).
  /// An analogous parametrization is performed of the second SDE with y = c +
  /// lambda1 * r.
  /// x get resolved by resolving lambda0 from the condition that |x-y| is the
  /// shortest distance between two skew lines.
  Acts::Vector3D ac    = c - a;
  double         qr    = q.dot(r);
  double         denom = q.dot(q) - qr * qr;

  // Check for numerical stability
  if (fabs(denom) > 10e-7)
    // Return lambda0
    return (ac.dot(r) * qr - ac.dot(q) * r.dot(r)) / denom;
  // lambda0 is in the interval [-1,0]. This return serves as error check.
  return 1.;
}

bool
Acts::TwoHitsSpacePointBuilder::recoverSpacePoint(
    Acts::TwoHitsSpacePointBuilder::SpacePointParameters& spaPoPa) const
{
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (m_cfg.stripLengthGapTolerance <= 0.) return false;
  spaPoPa.qmag = spaPoPa.q.mag();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spaPoPa.limitExtended
      = spaPoPa.limit + m_cfg.stripLengthGapTolerance / spaPoPa.qmag;
  // Check if m is just slightly outside
  if (fabs(spaPoPa.m) > spaPoPa.limitExtended) return false;
  // Calculate n if not performed previously
  if (spaPoPa.n == 0.)
    spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs);
  // Check if n is just slightly outside
  if (fabs(spaPoPa.n) > spaPoPa.limitExtended) return false;

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the space point can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other measurement, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory. The
  /// would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  double secOnFirstScale
      = spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.qmag * spaPoPa.qmag);
  // Check if both overshoots are in the same direction
  if (spaPoPa.m > 1. && spaPoPa.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spaPoPa.m - 1.;
    double nOvershoot
        = (spaPoPa.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot
        = (mOvershoot > nOvershoot) ? mOvershoot : nOvershoot;
    // Move m and n towards 0
    spaPoPa.m -= biggerOvershoot;
    spaPoPa.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // Check if both overshoots are in the same direction
  if (spaPoPa.m < -1. && spaPoPa.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spaPoPa.m + 1.);
    double nOvershoot
        = -(spaPoPa.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot
        = (mOvershoot > nOvershoot) ? mOvershoot : nOvershoot;
    // Move m and n towards 0
    spaPoPa.m += biggerOvershoot;
    spaPoPa.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // No solution could be found
  return false;
}

void
Acts::TwoHitsSpacePointBuilder::calculateSpacePoints()
{
  /// Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  Acts::TwoHitsSpacePointBuilder::SpacePointParameters spaPoPa;

  // Walk over every found candidate pair
  for (auto& hits : m_allCombSpacePoints) {
	assert(hits.hitModule.size() == 2);
    // If the space point is already calculated this can be skipped
    if (hits.spacePoint != Acts::Vector3D::Zero(3)) continue;

    // Calculate the ends of the SDEs
    const auto& ends1 = endsOfStrip(*(hits.hitModule[0]));
    const auto& ends2 = endsOfStrip(*(hits.hitModule[1]));

    /// The following algorithm is meant for finding the position on the first
    /// strip if there is a corresponding hit on the second strip. The
    /// resulting point is a point x on the first surfaces. This point is
    /// along a line between the points a (top end of the strip)
    /// and b (bottom end of the strip). The location can be parametrized as
    /// 	2 * x = (1 + m) a + (1 - m) b
    /// as function of the scalar m. m is a parameter in the interval
    /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
    /// from the vertex to the hit on the second strip y is needed to be a
    /// multiple k of the vector from vertex to the hit on the first strip x.
    /// As a consequence of this demand y = k * x needs to be on the
    /// connecting line between the top (c) and bottom (d) end of
    /// the second strip. If both hits correspond to each other, the condition
    /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
    /// needs to be fulfilled. Inserting the first equation into this
    /// equation leads to the condition for m as given in the following
    /// algorithm and therefore to the calculation of x.
    /// The same calculation can be repeated for y. Its corresponding
    /// parameter will be named n.

    spaPoPa.reset();
    spaPoPa.q = ends1.first - ends1.second;
    spaPoPa.r = ends2.first - ends2.second;

    // Fast skipping if a perpendicular projection should be used
    double resultPerpProj;
    if (m_cfg.usePerpProj
        && (resultPerpProj
            = calcPerpProj(ends1.first, ends2.first, spaPoPa.q, spaPoPa.r)
                <= 0.)) {
      hits.spacePoint = ends1.first + resultPerpProj * spaPoPa.q;
      continue;
    }

    spaPoPa.s  = ends1.first + ends1.second - 2 * m_cfg.vertex;
    spaPoPa.t  = ends2.first + ends2.second - 2 * m_cfg.vertex;
    spaPoPa.qs = spaPoPa.q.cross(spaPoPa.s);
    spaPoPa.rt = spaPoPa.r.cross(spaPoPa.t);
    spaPoPa.m  = -spaPoPa.s.dot(spaPoPa.rt) / spaPoPa.q.dot(spaPoPa.rt);

    // Set the limit for the parameter
    if (spaPoPa.limit == 1. && m_cfg.stripLengthTolerance != 0.)
      spaPoPa.limit = 1. + m_cfg.stripLengthTolerance;

    // Check if m and n can be resolved in the interval (-1, 1)
    if (fabs(spaPoPa.m) <= spaPoPa.limit
        && fabs(spaPoPa.n
                = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs))
            <= spaPoPa.limit)
      // Store the space point
      hits.spacePoint
          = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
    else
        /// If this point is reached then it was not possible to resolve both
        /// points such that they are on their SDEs
        /// The following code treats a possible recovery of points resolved
        /// slightly outside of the SDE.
        /// @note This procedure is an indirect variation of the vertex
        /// position.
        // Check if a recovery the point(s) and store them if successful
        if (recoverSpacePoint(spaPoPa))
      hits.spacePoint
          = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
  }
}

const std::vector<Acts::SpacePoint>&
Acts::TwoHitsSpacePointBuilder::spacePoints()
{
  return m_allCombSpacePoints;
}
