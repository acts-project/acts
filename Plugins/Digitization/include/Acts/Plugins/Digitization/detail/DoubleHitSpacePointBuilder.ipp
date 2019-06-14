// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <limits>
#include "Acts/Plugins/Digitization/DoubleHitSpacePointBuilder.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace {
/// @brief Storage container for variables related to the calculation of space
/// points
struct SpacePointParameters {
  /// Vector pointing from bottom to top end of first SDE
  Vector3D q;
  /// Vector pointing from bottom to top end of second SDE
  Vector3D r;
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3D s;
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3D t;
  /// Cross product between SpacePointParameters::q and
  /// SpacePointParameters::s
  Vector3D qs;
  /// Cross product between SpacePointParameters::r and
  /// SpacePointParameters::t
  Vector3D rt;
  /// Magnitude of SpacePointParameters::q
  double qmag = 0.;
  /// Parameter that determines the hit position on the first SDE
  double m = 0.;
  /// Parameter that determines the hit position on the second SDE
  double n = 0.;
  /// Regular limit of the absolut values of SpacePointParameters::m and
  /// SpacePointParameters::n
  double limit = 1.;
  /// Limit of SpacePointParameters::m and SpacePointParameters::n in case of
  /// variable vertex
  double limitExtended = 0.;
};

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two clusters
///
/// @param [in] pos1 position of the first cluster
/// @param [in] pos2 position the second cluster
/// @param [in] maxDistance Maximum distance between two clusters
/// @param [in] maxAngleTheta2 Maximum squared theta angle between two clusters
/// @param [in] maxAnglePhi2 Maximum squared phi angle between two clusters
///
/// @return The squared sum within configuration parameters, otherwise -1
double differenceOfClustersChecked(const Vector3D& pos1, const Vector3D& pos2,
                                   const Vector3D& posVertex,
                                   const double maxDistance,
                                   const double maxAngleTheta2,
                                   const double maxAnglePhi2) {
  // Check if measurements are close enough to each other
  if ((pos1 - pos2).norm() > maxDistance) {
    return -1.;
  }

  // Calculate the angles of the vectors
  double phi1, theta1, phi2, theta2;
  phi1 = VectorHelpers::phi(pos1 - posVertex);
  theta1 = VectorHelpers::theta(pos1 - posVertex);
  phi2 = VectorHelpers::phi(pos2 - posVertex);
  theta2 = VectorHelpers::theta(pos2 - posVertex);

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > maxAngleTheta2) {
    return -1.;
  }

  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > maxAnglePhi2) {
    return -1.;
  }

  // Return the squared distance between both vector
  return diffTheta2 + diffPhi2;
}

/// @brief This function finds the top and bottom end of a detector segment in
/// local coordinates
///
/// @param [in] local Local position of the cluster
/// @param [in] segment Segmentation of the detector element
///
/// @return Pair containing the top and bottom end
std::pair<Acts::Vector2D, Acts::Vector2D> findLocalTopAndBottomEnd(
    const Acts::Vector2D& local, const Acts::CartesianSegmentation* segment) {
  auto& binData = segment->binUtility().binningData();
  auto& boundariesX = binData[0].boundaries();
  auto& boundariesY = binData[1].boundaries();

  // Search the x-/y-bin of the cluster
  size_t binX = binData[0].searchLocal(local);
  size_t binY = binData[1].searchLocal(local);

  // Storage of the local top (first) and bottom (second) end
  std::pair<Acts::Vector2D, Acts::Vector2D> topBottomLocal;

  if (boundariesX[binX + 1] - boundariesX[binX] <
      boundariesY[binY + 1] - boundariesY[binY]) {
    // Set the top and bottom end of the strip in local coordinates
    topBottomLocal.first = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                            boundariesY[binY + 1]};
    topBottomLocal.second = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                             boundariesY[binY]};
  } else {
    // Set the top and bottom end of the strip in local coordinates
    topBottomLocal.first = {boundariesX[binX],
                            (boundariesY[binY] + boundariesY[binY + 1]) / 2};
    topBottomLocal.second = {boundariesX[binX + 1],
                             (boundariesY[binY] + boundariesY[binY + 1]) / 2};
  }
  return topBottomLocal;
}

/// @brief Calculates a space point whithout using the vertex
/// @note This is mostly to resolve space points from cosmic data
/// @param a vector to the top end of the first SDE
/// @param c vector to the top end of the second SDE
/// @param q vector from the bottom to the top end of the first SDE
/// @param r vector from the bottom to the top end of the second SDE
/// @return parameter that indicates the location of the space point; returns
/// 1. if it failed
/// @note The meaning of the parameter is explained in more detail in the
/// function body
double calcPerpendicularProjection(const Acts::Vector3D& a,
                                   const Acts::Vector3D& c,
                                   const Acts::Vector3D& q,
                                   const Acts::Vector3D& r) {
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
  Acts::Vector3D ac = c - a;
  double qr = q.dot(r);
  double denom = q.dot(q) - qr * qr;

  // Check for numerical stability
  if (fabs(denom) > 10e-7) {
    // Return lambda0
    return (ac.dot(r) * qr - ac.dot(q) * r.dot(r)) / denom;
  }
  // lambda0 is in the interval [-1,0]. This return serves as error check.
  return 1.;
}

/// @brief This function tests if a space point can be estimated by a more
/// tolerant treatment of construction. In fact, this function indirectly
/// allows shifts of the vertex.
///
/// @param [in] spaPoPa container that stores geometric parameters and rules of
/// the space point formation
/// @param [in] stripLengthGapTolerance Tolerance scaling factor of the gap
/// between strip detector elements
///
/// @return indicator if the test was successful
bool recoverSpacePoint(SpacePointParameters& spaPoPa,
                       double stripLengthGapTolerance) {
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (stripLengthGapTolerance <= 0.) {
    return false;
  }
  spaPoPa.qmag = spaPoPa.q.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spaPoPa.limitExtended =
      spaPoPa.limit + stripLengthGapTolerance / spaPoPa.qmag;
  // Check if m is just slightly outside
  if (fabs(spaPoPa.m) > spaPoPa.limitExtended) {
    return false;
  }
  // Calculate n if not performed previously
  if (spaPoPa.n == 0.) {
    spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs);
  }
  // Check if n is just slightly outside
  if (fabs(spaPoPa.n) > spaPoPa.limitExtended) {
    return false;
  }

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
  double secOnFirstScale =
      spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.qmag * spaPoPa.qmag);
  // Check if both overshoots are in the same direction
  if (spaPoPa.m > 1. && spaPoPa.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spaPoPa.m - 1.;
    double nOvershoot =
        (spaPoPa.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
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
    double nOvershoot =
        -(spaPoPa.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m += biggerOvershoot;
    spaPoPa.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // No solution could be found
  return false;
}

/// @brief This function performs a straight forward calculation of a space
/// point and returns whether it was succesful or not.
///
/// @param [in] stripEnds1 Top and bottom end of the first strip detector
/// element
/// @param [in] stripEnds1 Top and bottom end of the second strip detector
/// element
/// @param [in] posVertex Position of the vertex
/// @param [in, out] spaPoPa Data container of the calculations
/// @param [in] stripLengthTolerance Tolerance scaling factor on the strip
/// detector element length
///
/// @return Boolean statement whether the space point calculation was succesful
bool calculateSpacePoint(
    const std::pair<Acts::Vector3D, Acts::Vector3D>& stripEnds1,
    const std::pair<Acts::Vector3D, Acts::Vector3D>& stripEnds2,
    const Acts::Vector3D& posVertex, SpacePointParameters& spaPoPa,
    const double stripLengthTolerance) {
  /// The following algorithm is meant for finding the position on the first
  /// strip if there is a corresponding cluster on the second strip. The
  /// resulting point is a point x on the first surfaces. This point is
  /// along a line between the points a (top end of the strip)
  /// and b (bottom end of the strip). The location can be parametrized as
  /// 	2 * x = (1 + m) a + (1 - m) b
  /// as function of the scalar m. m is a parameter in the interval
  /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
  /// from the vertex to the cluster on the second strip y is needed to be a
  /// multiple k of the vector from vertex to the hit on the first strip x.
  /// As a consequence of this demand y = k * x needs to be on the
  /// connecting line between the top (c) and bottom (d) end of
  /// the second strip. If both clusters correspond to each other, the
  /// condition
  /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
  /// needs to be fulfilled. Inserting the first equation into this
  /// equation leads to the condition for m as given in the following
  /// algorithm and therefore to the calculation of x.
  /// The same calculation can be repeated for y. Its corresponding
  /// parameter will be named n.

  spaPoPa.s = stripEnds1.first + stripEnds1.second - 2 * posVertex;
  spaPoPa.t = stripEnds2.first + stripEnds2.second - 2 * posVertex;
  spaPoPa.qs = spaPoPa.q.cross(spaPoPa.s);
  spaPoPa.rt = spaPoPa.r.cross(spaPoPa.t);
  spaPoPa.m = -spaPoPa.s.dot(spaPoPa.rt) / spaPoPa.q.dot(spaPoPa.rt);

  // Set the limit for the parameter
  if (spaPoPa.limit == 1. && stripLengthTolerance != 0.) {
    spaPoPa.limit = 1. + stripLengthTolerance;
  }

  // Check if m and n can be resolved in the interval (-1, 1)
  return (fabs(spaPoPa.m) <= spaPoPa.limit &&
          fabs(spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) /
                           spaPoPa.r.dot(spaPoPa.qs)) <= spaPoPa.limit);
}
}  // namespace

///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///
template <typename Cluster>
Acts::SpacePointBuilder<Acts::SpacePoint<Cluster>>::SpacePointBuilder(
    DoubleHitSpacePointConfig cfg)
    : m_cfg(std::move(cfg)) {}

template <typename Cluster>
Acts::Vector2D Acts::SpacePointBuilder<Acts::SpacePoint<Cluster>>::localCoords(
    const Cluster& cluster) const {
  // Local position information
  auto par = cluster.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

template <typename Cluster>
Acts::Vector3D Acts::SpacePointBuilder<Acts::SpacePoint<Cluster>>::globalCoords(
    const GeometryContext& gctx, const Cluster& cluster) const {
  // Receive corresponding surface
  auto& clusterSurface = cluster.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(gctx, localCoords(cluster), mom, pos);

  return pos;
}

template <typename Cluster>
void Acts::SpacePointBuilder<Acts::SpacePoint<Cluster>>::makeClusterPairs(
    const GeometryContext& gctx,
    const std::vector<const Cluster*>& clustersFront,
    const std::vector<const Cluster*>& clustersBack,
    std::vector<std::pair<const Cluster*, const Cluster*>>& clusterPairs)
    const {
  // Return if no clusters are given in a vector
  if (clustersFront.empty() || clustersBack.empty()) {
    return;
  }

  // Declare helper variables
  double currentDiff;
  double diffMin;
  unsigned int clusterMinDist;

  // Walk through all clusters on both surfaces
  for (unsigned int iClustersFront = 0; iClustersFront < clustersFront.size();
       iClustersFront++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of clusters
    clusterMinDist = clustersBack.size();
    for (unsigned int iClustersBack = 0; iClustersBack < clustersBack.size();
         iClustersBack++) {
      // Calculate the distances between the hits
      currentDiff = differenceOfClustersChecked(
          globalCoords(gctx, *(clustersFront[iClustersFront])),
          globalCoords(gctx, *(clustersBack[iClustersBack])), m_cfg.vertex,
          m_cfg.diffDist, m_cfg.diffPhi2, m_cfg.diffTheta2);
      // Store the closest clusters (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        clusterMinDist = iClustersBack;
      }
    }

    // Store the best (=closest) result
    if (clusterMinDist < clustersBack.size()) {
      std::pair<const Cluster*, const Cluster*> clusterPair;
      clusterPair = std::make_pair(clustersFront[iClustersFront],
                                   clustersBack[clusterMinDist]);
      clusterPairs.push_back(clusterPair);
    }
  }
}

template <typename Cluster>
std::pair<Acts::Vector3D, Acts::Vector3D>
Acts::SpacePointBuilder<Acts::SpacePoint<Cluster>>::endsOfStrip(
    const GeometryContext& gctx, const Cluster& cluster) const {
  // Calculate the local coordinates of the cluster
  const Acts::Vector2D local = localCoords(cluster);

  // Receive the binning
  auto segment = dynamic_cast<const Acts::CartesianSegmentation*>(
      &(cluster.digitizationModule()->segmentation()));

  std::pair<Vector2D, Vector2D> topBottomLocal =
      findLocalTopAndBottomEnd(local, segment);

  // Calculate the global coordinates of the top and bottom end of the strip
  Acts::Vector3D topGlobal, bottomGlobal, mom;  // mom is a dummy variable
  const auto* sur = &cluster.referenceSurface();
  sur->localToGlobal(gctx, topBottomLocal.first, mom, topGlobal);
  sur->localToGlobal(gctx, topBottomLocal.second, mom, bottomGlobal);

  // Return the top and bottom end of the strip in global coordinates
  return std::make_pair(topGlobal, bottomGlobal);
}

template <typename Cluster>
void Acts::SpacePointBuilder<Acts::SpacePoint<Cluster>>::calculateSpacePoints(
    const GeometryContext& gctx,
    const std::vector<std::pair<const Cluster*, const Cluster*>>& clusterPairs,
    std::vector<Acts::SpacePoint<Cluster>>& spacePoints) const {
  /// Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  SpacePointParameters spaPoPa;

  // Walk over every found candidate pair
  for (const auto& cp : clusterPairs) {
    // Calculate the ends of the SDEs
    const auto& ends1 = endsOfStrip(gctx, *(cp.first));
    const auto& ends2 = endsOfStrip(gctx, *(cp.second));

    spaPoPa.q = ends1.first - ends1.second;
    spaPoPa.r = ends2.first - ends2.second;

    // Fast skipping if a perpendicular projection should be used
    double resultPerpProj;
    if (m_cfg.usePerpProj) {
      resultPerpProj = calcPerpendicularProjection(ends1.first, ends2.first,
                                                   spaPoPa.q, spaPoPa.r);
      if (resultPerpProj <= 0.) {
        Acts::SpacePoint<Cluster> sp;
        sp.clusterModule.push_back(cp.first);
        sp.clusterModule.push_back(cp.second);
        Vector3D pos = ends1.first + resultPerpProj * spaPoPa.q;
        // TODO: Clusters should deliver timestamp
        sp.spacePoint = {pos.x(), pos.y(), pos.z(), 0.};
        spacePoints.push_back(std::move(sp));
        continue;
      }
    }

    if (calculateSpacePoint(ends1, ends2, m_cfg.vertex, spaPoPa,
                            m_cfg.stripLengthTolerance)) {
      // Store the space point
      Acts::SpacePoint<Cluster> sp;
      sp.clusterModule.push_back(cp.first);
      sp.clusterModule.push_back(cp.second);
      Vector3D pos = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
      // TODO: Clusters should deliver timestamp
      sp.spacePoint = {pos.x(), pos.y(), pos.z(), 0.};
      spacePoints.push_back(std::move(sp));
    } else {
      /// If this point is reached then it was not possible to resolve both
      /// points such that they are on their SDEs
      /// The following code treats a possible recovery of points resolved
      /// slightly outside of the SDE.
      /// @note This procedure is an indirect variation of the vertex
      /// position.
      // Check if a recovery the point(s) and store them if successful
      if (recoverSpacePoint(spaPoPa, m_cfg.stripLengthGapTolerance)) {
        Acts::SpacePoint<Cluster> sp;
        sp.clusterModule.push_back(cp.first);
        sp.clusterModule.push_back(cp.second);
        Vector3D pos =
            0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
        // TODO: Clusters should deliver timestamp
        sp.spacePoint = {pos.x(), pos.y(), pos.z(), 0.};
        spacePoints.push_back(std::move(sp));
      }
    }
  }
}