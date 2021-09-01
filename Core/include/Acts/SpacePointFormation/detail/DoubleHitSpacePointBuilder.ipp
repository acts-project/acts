// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include <iostream> // just for tests
#include <cmath>
#include <limits>
#include <variant>

namespace Acts {
namespace detail {
/// @brief Storage container for variables related to the calculation of space
/// points
struct SpacePointParameters {
  /// Vector pointing from bottom to top end of first SDE
  Vector3 q;
  /// Vector pointing from bottom to top end of second SDE
  Vector3 r;
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3 s;
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3 t;
  /// Cross product between SpacePointParameters::q and
  /// SpacePointParameters::s
  Vector3 qs;
  /// Cross product between SpacePointParameters::r and
  /// SpacePointParameters::t
  Vector3 rt;
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

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two measurements
///
/// @param [in] pos1 position of the first measurement
/// @param [in] pos2 position the second measurement
/// @param [in] maxDistance Maximum distance between two measurements
/// @param [in] maxAngleTheta2 Maximum squared theta angle between two
/// measurements
/// @param [in] maxAnglePhi2 Maximum squared phi angle between two measurements
///
/// @return The squared sum within configuration parameters, otherwise -1
inline double differenceOfMeasurementsChecked(const Vector3& pos1,
                                              const Vector3& pos2,
                                              const Vector3& posVertex,
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
/// @param [in] local Local position of the Measurement
/// @param [in] segment Segmentation of the detector element
///
/// @return Pair containing the top and bottom end
inline std::pair<Vector2, Vector2> findLocalTopAndBottomEnd(
    const Vector2& local, const CartesianSegmentation* segment) {
      std::cout << "findLocalTopAndBottomEnd" << std::endl;
  auto& binData = segment->binUtility().binningData();
  auto& boundariesX = binData[0].boundaries();
  auto& boundariesY = binData[1].boundaries();

  // Search the x-/y-bin of the Measurement
  size_t binX = binData[0].searchLocal(local);
  size_t binY = binData[1].searchLocal(local);

  // Storage of the local top (first) and bottom (second) end
  std::pair<Vector2, Vector2> topBottomLocal;
  // topBottomLocal.first;
  // topBottomLocal.second;
  Vector2 topLocal = Vector2(0,0);
  Vector2 bottomLocal = Vector2(0,0);
  //std::cout << "boundary size p" << boundariesX.size() << " "
            //<< boundariesY.size() << std::endl;

  //std::cout << boundariesX.at(0) << std::endl;
  //std::cout << " " << boundariesX[1] << std::endl;
  //std::cout << boundariesY[0] << " " << boundariesY[1] << std::endl;
  ///std::cout << binX << binY << std::endl;
  if (boundariesX[binX + 1] - boundariesX[binX] < boundariesY[binY + 1] - boundariesY[binY]) {
    // Set the top and bottom end of the strip in local coordinates
  //  std::cout << "check01" << std::endl;
  //  std::cout << boundariesX[0] << std::endl;
  //  std::cout << boundariesX[1] << std::endl;
  //  std::cout << boundariesY[0] << std::endl;
  //  std::cout << boundariesY[1] << std::endl;
    topLocal = Vector2((boundariesX[binX] + boundariesX[binX + 1]) / 2, boundariesY[binY]);
    bottomLocal = Vector2((boundariesX[binX] + boundariesX[binX + 1]) / 2,
                          boundariesY[binY + 1]);
//    std::cout << "check1" << std::endl;
    // topBottomLocal = std::make_pair(topLocal,bottomLocal);
//    std::cout << "check1.1" << std::endl;
    // topBottomLocal.first = Vector2((boundariesX[binX] + boundariesX[binX +
    // 1]) / 2,
    //                        boundariesY[binY + 1]);
    //                        std::cout << "check0" << std::endl;
    // topBottomLocal.second = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
    //                       boundariesY[binY]};
    //std::cout << "check1" << std::endl;
  } else {
    // Set the top and bottom end of the strip in local coordinates
    //std::cout << "check2" << std::endl;
    /*     topBottomLocal.first = {boundariesX[binX],
                                (boundariesY[binY] + boundariesY[binY + 1]) /
       2}; std::cout << "check3" << std::endl; topBottomLocal.second =
       {boundariesX[binX + 1],
                                 (boundariesY[binY] + boundariesY[binY + 1]) /
       2}; */
  }

  // std::pair<Vector2, Vector2> topBottomLocal =
  // std::make_pair(topLocal,bottomLocal); std::cout  << "topbottom local" <<
  // topBottomLocal.first << " " << topBottomLocal.second << std::endl;
  std::cout << topLocal[0] << " " << topLocal[1] << std::endl;
  //std::cout << bottomLocal[0] << " " << bottomLocal[1] << std::endl;
  //return {topLocal, bottomLocal};
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
inline double calcPerpendicularProjection(const Vector3& a, const Vector3& c,
                                          const Vector3& q, const Vector3& r) {
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
  Vector3 ac = c - a;
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
inline bool recoverSpacePoint(SpacePointParameters& spaPoPa,
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
inline bool calculateSpacePoint(const std::pair<Vector3, Vector3>& stripEnds1,
                                const std::pair<Vector3, Vector3>& stripEnds2,
                                const Vector3& posVertex,
                                SpacePointParameters& spaPoPa,
                                const double stripLengthTolerance) {
  /// The following algorithm is meant for finding the position on the first
  /// strip if there is a corresponding Measurement on the second strip. The
  /// resulting point is a point x on the first surfaces. This point is
  /// along a line between the points a (top end of the strip)
  /// and b (bottom end of the strip). The location can be parametrized as
  /// 	2 * x = (1 + m) a + (1 - m) b
  /// as function of the scalar m. m is a parameter in the interval
  /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
  /// from the vertex to the Measurement on the second strip y is needed to be a
  /// multiple k of the vector from vertex to the hit on the first strip x.
  /// As a consequence of this demand y = k * x needs to be on the
  /// connecting line between the top (c) and bottom (d) end of
  /// the second strip. If both measurements correspond to each other, the
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

  return true;
}
}  // namespace detail
}  // namespace Acts

///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///
template <typename spacepoint_t, typename cluster_t>
Acts::DoubleHitSpacePointBuilder<spacepoint_t, cluster_t>::
    DoubleHitSpacePointBuilder(DoubleHitSpacePointBuilderConfig cfg)
    : m_cfg(std::move(cfg)) {}

template <typename spacepoint_t, typename cluster_t>
Acts::Vector2
Acts::DoubleHitSpacePointBuilder<spacepoint_t, cluster_t>::localCoords(
    const cluster_t& cluster) const {
  // Local position information
  const auto meas = cluster.measurement();
  // auto par = meas.parameters();

  Acts::Vector2 localPos = std::visit(
      [](const auto& x) {
        // double local = std::visit([](const auto &x){
        Acts::BoundVector par = x.expander() * x.parameters();
        // Acts::Vector2 pos(par[0],par[1]);
        return Acts::Vector2(par[0], par[1]);
      },
      meas);
  // return par[1];},meas);

  // Acts::BoundVector par = std::visit([](const auto& x) { return
  // x.parameters(); }, meas); Acts::Vector2
  // local(par[Acts::BoundIndices::eBoundLoc0],
  // par[Acts::BoundIndices::eBoundLoc1]);
  // Acts::Vector2 local(par[0],par[1]);

  // return local;
  // Acts::Vector2 localPos(0,0);
  return localPos;
}

template <typename spacepoint_t, typename cluster_t>
// Acts::Vector3
std::pair<Acts::Vector3, Acts::Vector2>
Acts::DoubleHitSpacePointBuilder<spacepoint_t, cluster_t>::globalCoords(
    const Acts::GeometryContext& gctx, const cluster_t& clus) const {
  // // Receive corresponding surface
  // auto& measurementSurface = measurement;
  auto meas = clus.measurement();
  auto slink = std::visit([](const auto& x) { return x.sourceLink(); }, meas);
  // auto slink = meas.sourceLink();
  // auto slink  = meas.measurement().sourceLink();
  const auto geoId = slink.geometryId();


  const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
if (surface == nullptr) std::cout << "surface is null" << std::endl;
  
  auto [localPos, localCov] = std::visit(
      [](const auto& measurement) {
        auto expander = measurement.expander();
        // auto indices = measurement.indices();
        Acts::BoundVector par = expander * measurement.parameters();
        //std::cout << "measurement parameters" << std::endl << par << std::endl;
        Acts::BoundSymMatrix cov =
            expander * measurement.covariance() * expander.transpose();
        // extract local position
        Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
        // extract local position covariance.
        Acts::SymMatrix2 lcov =
            cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);
  // std::cout << "local pos:" << std::endl << localPos << std::endl;
  // transform local position to global coordinates
  Acts::Vector3 globalFakeMom(1, 1, 1);
  //std::cout << "tmp" << std::endl;
  //std::cout << geoId << std::endl;
  auto stype = surface->type();
  std::cout << "surface type " << stype << std::endl;

  Acts::Vector3 globalPos = surface->localToGlobal(gctx, localPos, globalFakeMom);
  Acts::RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, globalFakeMom);

  auto x = globalPos[Acts::ePos0];
  auto y = globalPos[Acts::ePos1];
  auto scale = 2 / std::hypot(x, y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Acts::Vector2(var[0], var[1]);
  return std::make_pair(globalPos, gcov);
}

template <typename spacepoint_t, typename cluster_t>
void Acts::DoubleHitSpacePointBuilder<spacepoint_t, cluster_t>::
    makeMeasurementPairs(
        const Acts::GeometryContext& gctx,
        const std::vector<const cluster_t*>& clustersFront,
        const std::vector<const cluster_t*>& clustersBack,
        std::vector<std::pair<const cluster_t*, const cluster_t*>>&
            clusterPairs) const {
  // Return if no Measurements are given in a vector
  if (clustersFront.empty() || clustersBack.empty()) {
    return;
  }
  // std::cout << "clusters front in DHPB" << std::endl;
  // std::cout << "size " << clustersFront.size() << std::endl;
  // auto clus = clustersFront[0];
  // auto meas = clus.measurement();

  // auto slink = std::visit([](const auto& x) { return x.sourceLink(); },
  // meas); std::cout << slink.geometryId() << std::endl; auto [gPos,gCov] =
  // globalCoords(gctx,meas); std::cout << gPos << std::endl;

  // Declare helper variables
  double currentDiff;
  double diffMin;
  unsigned int clusterMinDist;

  // Walk through all Measurements on both surfaces
  for (unsigned int iClustersFront = 0; iClustersFront < clustersFront.size();
       iClustersFront++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of Measurements
    clusterMinDist = clustersBack.size();
    for (unsigned int iClustersBack = 0; iClustersBack < clustersBack.size();
         iClustersBack++) {
      // Calculate the distances between the hits


      // 
      auto clus = *clustersFront[iClustersFront];
      if (clustersFront[iClustersFront] == nullptr) std::cout << "cluster is null" << std::endl;
      else std::cout << "cluster is not null" << std::endl;
auto meas = clus.measurement();
int cidx = clus.index();
std::cout << "cluster index " << cidx << std::endl;

  auto slink = std::visit([](const auto& x) { 
    auto sl = x.sourceLink();
    auto slid = sl.index();
    std::cout << "slink index " << slid << std::endl;
    auto ggg = sl.geometryId();
    std::cout << "ggg " << ggg << std::endl;
    return x.sourceLink(); 
    }, meas);
  // auto slink = meas.sourceLink();
  // auto slink  = meas.measurement().sourceLink();
  const auto geoId = slink.geometryId();
  std::cout << iClustersFront << " " << geoId << std::endl;
      ///
      auto gpos_front = globalCoords(gctx, *clustersFront[iClustersFront]);
      auto gpos_back = globalCoords(gctx, *clustersFront[iClustersFront]);
      currentDiff = detail::differenceOfMeasurementsChecked(
          gpos_front.first, gpos_back.first, m_cfg.vertex, m_cfg.diffDist,
          m_cfg.diffPhi2, m_cfg.diffTheta2);
      // Store the closest Measurements (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        clusterMinDist = iClustersBack;
      }
    }

    // Store the best (=closest) result
    if (clusterMinDist < clustersBack.size()) {
      std::pair<const cluster_t*, const cluster_t*> clusterPair;
      clusterPair = std::make_pair(clustersFront[iClustersFront],
                                   clustersBack[clusterMinDist]);
      clusterPairs.push_back(clusterPair);
    }
  }
}

template <typename spacepoint_t, typename cluster_t>
std::pair<Acts::Vector3, Acts::Vector3>
Acts::DoubleHitSpacePointBuilder<spacepoint_t, cluster_t>::endsOfStrip(
    const Acts::GeometryContext& gctx, const cluster_t& cluster) const {
  // Calculate the local coordinates of the Measurement
  const Acts::Vector2 local = localCoords(cluster);
  //std::cout << "local coords :" << local[0] << " " << local[1] << std::endl;

  // Receive the binning
  const auto segment = dynamic_cast<const Acts::CartesianSegmentation*>(
      &(cluster.segmentation()));
  //&(measurement.digitizationModule()->segmentation()));
  
  //std::cout << std::endl;
//detail::findLocalTopAndBottomEnd(local, segment);
///std::cout << "check2" << std::endl;
  std::pair<Vector2, Vector2> topBottomLocal =
      detail::findLocalTopAndBottomEnd(local, segment);
  std::cout << "topbottom local calculated " ;
  // Calculate the global coordinates of the top and bottom end of the strip
  
  
  Acts::Vector3 fakeMom(1,1,1);
  fakeMom << 3., 2., 1.;  // mom is a dummy variable
  // const auto* sur = &measurement.referenceObject();
  auto meas = cluster.measurement();
  auto slink = std::visit([](const auto& x) { return x.sourceLink(); }, meas);
  const auto geoId = slink.geometryId();
  // const auto* sur = &meas.referenceObject();
  const Acts::Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
  Acts::Vector3 topGlobal =
      surface->localToGlobal(gctx, topBottomLocal.first, fakeMom);
  Acts::Vector3 bottomGlobal =
      surface->localToGlobal(gctx, topBottomLocal.second, fakeMom);

  // Return the top and bottom end of the strip in global coordinates
  return std::make_pair(topGlobal, bottomGlobal);
  //Acts::Vector3 tmp(0,0,0);
  //return std::make_pair(tmp,tmp);
}

template <typename spacepoint_t, typename cluster_t>
void Acts::DoubleHitSpacePointBuilder<spacepoint_t, cluster_t>::
    calculateSpacePoints(
        const Acts::GeometryContext& gctx,
        const std::vector<std::pair<const cluster_t*, const cluster_t*>>&
            measurementPairs,
        std::vector<spacepoint_t>& spacePoints) const {
  /// Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  detail::SpacePointParameters spaPoPa;

  // Walk over every found candidate pair
  for (const auto& cp : measurementPairs) {
    // Calculate the ends of the SDEs
    const auto& ends1 = endsOfStrip(gctx, *(cp.first));
    const auto& ends2 = endsOfStrip(gctx, *(cp.second));

    spaPoPa.q = ends1.first - ends1.second;
    spaPoPa.r = ends2.first - ends2.second;

    // Fast skipping if a perpendicular projection should be used
    double resultPerpProj;
    if (m_cfg.usePerpProj) {
      resultPerpProj = detail::calcPerpendicularProjection(
          ends1.first, ends2.first, spaPoPa.q, spaPoPa.r);
      if (resultPerpProj <= 0.) {
        // spacepoint_t sp;
        // sp.MeasurementModule.push_back(cp.first);
        // sp.MeasurementModule.push_back(cp.second);
        Vector3 pos = ends1.first + resultPerpProj * spaPoPa.q;
        double varRho = -1.;  // TODO impriment variance rho and z
        double varZ = -1;
        size_t measurementIndex = 0;  // should be indices
        auto sp = spacepoint_t(pos, varRho, varZ, measurementIndex);
        spacePoints.push_back(std::move(sp));
        continue;
      }
    }

    if (calculateSpacePoint(ends1, ends2, m_cfg.vertex, spaPoPa,
                            m_cfg.stripLengthTolerance)) {
      // Store the space point
      // spacepoint_t sp;
      // sp.MeasurementModule.push_back(cp.first);
      // sp.MeasurementModule.push_back(cp.second);
      Vector3 pos = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
      double varRho = -1.;  // TODO impriment variance rho and z
      double varZ = -1;
      size_t measurementIndex = 0;  // should be indices
      auto sp = spacepoint_t(pos, varRho, varZ, measurementIndex);
      spacePoints.push_back(std::move(sp));
      // TODO: Measurements should deliver timestamp
      /// sp.vector = pos;
      // spacePoints.push_back(std::move(sp));
    } else {
      /// If this point is reached then it was not possible to resolve both
      /// points such that they are on their SDEs
      /// The following code treats a possible recovery of points resolved
      /// slightly outside of the SDE.
      /// @note This procedure is an indirect variation of the vertex
      /// position.
      // Check if a recovery the point(s) and store them if successful
      if (detail::recoverSpacePoint(spaPoPa, m_cfg.stripLengthGapTolerance)) {
        // spacepoint_t sp;
        // sp.MeasurementModule.push_back(cp.first);
        // sp.MeasurementModule.push_back(cp.second);
        Vector3 pos =
            0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
        double varRho = -1.;  // TODO impriment variance rho and z
        double varZ = -1;
        size_t measurementIndex = 0;  // should be indices
        auto sp = spacepoint_t(pos, varRho, varZ, measurementIndex);
        spacePoints.push_back(std::move(sp));
        // sp.vector = pos;
        // spacePoints.push_back(std::move(sp));
      }
    }
  }
}
