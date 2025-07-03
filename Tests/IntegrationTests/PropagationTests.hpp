// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <utility>

inline std::shared_ptr<const Acts::TrackingGeometry> createDenseBlock(
    const Acts::GeometryContext& geoCtx) {
  using namespace Acts;
  using namespace UnitLiterals;

  CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0., 0., 0.};
  vConf.length = {4_m, 4_m, 4_m};
  vConf.volumeMaterial =
      std::make_shared<const HomogeneousVolumeMaterial>(Test::makeBeryllium());
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg.push_back(vConf);
  conf.position = {0., 0., 0.};
  conf.length = {4_m, 4_m, 4_m};
  CuboidVolumeBuilder cvb(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });

  return TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);
}

inline std::tuple<std::shared_ptr<const Acts::TrackingGeometry>,
                  std::vector<const Acts::Surface*>>
createDenseTelescope(const Acts::GeometryContext& geoCtx,
                     Acts::Material material, double thickness) {
  using namespace Acts;
  using namespace UnitLiterals;

  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {4_m, 2_m, 2_m};

  {
    CuboidVolumeBuilder::VolumeConfig start;
    start.position = {-1_m, 0, 0};
    start.length = {1_m, 1_m, 1_m};
    start.name = "start";

    conf.volumeCfg.push_back(start);
  }

  if (thickness < 1_m) {
    CuboidVolumeBuilder::VolumeConfig gap;
    gap.position = {-0.25 * (1_m + thickness), 0, 0};
    gap.length = {0.5 * (1_m - thickness), 1_m, 1_m};
    gap.name = "gap1";

    conf.volumeCfg.push_back(gap);
  }

  {
    CuboidVolumeBuilder::VolumeConfig dense;
    dense.position = {0, 0, 0};
    dense.length = {thickness, 1_m, 1_m};
    dense.volumeMaterial =
        std::make_shared<const HomogeneousVolumeMaterial>(material);
    dense.name = "dense";

    conf.volumeCfg.push_back(dense);
  }

  if (thickness < 1_m) {
    CuboidVolumeBuilder::VolumeConfig gap;
    gap.position = {0.25 * (1_m + thickness), 0, 0};
    gap.length = {0.5 * (1_m - thickness), 1_m, 1_m};
    gap.name = "gap2";

    conf.volumeCfg.push_back(gap);
  }

  {
    CuboidVolumeBuilder::SurfaceConfig surface;
    surface.position = {1.5_m, 0, 0};
    surface.rotation =
        Eigen::AngleAxisd(std::numbers::pi / 2, Eigen::Vector3d::UnitY())
            .matrix();
    surface.rBounds = std::make_shared<RectangleBounds>(1_m, 1_m);

    CuboidVolumeBuilder::LayerConfig layer;
    layer.surfaceCfg.push_back(surface);

    CuboidVolumeBuilder::VolumeConfig end;
    end.position = {1_m, 0, 0};
    end.length = {1_m, 1_m, 1_m};
    end.layerCfg.push_back(layer);
    end.name = "end";

    conf.volumeCfg.push_back(end);
  }

  CuboidVolumeBuilder cvb(conf);

  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  auto detector = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);

  std::vector<const Surface*> surfaces;
  detector->visitSurfaces(
      [&](const Surface* surface) { surfaces.push_back(surface); });

  return {std::move(detector), std::move(surfaces)};
}

// parameter construction helpers

/// Construct (initial) curvilinear parameters.
inline Acts::BoundTrackParameters makeParametersCurvilinear(double phi,
                                                            double theta,
                                                            double absMom,
                                                            double charge) {
  using namespace Acts;
  using namespace UnitLiterals;

  Vector4 pos4 = Vector4::Zero();
  auto particleHypothesis = ParticleHypothesis::pionLike(std::abs(charge));
  return BoundTrackParameters::createCurvilinear(
      pos4, phi, theta, particleHypothesis.qOverP(absMom, charge), std::nullopt,
      particleHypothesis);
}

/// Construct (initial) curvilinear parameters with covariance.
inline Acts::BoundTrackParameters makeParametersCurvilinearWithCovariance(
    double phi, double theta, double absMom, double charge) {
  using namespace Acts;
  using namespace UnitLiterals;

  BoundVector stddev = BoundVector::Zero();
  // TODO use momentum-dependent resolutions
  stddev[eBoundLoc0] = 15_um;
  stddev[eBoundLoc1] = 80_um;
  stddev[eBoundTime] = 20_ps;
  stddev[eBoundPhi] = 20_mrad;
  stddev[eBoundTheta] = 30_mrad;
  stddev[eBoundQOverP] = 1_e / 10_GeV;
  BoundSquareMatrix corr = BoundSquareMatrix::Identity();
  corr(eBoundLoc0, eBoundLoc1) = corr(eBoundLoc1, eBoundLoc0) = 0.125;
  corr(eBoundLoc0, eBoundPhi) = corr(eBoundPhi, eBoundLoc0) = 0.25;
  corr(eBoundLoc1, eBoundTheta) = corr(eBoundTheta, eBoundLoc1) = -0.25;
  corr(eBoundTime, eBoundQOverP) = corr(eBoundQOverP, eBoundTime) = 0.125;
  corr(eBoundPhi, eBoundTheta) = corr(eBoundTheta, eBoundPhi) = -0.25;
  corr(eBoundPhi, eBoundQOverP) = corr(eBoundPhi, eBoundQOverP) = -0.125;
  corr(eBoundTheta, eBoundQOverP) = corr(eBoundTheta, eBoundQOverP) = 0.5;
  BoundSquareMatrix cov = stddev.asDiagonal() * corr * stddev.asDiagonal();

  Vector4 pos4 = Vector4::Zero();
  auto particleHypothesis = ParticleHypothesis::pionLike(std::abs(charge));
  return BoundTrackParameters::createCurvilinear(
      pos4, phi, theta, particleHypothesis.qOverP(absMom, charge), cov,
      particleHypothesis);
}

/// Construct (initial) neutral curvilinear parameters.
inline Acts::BoundTrackParameters makeParametersCurvilinearNeutral(
    double phi, double theta, double absMom) {
  using namespace Acts;
  using namespace UnitLiterals;

  Vector4 pos4 = Vector4::Zero();
  return BoundTrackParameters::createCurvilinear(
      pos4, phi, theta, 1 / absMom, std::nullopt, ParticleHypothesis::pion0());
}

// helpers to compare track parameters

/// Check that two parameters object are consistent within the tolerances.
///
/// \warning Does not check that they are defined on the same surface.
inline void checkParametersConsistency(const Acts::BoundTrackParameters& cmp,
                                       const Acts::BoundTrackParameters& ref,
                                       const Acts::GeometryContext& geoCtx,
                                       double epsPos, double epsDir,
                                       double epsMom) {
  using namespace Acts;

  // check stored parameters
  CHECK_CLOSE_ABS(cmp.template get<eBoundLoc0>(),
                  ref.template get<eBoundLoc0>(), epsPos);
  CHECK_CLOSE_ABS(cmp.template get<eBoundLoc1>(),
                  ref.template get<eBoundLoc1>(), epsPos);
  CHECK_CLOSE_ABS(cmp.template get<eBoundTime>(),
                  ref.template get<eBoundTime>(), epsPos);
  // check phi equivalence with circularity
  CHECK_SMALL(detail::radian_sym(cmp.template get<eBoundPhi>() -
                                 ref.template get<eBoundPhi>()),
              epsDir);
  CHECK_CLOSE_ABS(cmp.template get<eBoundTheta>(),
                  ref.template get<eBoundTheta>(), epsDir);
  CHECK_CLOSE_ABS(cmp.template get<eBoundQOverP>(),
                  ref.template get<eBoundQOverP>(), epsMom);
  // check derived parameters
  CHECK_CLOSE_ABS(cmp.position(geoCtx), ref.position(geoCtx), epsPos);
  CHECK_CLOSE_ABS(cmp.time(), ref.time(), epsPos);
  CHECK_CLOSE_ABS(cmp.direction(), ref.direction(), epsDir);
  CHECK_CLOSE_ABS(cmp.absoluteMomentum(), ref.absoluteMomentum(), epsMom);
  // charge should be identical not just similar
  BOOST_CHECK_EQUAL(cmp.charge(), ref.charge());
}

/// Check that two parameters covariances are consistent within the tolerances.
///
/// \warning Does not check that the parameters value itself are consistent.
inline void checkCovarianceConsistency(const Acts::BoundTrackParameters& cmp,
                                       const Acts::BoundTrackParameters& ref,
                                       double relativeTolerance) {
  // either both or none have covariance set
  if (cmp.covariance().has_value()) {
    // comparison parameters have covariance but the reference does not
    BOOST_CHECK(ref.covariance().has_value());
  }
  if (ref.covariance().has_value()) {
    // reference parameters have covariance but the comparison does not
    BOOST_CHECK(cmp.covariance().has_value());
  }
  if (cmp.covariance().has_value() && ref.covariance().has_value()) {
    CHECK_CLOSE_COVARIANCE(cmp.covariance().value(), ref.covariance().value(),
                           relativeTolerance);
  }
}

// helpers to construct target surfaces from track states

/// Construct the transformation from the curvilinear to the global coordinates.
inline Acts::Transform3 createCurvilinearTransform(
    const Acts::BoundTrackParameters& params,
    const Acts::GeometryContext& geoCtx) {
  using namespace Acts;

  Vector3 unitW = params.direction();
  auto [unitU, unitV] = createCurvilinearUnitVectors(unitW);

  RotationMatrix3 rotation = RotationMatrix3::Zero();
  rotation.col(0) = unitU;
  rotation.col(1) = unitV;
  rotation.col(2) = unitW;
  Translation3 offset(params.position(geoCtx));
  Transform3 toGlobal = offset * rotation;

  return toGlobal;
}

/// Construct a z-cylinder centered at zero with the track on its surface.
struct ZCylinderSurfaceBuilder {
  std::shared_ptr<Acts::CylinderSurface> operator()(
      const Acts::BoundTrackParameters& params,
      const Acts::GeometryContext& geoCtx) {
    using namespace Acts;

    auto radius = params.position(geoCtx).template head<2>().norm();
    auto halfz = std::numeric_limits<double>::max();
    return Surface::makeShared<CylinderSurface>(Transform3::Identity(), radius,
                                                halfz);
  }
};

/// Construct a disc at track position with plane normal along track tangent.
struct DiscSurfaceBuilder {
  std::shared_ptr<Acts::DiscSurface> operator()(
      const Acts::BoundTrackParameters& params,
      const Acts::GeometryContext& geoCtx) {
    using namespace Acts;
    using namespace UnitLiterals;

    auto cl = createCurvilinearTransform(params, geoCtx);
    // shift the origin of the plane so the local particle position does not
    // sit directly at the rho=0,phi=undefined singularity
    // TODO this is a hack do avoid issues with the numerical covariance
    //      transport that does not work well at rho=0,
    Vector3 localOffset = Vector3::Zero();
    localOffset[ePos0] = 1_cm;
    localOffset[ePos1] = -1_cm;
    Vector3 globalOriginDelta = cl.linear() * localOffset;
    cl.pretranslate(globalOriginDelta);

    return Surface::makeShared<DiscSurface>(cl);
  }
};

/// Construct a plane at track position with plane normal along track tangent.
struct PlaneSurfaceBuilder {
  std::shared_ptr<Acts::PlaneSurface> operator()(
      const Acts::BoundTrackParameters& params,
      const Acts::GeometryContext& geoCtx) {
    using namespace Acts;

    return Surface::makeShared<PlaneSurface>(
        createCurvilinearTransform(params, geoCtx));
  }
};

/// Construct a z-straw at the track position.
struct ZStrawSurfaceBuilder {
  std::shared_ptr<Acts::StrawSurface> operator()(
      const Acts::BoundTrackParameters& params,
      const Acts::GeometryContext& geoCtx) {
    using namespace Acts;

    return Surface::makeShared<StrawSurface>(
        Transform3(Translation3(params.position(geoCtx))));
  }
};

// helper functions to run the propagation with additional checks

/// Propagate the initial parameters for the given pathlength in space.
///
/// Use a negative path length to indicate backward propagation.
template <typename propagator_t,
          typename options_t = typename propagator_t::template Options<>>
inline std::pair<Acts::BoundTrackParameters, double> transportFreely(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::BoundTrackParameters& initialParams, double pathLength) {
  using namespace Acts;
  using namespace UnitLiterals;

  // setup propagation options
  options_t options(geoCtx, magCtx);
  options.direction = Direction::fromScalar(pathLength);
  options.pathLimit = pathLength;
  options.surfaceTolerance = 1_nm;
  options.stepping.stepTolerance = 1_nm;

  auto result = propagator.propagate(initialParams, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  return {*result.value().endParameters, result.value().pathLength};
}

/// Propagate the initial parameters to the target surface.
template <typename propagator_t,
          typename options_t = typename propagator_t::template Options<>>
inline std::pair<Acts::BoundTrackParameters, double> transportToSurface(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::BoundTrackParameters& initialParams,
    const Acts::Surface& targetSurface, double pathLimit) {
  using namespace Acts;
  using namespace UnitLiterals;

  // setup propagation options
  options_t options(geoCtx, magCtx);
  options.direction = Direction::Forward();
  options.pathLimit = pathLimit;
  options.surfaceTolerance = 1_nm;
  options.stepping.stepTolerance = 1_nm;

  auto result = propagator.propagate(initialParams, targetSurface, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  return {*result.value().endParameters, result.value().pathLength};
}

// self-consistency tests for a single propagator

/// Propagate the initial parameters the given path length along its
/// trajectory and then propagate the final parameters back. Verify that the
/// propagated parameters match the initial ones.
template <typename propagator_t,
          typename options_t = typename propagator_t::template Options<>>
inline void runForwardBackwardTest(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::BoundTrackParameters& initialParams, double pathLength,
    double epsPos, double epsDir, double epsMom) {
  // propagate parameters Acts::Direction::Forward()
  auto [fwdParams, fwdPathLength] = transportFreely<propagator_t, options_t>(
      propagator, geoCtx, magCtx, initialParams, pathLength);
  CHECK_CLOSE_ABS(fwdPathLength, pathLength, epsPos);
  // propagate propagated parameters back again
  auto [bwdParams, bwdPathLength] = transportFreely<propagator_t, options_t>(
      propagator, geoCtx, magCtx, fwdParams, -pathLength);
  CHECK_CLOSE_ABS(bwdPathLength, -pathLength, epsPos);
  // check that initial and back-propagated parameters match
  checkParametersConsistency(initialParams, bwdParams, geoCtx, epsPos, epsDir,
                             epsMom);
}

/// Propagate the initial parameters once for the given path length and
/// use the propagated parameters to define a target surface. Propagate the
/// initial parameters again to the target surface. Verify that the surface has
/// been found and the parameters are consistent.
template <typename propagator_t, typename surface_builder_t,
          typename options_t = typename propagator_t::template Options<>>
inline void runToSurfaceTest(const propagator_t& propagator,
                             const Acts::GeometryContext& geoCtx,
                             const Acts::MagneticFieldContext& magCtx,
                             const Acts::BoundTrackParameters& initialParams,
                             double pathLength,
                             surface_builder_t&& buildTargetSurface,
                             double epsPos, double epsDir, double epsMom) {
  // free propagation for the given path length
  auto [freeParams, freePathLength] = transportFreely<propagator_t, options_t>(
      propagator, geoCtx, magCtx, initialParams, pathLength);
  CHECK_CLOSE_ABS(freePathLength, pathLength, epsPos);
  // build a target surface at the propagated position
  auto surface = buildTargetSurface(freeParams, geoCtx);
  BOOST_CHECK(surface);

  // bound propagation onto the target surface
  // increase path length limit to ensure the surface can be reached
  auto [surfParams, surfPathLength] =
      transportToSurface<propagator_t, options_t>(propagator, geoCtx, magCtx,
                                                  initialParams, *surface,
                                                  1.5 * pathLength);
  CHECK_CLOSE_ABS(surfPathLength, pathLength, epsPos);

  // check that the to-surface propagation matches the defining free parameters
  CHECK_CLOSE_ABS(surfParams.position(geoCtx), freeParams.position(geoCtx),
                  epsPos);
  CHECK_CLOSE_ABS(surfParams.time(), freeParams.time(), epsPos);
  CHECK_CLOSE_ABS(surfParams.direction(), freeParams.direction(), epsDir);
  CHECK_CLOSE_ABS(surfParams.absoluteMomentum(), freeParams.absoluteMomentum(),
                  epsMom);
  CHECK_CLOSE_ABS(surfPathLength, freePathLength, epsPos);
}

// consistency tests between two propagators

/// Propagate the initial parameters along their trajectory for the given path
/// length using two different propagators and verify consistent output.
template <typename cmp_propagator_t, typename ref_propagator_t>
inline void runForwardComparisonTest(
    const cmp_propagator_t& cmpPropagator,
    const ref_propagator_t& refPropagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::BoundTrackParameters& initialParams, double pathLength,
    double epsPos, double epsDir, double epsMom, double tolCov) {
  // propagate twice using the two different propagators
  auto [cmpParams, cmpPath] = transportFreely<cmp_propagator_t>(
      cmpPropagator, geoCtx, magCtx, initialParams, pathLength);
  auto [refParams, refPath] = transportFreely<ref_propagator_t>(
      refPropagator, geoCtx, magCtx, initialParams, pathLength);
  // check parameter comparison
  checkParametersConsistency(cmpParams, refParams, geoCtx, epsPos, epsDir,
                             epsMom);
  checkCovarianceConsistency(cmpParams, refParams, tolCov);
  CHECK_CLOSE_ABS(cmpPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(refPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(cmpPath, refPath, epsPos);
}

/// Propagate the initial parameters along their trajectory for the given path
/// length using the reference propagator. Use the propagated track parameters
/// to define a target plane. Propagate the initial parameters using two
/// different propagators and verify consistent output.
template <typename cmp_propagator_t, typename ref_propagator_t,
          typename surface_builder_t>
inline void runToSurfaceComparisonTest(
    const cmp_propagator_t& cmpPropagator,
    const ref_propagator_t& refPropagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::BoundTrackParameters& initialParams, double pathLength,
    surface_builder_t&& buildTargetSurface, double epsPos, double epsDir,
    double epsMom, double tolCov) {
  // free propagation with the reference propagator for the given path length
  auto [freeParams, freePathLength] = transportFreely<ref_propagator_t>(
      refPropagator, geoCtx, magCtx, initialParams, pathLength);
  CHECK_CLOSE_ABS(freePathLength, pathLength, epsPos);

  // build a target surface at the propagated position
  auto surface = buildTargetSurface(freeParams, geoCtx);
  BOOST_CHECK(surface);

  // propagate twice to the surface using the two different propagators
  // increase path length limit to ensure the surface can be reached
  auto [cmpParams, cmpPath] = transportToSurface<cmp_propagator_t>(
      cmpPropagator, geoCtx, magCtx, initialParams, *surface, 1.5 * pathLength);
  auto [refParams, refPath] = transportToSurface<ref_propagator_t>(
      refPropagator, geoCtx, magCtx, initialParams, *surface, 1.5 * pathLength);
  // check parameter comparison
  checkParametersConsistency(cmpParams, refParams, geoCtx, epsPos, epsDir,
                             epsMom);
  checkCovarianceConsistency(cmpParams, refParams, tolCov);
  CHECK_CLOSE_ABS(cmpPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(refPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(cmpPath, refPath, epsPos);
}

template <typename propagator_t>
inline void runDenseForwardTest(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx, double p, double q,
    const Acts::Surface& targetSurface, const Acts::Material& material,
    double thickness, const Acts::Logger& logger = Acts::getDummyLogger()) {
  using namespace Acts;
  using namespace UnitLiterals;

  const auto particleHypothesis = ParticleHypothesis::muon();
  const float mass = particleHypothesis.mass();
  const float absQ = std::abs(q);
  const double initialQOverP = particleHypothesis.qOverP(p, q);
  const auto initialParams = BoundTrackParameters::createCurvilinear(
      Vector4(-1.5_m, 0, 0, 0), Vector3::UnitX(), initialQOverP,
      BoundVector::Constant(1e-16).asDiagonal(), particleHypothesis);

  typename propagator_t::template Options<> options(geoCtx, magCtx);
  options.maxSteps = 10000;
  options.stepping.maxStepSize = 100_mm;
  options.stepping.dense.meanEnergyLoss = true;
  options.stepping.dense.momentumCutOff = 1_MeV;

  const auto result =
      propagator.propagate(initialParams, targetSurface, options);

  BOOST_CHECK(result.ok());
  CHECK_CLOSE_REL(3_m, result->pathLength, 1e-6);
  BOOST_CHECK(result->endParameters);

  BoundTrackParameters endParams = result->endParameters.value();

  BOOST_CHECK(endParams.covariance());
  CHECK_CLOSE_ABS(initialParams.position(geoCtx) + Acts::Vector3(3_m, 0, 0),
                  endParams.position(geoCtx), 1e-6);
  CHECK_CLOSE_ABS(initialParams.direction(), endParams.direction(), 1e-6);

  const auto& cov = endParams.covariance().value();

  const double endQOverP = endParams.qOverP();
  const double endP = endParams.absoluteMomentum();
  const double endEloss =
      Acts::fastHypot(p, mass) - Acts::fastHypot(endP, mass);
  const double endErrX = std::sqrt(cov(eBoundLoc0, eBoundLoc0));
  const double endErrY = std::sqrt(cov(eBoundLoc1, eBoundLoc1));
  const double endErrQOverP = std::sqrt(cov(eBoundQOverP, eBoundQOverP));
  const double endErrP =
      (absQ / std::pow(endParams.qOverP(), 2)) * endErrQOverP;
  const double endErrE = (p / Acts::fastHypot(p, mass)) * endErrP;
  const double endErrTheta = std::sqrt(cov(eBoundTheta, eBoundTheta));
  const double endErrPhi = std::sqrt(cov(eBoundPhi, eBoundPhi));

  ACTS_VERBOSE("input q/p = " << initialQOverP);
  ACTS_VERBOSE("input p = " << p);
  ACTS_VERBOSE("output q/p = " << endQOverP);
  ACTS_VERBOSE("output p = " << endP);
  ACTS_VERBOSE("output t = " << endParams.time());
  ACTS_VERBOSE("output eloss = " << endEloss);
  ACTS_VERBOSE("output err x = " << endErrX);
  ACTS_VERBOSE("output err y = " << endErrY);
  ACTS_VERBOSE("output err q/p = " << endErrQOverP);
  ACTS_VERBOSE("output err p = " << endErrP);
  ACTS_VERBOSE("output err E = " << endErrE);
  ACTS_VERBOSE("output err theta = " << endErrTheta);
  ACTS_VERBOSE("output err phi = " << endErrPhi);

  BOOST_CHECK_CLOSE(endErrX, endErrY, 1e-6);
  BOOST_CHECK_CLOSE(endErrTheta, endErrPhi, 1e-6);

  const float elossInitial = computeEnergyLossMean(
      MaterialSlab(material, thickness), particleHypothesis.absolutePdg(), mass,
      initialQOverP, absQ);
  const float elossFinal = computeEnergyLossMean(
      MaterialSlab(material, thickness), particleHypothesis.absolutePdg(), mass,
      endQOverP, absQ);

  ACTS_VERBOSE("ref eloss initial = " << elossInitial);
  ACTS_VERBOSE("ref eloss final = " << elossFinal);

  // energy loss will be smaller than elossInitial because of the increasing
  // radiation losses
  BOOST_CHECK_LT(endEloss, elossInitial);
  BOOST_CHECK_GT(endEloss, elossFinal);

  const float theta0Initial = computeMultipleScatteringTheta0(
      MaterialSlab(material, thickness), particleHypothesis.absolutePdg(), mass,
      initialQOverP, absQ);
  const float theta0Final = computeMultipleScatteringTheta0(
      MaterialSlab(material, thickness), particleHypothesis.absolutePdg(), mass,
      endQOverP, absQ);

  ACTS_VERBOSE("ref theta0 initial = " << theta0Initial);
  ACTS_VERBOSE("ref theta0 final = " << theta0Final);

  // angle errors will be larger than theta0Initial because of the energy loss
  BOOST_CHECK_GE(endErrTheta, theta0Initial);
  BOOST_CHECK_LT(endErrTheta, theta0Final);
  CHECK_CLOSE_REL(endErrTheta, 0.5 * (theta0Initial + theta0Final), 0.1);

  // estimates the positional uncertainty after crossing the material
  // block by integrating the angle variance over substeps
  const auto estimateVarX = [&](double qOverP, std::size_t substeps) {
    double previousTheta0 = 0;
    double varAngle = 0;
    double varPosition = 0;
    double covAnglePosition = 0;

    const auto step = [&](double substep, double theta0) {
      double deltaVarTheta = square(theta0) - square(previousTheta0);
      double deltaVarPos = varAngle * square(substep) +
                           2 * covAnglePosition * substep +
                           deltaVarTheta * (square(substep) / 3);
      double deltaCovAnglePosition =
          varAngle * substep + deltaVarTheta * substep / 2;
      previousTheta0 = theta0;
      varAngle += deltaVarTheta;
      varPosition += deltaVarPos;
      covAnglePosition += deltaCovAnglePosition;
    };

    // step through the material block
    for (std::size_t i = 0; i < substeps; ++i) {
      const float theta0 = computeMultipleScatteringTheta0(
          MaterialSlab(material, (i + 1) * thickness / substeps),
          particleHypothesis.absolutePdg(), mass, qOverP, absQ);
      step(thickness / substeps, theta0);
    }

    // step to the target surface
    step(1_m, previousTheta0);

    return varPosition;
  };

  // as a lower bound for sigma_x we use the initial theta0
  const double lowerBoundVarX = estimateVarX(initialQOverP, 10);
  // as an upper bound for sigma_x we use the final theta0
  const double upperBoundVarX = estimateVarX(endQOverP, 10);

  const double lowerBoundSigmaX = std::sqrt(lowerBoundVarX);
  const double upperBoundSigmaX = std::sqrt(upperBoundVarX);

  ACTS_VERBOSE("ref lower bound sigma x = " << lowerBoundSigmaX);
  ACTS_VERBOSE("ref upper bound sigma x = " << upperBoundSigmaX);

  BOOST_CHECK_GT(endErrX, lowerBoundSigmaX);
  BOOST_CHECK_LT(endErrX, upperBoundSigmaX);
  CHECK_CLOSE_REL(endErrX, 0.5 * (lowerBoundSigmaX + upperBoundSigmaX), 0.2);

  const float qOverPSigma = computeEnergyLossLandauSigmaQOverP(
      MaterialSlab(material, thickness), mass, initialQOverP, absQ);

  ACTS_VERBOSE("ref sigma q/p = " << qOverPSigma);

  CHECK_CLOSE_REL(endErrQOverP, qOverPSigma, 1e-4);
}
