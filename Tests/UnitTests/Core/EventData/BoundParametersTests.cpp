// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ParametersTestHelper.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

/// @brief Unit test for parameters at a plane
///
BOOST_DATA_TEST_CASE(
    bound_to_plane_test,
    bdata::random((bdata::seed = 1240,
                   bdata::distribution =
                       std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 2351,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 3412,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((
            bdata::seed = 5732,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 8941,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 1295,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::xrange(100),
    x, y, z, a, b, c, index) {
  (void)index;
  Vector3D center{x, y, z};
  auto transform = std::make_shared<Transform3D>();
  transform->setIdentity();
  RotationMatrix3D rot;
  rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY()) *
        AngleAxis3D(c, Vector3D::UnitZ());
  transform->prerotate(rot);
  transform->pretranslate(center);
  // create the surfacex
  auto bounds = std::make_shared<RectangleBounds>(100., 100.);
  auto pSurface = Surface::makeShared<PlaneSurface>(transform, bounds);

  // now create parameters on this surface
  // l_x, l_y, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {
      {-0.1234, 9.8765, 0.45, 0.888, 0.001, 21.}};
  SingleBoundTrackParameters<ChargedPolicy>::ParVector_t pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  const double phi = pars_array[2];
  const double theta = pars_array[3];
  double p = fabs(1. / pars_array[4]);
  Vector3D direction(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  Vector3D mom = p * direction;
  // the global position
  Vector3D pos =
      center + pars_array[0] * rot.col(0) + pars_array[1] * rot.col(1);
  // constructor from parameter vector
  BoundParameters ataPlane_from_pars(tgContext, std::nullopt, pars, pSurface);
  consistencyCheck(ataPlane_from_pars, pos, mom, 1., 21., pars_array);
  // constructor from global parameters
  BoundParameters ataPlane_from_global(tgContext, std::nullopt, pos, mom, 1.,
                                       21., pSurface);
  consistencyCheck(ataPlane_from_global, pos, mom, 1., 21., pars_array);
  // constructor for neutral parameters
  NeutralBoundTrackParameters n_ataPlane_from_pars(tgContext, std::nullopt,
                                                   pars, pSurface);
  consistencyCheck(n_ataPlane_from_pars, pos, mom, 0., 21., pars_array);
  // constructor for neutral global parameters
  NeutralBoundTrackParameters n_ataPlane_from_global(tgContext, std::nullopt,
                                                     pos, mom, 21., pSurface);
  consistencyCheck(n_ataPlane_from_global, pos, mom, 0., 21., pars_array);

  // check shared ownership of same surface
  BOOST_CHECK_EQUAL(&ataPlane_from_pars.referenceSurface(), pSurface.get());
  BOOST_CHECK_EQUAL(&ataPlane_from_pars.referenceSurface(),
                    &ataPlane_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(&n_ataPlane_from_pars.referenceSurface(), pSurface.get());
  BOOST_CHECK_EQUAL(&n_ataPlane_from_pars.referenceSurface(),
                    &n_ataPlane_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(pSurface.use_count(), 5u);

  // check that the reference frame is the rotation matrix
  CHECK_CLOSE_REL(ataPlane_from_pars.referenceFrame(tgContext), rot, 1e-6);

  /// modification test via setter functions
  double ux = 0.3;
  double uy = 0.4;

  ataPlane_from_pars.set<Acts::eLOC_X>(tgContext, ux);
  ataPlane_from_pars.set<Acts::eLOC_Y>(tgContext, uy);
  // we should have a new updated position
  Vector3D lPosition3D(ux, uy, 0.);
  Vector3D uposition = rot * lPosition3D + center;
  CHECK_CLOSE_REL(uposition, ataPlane_from_pars.position(), 1e-6);

  double uphi = 1.2;
  double utheta = 0.2;
  double uqop = 0.025;
  double ut = 1337.;

  ataPlane_from_pars.set<Acts::ePHI>(tgContext, uphi);
  ataPlane_from_pars.set<Acts::eTHETA>(tgContext, utheta);
  ataPlane_from_pars.set<Acts::eQOP>(tgContext, uqop);
  ataPlane_from_pars.set<Acts::eT>(tgContext, ut);
  // we should have a new updated momentum
  Vector3D umomentum = 40. * Vector3D(cos(uphi) * sin(utheta),
                                      sin(uphi) * sin(utheta), cos(utheta));

  CHECK_CLOSE_REL(umomentum, ataPlane_from_pars.momentum(), 1e-6);
}

/// @brief Unit test for parameters at a disc
///
BOOST_DATA_TEST_CASE(
    bound_to_disc_test,
    bdata::random((bdata::seed = 9810,
                   bdata::distribution =
                       std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 1221,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 12132,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((
            bdata::seed = 16783,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 13984,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 77615,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::xrange(100),
    x, y, z, a, b, c, index) {
  (void)index;
  Vector3D center{x, y, z};
  auto transform = std::make_shared<Transform3D>();
  transform->setIdentity();
  RotationMatrix3D rot;
  rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY()) *
        AngleAxis3D(c, Vector3D::UnitZ());
  transform->prerotate(rot);
  transform->pretranslate(center);

  auto bounds = std::make_shared<RadialBounds>(100., 1200.);
  auto dSurface = Surface::makeShared<DiscSurface>(transform, bounds);

  // now create parameters on this surface
  // r, phi, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {{125., 0.345, 0.45, 0.888, 0.001, 21.}};
  SingleBoundTrackParameters<ChargedPolicy>::ParVector_t pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  const double phi = pars_array[2];
  const double theta = pars_array[3];
  double p = fabs(1. / pars_array[4]);
  Vector3D direction(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  Vector3D mom = p * direction;
  Vector3D pos = (pars_array[0] * cos(pars_array[1])) * rot.col(0) +
                 (pars_array[0] * sin(pars_array[1])) * rot.col(1) + center;
  // constructor from parameter vector
  BoundParameters ataDisc_from_pars(tgContext, std::nullopt, pars, dSurface);
  consistencyCheck(ataDisc_from_pars, pos, mom, 1., 21., pars_array);
  // constructor from global parameters
  BoundParameters ataDisc_from_global(tgContext, std::nullopt, pos, mom, 1.,
                                      21., dSurface);
  consistencyCheck(ataDisc_from_global, pos, mom, 1., 21., pars_array);
  // constructor for neutral parameters
  NeutralBoundTrackParameters n_ataDisc_from_pars(tgContext, std::nullopt, pars,
                                                  dSurface);
  consistencyCheck(n_ataDisc_from_pars, pos, mom, 0., 21., pars_array);
  // constructor for neutral global parameters
  NeutralBoundTrackParameters n_ataDisc_from_global(tgContext, std::nullopt,
                                                    pos, mom, 21., dSurface);
  consistencyCheck(n_ataDisc_from_global, pos, mom, 0., 21., pars_array);

  // check shared ownership of same surface
  BOOST_CHECK_EQUAL(&ataDisc_from_pars.referenceSurface(), dSurface.get());
  BOOST_CHECK_EQUAL(&ataDisc_from_pars.referenceSurface(),
                    &ataDisc_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(&n_ataDisc_from_pars.referenceSurface(), dSurface.get());
  BOOST_CHECK_EQUAL(&n_ataDisc_from_pars.referenceSurface(),
                    &n_ataDisc_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(dSurface.use_count(), 5u);

  // check that the reference frame is the
  // rotation matrix of the surface
  const auto& dRotation =
      dSurface->transform(tgContext).matrix().block<3, 3>(0, 0);
  CHECK_CLOSE_REL(ataDisc_from_pars.referenceFrame(tgContext), dRotation, 1e-6);
}

/// @brief Unit test for parameters at a cylinder
///
BOOST_DATA_TEST_CASE(
    bound_to_cylinder_test,
    bdata::random((bdata::seed = 39810,
                   bdata::distribution =
                       std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 21221,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 62132,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((
            bdata::seed = 91683,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 39847,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 72615,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::xrange(100),
    x, y, z, a, b, c, index) {
  (void)index;

  Vector3D center{x, y, z};
  auto transform = std::make_shared<Transform3D>();
  transform->setIdentity();
  RotationMatrix3D rot;
  rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY()) *
        AngleAxis3D(c, Vector3D::UnitZ());
  transform->prerotate(rot);
  transform->pretranslate(center);

  auto bounds = std::make_shared<CylinderBounds>(100., 1200.);
  std::shared_ptr<const Surface> cSurface =
      Surface::makeShared<CylinderSurface>(transform, bounds);

  // now create parameters on this surface
  // rPhi, a, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {{125., 343., 0.45, 0.888, 0.001, 21.}};
  SingleBoundTrackParameters<ChargedPolicy>::ParVector_t pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  const double phi = pars_array[2];
  const double theta = pars_array[3];
  double p = fabs(1. / pars_array[4]);
  Vector3D direction(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  Vector3D mom = p * direction;

  // 3D position in local frame
  double r = bounds->get(CylinderBounds::eR);
  const double phi_l = pars_array[0] / r;
  Vector3D pos = (r * cos(phi_l)) * rot.col(0) + (r * sin(phi_l)) * rot.col(1) +
                 (pars_array[1]) * rot.col(2) + center;

  // constructor from parameter vector
  BoundParameters ataCylinder_from_pars(tgContext, std::nullopt, pars,
                                        cSurface);
  consistencyCheck(ataCylinder_from_pars, pos, mom, 1., 21., pars_array);
  // constructor from global parameters
  BoundParameters ataCylinder_from_global(tgContext, std::nullopt, pos, mom, 1.,
                                          21., cSurface);
  consistencyCheck(ataCylinder_from_global, pos, mom, 1., 21., pars_array);
  // constructor for neutral parameters
  NeutralBoundTrackParameters n_ataCylinder_from_pars(tgContext, std::nullopt,
                                                      pars, cSurface);
  consistencyCheck(n_ataCylinder_from_pars, pos, mom, 0., 21., pars_array);
  // constructor for neutral global parameters
  NeutralBoundTrackParameters n_ataCylinder_from_global(
      tgContext, std::nullopt, pos, mom, 21., cSurface);
  consistencyCheck(n_ataCylinder_from_global, pos, mom, 0., 21., pars_array);

  // check shared ownership of same surface
  BOOST_CHECK_EQUAL(&ataCylinder_from_pars.referenceSurface(), cSurface.get());
  BOOST_CHECK_EQUAL(&ataCylinder_from_pars.referenceSurface(),
                    &ataCylinder_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(&n_ataCylinder_from_pars.referenceSurface(),
                    cSurface.get());
  BOOST_CHECK_EQUAL(&n_ataCylinder_from_pars.referenceSurface(),
                    &n_ataCylinder_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(cSurface.use_count(), 5u);

  auto pPosition = ataCylinder_from_pars.position();
  // the reference frame is
  // transverse plane to the cylinder at the intersect
  Vector3D normal_at_intersect = cSurface->normal(tgContext, pPosition);
  Vector3D transverse_y = rot.col(2);
  Vector3D transverse_x = transverse_y.cross(normal_at_intersect);
  RotationMatrix3D refframe;
  refframe.col(0) = transverse_x;
  refframe.col(1) = transverse_y;
  refframe.col(2) = normal_at_intersect;
  // check if the manually constructed reference frame is the provided one
  CHECK_CLOSE_REL(ataCylinder_from_pars.referenceFrame(tgContext), refframe,
                  1e-6);
}

/// @brief Unit test for parameters at the perigee
///
BOOST_DATA_TEST_CASE(
    bound_to_perigee_test,
    bdata::random(
        (bdata::seed = 3980,
         bdata::distribution = std::uniform_real_distribution<>(-10., 10.))) ^
        bdata::random((bdata::seed = 2221,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-10., 10.))) ^
        bdata::random((bdata::seed = 2132,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-10., 10.))) ^
        bdata::random((
            bdata::seed = 9183,
            bdata::distribution = std::uniform_real_distribution<>(0., 0.05))) ^
        bdata::random((
            bdata::seed = 3947,
            bdata::distribution = std::uniform_real_distribution<>(0., 0.05))) ^
        bdata::xrange(100),
    x, y, z, a, b, index) {
  (void)index;
  Vector3D center{x, y, z};
  auto transform = std::make_shared<Transform3D>();
  transform->setIdentity();
  RotationMatrix3D rot;
  rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY());
  transform->prerotate(rot);
  transform->pretranslate(center);

  // the straw surface
  std::shared_ptr<const Surface> pSurface = Surface::makeShared<PerigeeSurface>(
      std::make_shared<const Transform3D>(*transform));

  // now create parameters on this surface
  // d0, z0, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {{-0.7321, 22.5, 0.45, 0.888, 0.001, 21.}};
  SingleBoundTrackParameters<ChargedPolicy>::ParVector_t pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  BoundParameters ataPerigee_from_pars(tgContext, std::nullopt, pars, pSurface);
  auto pos = ataPerigee_from_pars.position();
  auto mom = ataPerigee_from_pars.momentum();
  consistencyCheck(ataPerigee_from_pars, pos, mom, 1., 21., pars_array);
  // constructor from global parameters
  BoundParameters ataPerigee_from_global(tgContext, std::nullopt, pos, mom, 1.,
                                         21., pSurface);
  consistencyCheck(ataPerigee_from_global, pos, mom, 1., 21., pars_array);
  // constructor for neutral parameters
  NeutralBoundTrackParameters n_ataPerigee_from_pars(tgContext, std::nullopt,
                                                     pars, pSurface);
  consistencyCheck(n_ataPerigee_from_pars, pos, mom, 0., 21., pars_array);
  // constructor for neutral global parameters
  NeutralBoundTrackParameters n_ataPerigee_from_global(tgContext, std::nullopt,
                                                       pos, mom, 21., pSurface);
  consistencyCheck(n_ataPerigee_from_global, pos, mom, 0., 21., pars_array);

  // check shared ownership of same surface
  BOOST_CHECK_EQUAL(&ataPerigee_from_pars.referenceSurface(), pSurface.get());
  BOOST_CHECK_EQUAL(&ataPerigee_from_pars.referenceSurface(),
                    &ataPerigee_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(&n_ataPerigee_from_pars.referenceSurface(), pSurface.get());
  BOOST_CHECK_EQUAL(&n_ataPerigee_from_pars.referenceSurface(),
                    &n_ataPerigee_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(pSurface.use_count(), 5u);
}

/// @brief Unit test for parameters at a line
///
BOOST_DATA_TEST_CASE(
    bound_to_line_test,
    bdata::random((bdata::seed = 73980,
                   bdata::distribution =
                       std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 21221,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((bdata::seed = 62992,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-1000., 1000.))) ^
        bdata::random((
            bdata::seed = 900683,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 5439847,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::random((
            bdata::seed = 1972615,
            bdata::distribution = std::uniform_real_distribution<>(0., M_PI))) ^
        bdata::xrange(100),
    x, y, z, a, b, c, index) {
  using namespace Acts::UnitLiterals;

  (void)index;

  Vector3D center{x, y, z};
  auto transform = std::make_shared<Transform3D>();
  transform->setIdentity();
  RotationMatrix3D rot;
  rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY()) *
        AngleAxis3D(c, Vector3D::UnitZ());
  transform->prerotate(rot);
  transform->pretranslate(center);

  // the straw surface
  auto sSurface = Surface::makeShared<StrawSurface>(transform, 2_mm, 1_m);

  // now create parameters on this surface
  // r, z, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {{0.2321, 22.5, 0.45, 0.888, 0.001, 21.}};
  SingleBoundTrackParameters<ChargedPolicy>::ParVector_t pars;
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  // constructor from parameter vector
  BoundParameters ataLine_from_pars(tgContext, std::nullopt, pars, sSurface);
  auto pos = ataLine_from_pars.position();
  auto mom = ataLine_from_pars.momentum();
  consistencyCheck(ataLine_from_pars, pos, mom, 1., 21., pars_array);
  // constructor from global parameters
  BoundParameters ataLine_from_global(tgContext, std::nullopt, pos, mom, 1.,
                                      21., sSurface);
  consistencyCheck(ataLine_from_global, pos, mom, 1., 21., pars_array);
  // constructor for neutral parameters
  NeutralBoundTrackParameters n_ataLine_from_pars(tgContext, std::nullopt, pars,
                                                  sSurface);
  consistencyCheck(n_ataLine_from_pars, pos, mom, 0., 21., pars_array);
  // constructor for neutral global parameters
  NeutralBoundTrackParameters n_ataLine_from_global(tgContext, std::nullopt,
                                                    pos, mom, 21., sSurface);
  consistencyCheck(n_ataLine_from_global, pos, mom, 0., 21., pars_array);

  // check shared ownership of same surface
  BOOST_CHECK_EQUAL(&ataLine_from_pars.referenceSurface(), sSurface.get());
  BOOST_CHECK_EQUAL(&ataLine_from_pars.referenceSurface(),
                    &ataLine_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(&n_ataLine_from_pars.referenceSurface(), sSurface.get());
  BOOST_CHECK_EQUAL(&n_ataLine_from_pars.referenceSurface(),
                    &n_ataLine_from_global.referenceSurface());
  BOOST_CHECK_EQUAL(sSurface.use_count(), 5u);
}
}  // namespace Test
}  // namespace Acts
