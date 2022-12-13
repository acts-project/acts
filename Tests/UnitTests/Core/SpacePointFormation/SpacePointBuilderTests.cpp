// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/TestSpacePoint.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <limits>
#include <variant>
namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;

using TestMeasurement = Acts::BoundVariantMeasurement;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;
// Construct initial track parameters.
CurvilinearTrackParameters makeParameters(double phi, double theta, double p,
                                          double q) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // Let the particle start from the origin
  Vector4 mPos4(-3_m, 0., 0., 0.);
  return CurvilinearTrackParameters(mPos4, phi, theta, p, q, cov);
}

const Vector2 getLocalPos(const TestMeasurement* meas) {
  return std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        BoundVector par = expander * x.parameters();
        Vector2 local(par[BoundIndices::eBoundLoc0],
                      par[BoundIndices::eBoundLoc1]);
        return local;
      },
      *meas);
}

std::pair<Vector3, Vector3> stripEnds(
    const std::shared_ptr<const TrackingGeometry>& geo,
    const GeometryContext& gctx, const TestMeasurement* meas,
    const double stripFrac = 0.4) {
  const auto lpos = getLocalPos(meas);
  Vector3 globalFakeMom(1, 1, 1);
  const SourceLink& slink =
      std::visit([](const auto& x) { return x.sourceLink(); }, *meas);
  const auto geoId = slink.geometryId();
  const Surface* surface = geo->findSurface(geoId);

  const double stripLength = 40.;
  const double end1x = lpos[0] + stripLength * stripFrac;
  const double end1y = lpos[1];
  const double end2x = lpos[0] - stripLength * (1 - stripFrac);
  const double end2y = lpos[1];
  const Vector2 lpos1(end1x, end1y);
  const Vector2 lpos2(end2x, end2y);

  auto gPos1 = surface->localToGlobal(gctx, lpos1, globalFakeMom);
  auto gPos2 = surface->localToGlobal(gctx, lpos2, globalFakeMom);

  return std::make_pair(gPos1, gPos2);
}

// Create a test context
GeometryContext tgContext = GeometryContext();

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

// detector geometry
CubicTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();

// detector resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolution resStrip = {MeasurementType::eLoc01,
                                        {100_um, 100_um}};
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier().setVolume(2), resPixel},
    {GeometryIdentifier().setVolume(3).setLayer(2), resStrip},
    {GeometryIdentifier().setVolume(3).setLayer(4), resStrip},
    {GeometryIdentifier().setVolume(3).setLayer(6), resStrip},
    {GeometryIdentifier().setVolume(3).setLayer(8), resStrip},
};

// Construct a straight-line propagator.
static StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{std::move(geo)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator{cfg};
  Acts::StraightLineStepper stepper;
  return StraightPropagator(stepper, std::move(navigator));
}

// simulation propagator
const auto measPropagator = makeStraightPropagator(geometry);

std::default_random_engine rng(42);

BOOST_DATA_TEST_CASE(SpacePointBuilder_basic, bdata::xrange(1), index) {
  (void)index;

  double phi = 5._degree;
  double theta = 95._degree;
  double p = 50._GeV;
  double q = 1;

  Acts::Navigator navigator({
      geometry,
      true,  // sensitive
      true,  // material
      false  // passive
  });
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 2._T));
  ConstantFieldStepper stepper(std::move(field));

  ConstantFieldPropagator propagator(std::move(stepper), std::move(navigator));
  auto start = makeParameters(phi, theta, p, q);

  auto measurements =
      createMeasurements(propagator, geoCtx, magCtx, start, resolutions, rng);

  auto sourceLinks = measurements.sourceLinks;

  std::vector<const TestMeasurement*> frontMeasurements;
  std::vector<const TestMeasurement*> backMeasurements;
  std::vector<const TestMeasurement*> singleHitMeasurements;

  std::vector<const Vector3*> frontStripEnds;
  std::vector<const Vector3*> backStripEnds;

  for (auto& sl : sourceLinks) {
    const auto geoId = sl.geometryId();
    const auto volumeId = geoId.volume();
    if (volumeId == 2) {  // pixel type detector

      const TestMeasurement* meas = new TestMeasurement(
          makeMeasurement(SourceLink{sl}, sl.parameters, sl.covariance,
                          sl.indices[0], sl.indices[1]));
      singleHitMeasurements.emplace_back(meas);
    } else if (volumeId == 3) {  // strip type detector

      const auto layerId = geoId.layer();

      const TestMeasurement* meas = new TestMeasurement(
          makeMeasurement(SourceLink{sl}, sl.parameters, sl.covariance,
                          sl.indices[0], sl.indices[1]));

      if (layerId == 2 || layerId == 6) {
        frontMeasurements.emplace_back(meas);
      } else if (layerId == 4 || layerId == 8) {
        backMeasurements.emplace_back(meas);
      }
    }  // volume 3 (strip detector)
  }

  BOOST_CHECK_EQUAL(frontMeasurements.size(), 2);
  BOOST_CHECK_EQUAL(backMeasurements.size(), 2);

  Acts::Vector3 vertex = Vector3(-3_m, 0., 0.);

  auto spConstructor =
      [](const Acts::Vector3& pos, const Acts::Vector2& cov,
         boost::container::static_vector<Acts::SourceLink, 2> slinks)
      -> TestSpacePoint {
    return TestSpacePoint(pos, cov[0], cov[1], std::move(slinks));
  };

  auto spBuilderConfig = SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = geometry;
  spBuilderConfig.vertex = vertex;

  auto spBuilder =
      Acts::SpacePointBuilder<TestSpacePoint>(spBuilderConfig, spConstructor);

  // for cosmic  without vertex constraint, usePerpProj = true
  auto spBuilderConfig_perp = SpacePointBuilderConfig();
  spBuilderConfig_perp.trackingGeometry = geometry;
  spBuilderConfig_perp.vertex = vertex;
  spBuilderConfig_perp.usePerpProj = true;

  auto spBuilder_perp = Acts::SpacePointBuilder<TestSpacePoint>(
      spBuilderConfig_perp, spConstructor);

  TestSpacePointContainer spacePoints;
  TestSpacePointContainer spacePoints_extra;

  for (auto& meas : singleHitMeasurements) {
    std::vector<const TestMeasurement*> measVect;
    measVect.emplace_back(meas);

    SpacePointBuilderOptions spOpt;
    spBuilder.buildSpacePoint(geoCtx, measVect, spOpt,
                              std::back_inserter(spacePoints));
  }
  BOOST_CHECK_EQUAL(spacePoints.size(), 2);
  std::vector<std::pair<const TestMeasurement*, const TestMeasurement*>>
      measPairs;

  // strip SP building

  spBuilder.makeMeasurementPairs(tgContext, frontMeasurements, backMeasurements,
                                 measPairs);

  BOOST_CHECK_EQUAL(measPairs.size(), 2);

  for (auto& measPair : measPairs) {
    const auto meas1 = measPair.first;
    const auto meas2 = measPair.second;

    const std::pair<Vector3, Vector3> end1 = stripEnds(geometry, geoCtx, meas1);
    const std::pair<Vector3, Vector3> end2 = stripEnds(geometry, geoCtx, meas2);
    std::shared_ptr<const TestSpacePoint> spacePoint = nullptr;
    std::pair<const std::pair<Vector3, Vector3>,
              const std::pair<Vector3, Vector3>>
        strippair = std::make_pair(end1, end2);
    std::vector<const TestMeasurement*> measVect;
    measVect.emplace_back(meas1);
    measVect.emplace_back(meas2);

    SpacePointBuilderOptions spOpt{strippair};
    // nominal strip sp building
    spBuilder.buildSpacePoint(geoCtx, measVect, spOpt,
                              std::back_inserter(spacePoints));

    // sp building without vertex constraint
    spBuilder_perp.buildSpacePoint(geoCtx, measVect, spOpt,
                                   std::back_inserter(spacePoints));

    // put measurements slightly outside strips to test recovery
    const std::pair<Vector3, Vector3> end3 =
        stripEnds(geometry, geoCtx, meas1, 1.01);
    const std::pair<Vector3, Vector3> end4 =
        stripEnds(geometry, geoCtx, meas2, 1.02);
    // the other side of the strips
    const std::pair<Vector3, Vector3> end5 =
        stripEnds(geometry, geoCtx, meas1, -0.01);
    const std::pair<Vector3, Vector3> end6 =
        stripEnds(geometry, geoCtx, meas2, -0.02);

    auto spBuilderConfig_badStrips = SpacePointBuilderConfig();

    spBuilderConfig_badStrips.trackingGeometry = geometry;
    spBuilderConfig_badStrips.vertex = vertex;
    spBuilderConfig_badStrips.stripLengthTolerance = 0.0001;
    spBuilderConfig_badStrips.stripLengthGapTolerance = 50.;
    auto spBuilder_badStrips = Acts::SpacePointBuilder<TestSpacePoint>(
        spBuilderConfig_badStrips, spConstructor);
    // sp building with the recovery method
    SpacePointBuilderOptions spOpt_badStrips1{std::make_pair(end3, end4)};
    spBuilder_badStrips.buildSpacePoint(geoCtx, measVect, spOpt_badStrips1,
                                        std::back_inserter(spacePoints_extra));

    SpacePointBuilderOptions spOpt_badStrips2{std::make_pair(end5, end6)};
    spBuilder_badStrips.buildSpacePoint(geoCtx, measVect, spOpt_badStrips2,
                                        std::back_inserter(spacePoints_extra));
  }

  for (auto& sp : spacePoints) {
    std::cout << "space point (" << sp.x() << " " << sp.y() << " " << sp.z()
              << ") var (r,z): " << sp.varianceR() << " " << sp.varianceZ()
              << std::endl;
  }
  std::cout << "space points produced with bad strips:" << std::endl;
  for (auto& sp : spacePoints_extra) {
    std::cout << "space point (" << sp.x() << " " << sp.y() << " " << sp.z()
              << ") var (r,z): " << sp.varianceR() << " " << sp.varianceZ()
              << std::endl;
  }

  BOOST_CHECK_EQUAL(spacePoints.size(), 6);
}

}  // end of namespace Test
}  // namespace Acts
