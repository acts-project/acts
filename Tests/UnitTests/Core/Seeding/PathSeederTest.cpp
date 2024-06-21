// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/detail/CuboidalDetectorHelper.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Seeding/PathSeeder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Utilities/Logger.hpp"

BOOST_AUTO_TEST_SUITE(PathSeeder)

using namespace Acts;
using namespace Acts::UnitLiterals;

using Axis = Acts::detail::EquidistantAxis;

GeometryContext gctx;

// Parameters for the geometry
const ActsScalar halfY = 10.;
const ActsScalar halfZ = 10.;
const ActsScalar deltaX = 10.;
const ActsScalar deltaYZ = 1.;

// Intersection finding to get the
// region of interest for seeding
class NoFieldIntersectionFinder {
 public:
  ActsScalar tol = 1e-4;

  std::vector<const Surface*> m_surfaces;

  // Find the intersections along the path
  // and return them in the order of the path
  // length
  std::vector<std::pair<GeometryIdentifier, Vector3>> operator()(
      const GeometryContext& geoCtx, const Vector3& position,
      const Vector3& direction,
      [[maybe_unused]] const ActsScalar& Pmag = 0) const {
    std::vector<std::pair<GeometryIdentifier, Vector3>> sIntersections;
    // Intersect the surfaces
    for (auto& surface : m_surfaces) {
      // Get the intersection
      auto sMultiIntersection = surface->intersect(geoCtx, position, direction,
                                                   BoundaryCheck(true), tol);

      // Take the closest
      auto closestForward = sMultiIntersection.closestForward();

      // Store if the intersection is reachable
      if (closestForward.status() == IntersectionStatus::reachable &&
          closestForward.pathLength() > 0.0) {
        sIntersections.push_back(
            {closestForward.object()->geometryId(), closestForward.position()});
        continue;
      }
    }
    return sIntersections;
  }
};

// Grid to store the source links for
// the seeding fast lookup
class SourceLinkGrid : public ISourceLinkGrid<Axis> {
 public:
  using eGrid = Grid<std::vector<SourceLink>, Axis, Axis>;

  /// Lookup table collection
  std::unordered_map<int, eGrid> m_lookupTables;

  /// Surface accessor
  SourceLinkSurfaceAccessor m_surfaceAccessor;

  void initialize(const GeometryContext& geoCtx,
                  std::vector<SourceLink> sourceLinks) override {
    // Lookup table for each layer
    std::unordered_map<int, eGrid> lookupTable;

    // Construct a binned grid for each layer
    for (int i : {14, 15, 16, 17}) {
      Axis xAxis(-halfY, halfY, 50);
      Axis yAxis(-halfZ, halfZ, 50);

      eGrid grid(std::make_tuple(xAxis, yAxis));
      lookupTable.insert({i, grid});
    }
    // Fill the grid with source links
    for (auto& sl : sourceLinks) {
      auto ssl = sl.get<Acts::detail::Test::TestSourceLink>();
      auto id = ssl.m_geometryId;

      // Grid works with global positions
      Acts::Vector3 globalPos = m_surfaceAccessor(sl)->localToGlobal(
          geoCtx, ssl.parameters, Acts::Vector3{0, 1, 0});

      auto bin = lookupTable.at(id.sensitive())
                     .localBinsFromPosition(
                         Acts::Vector2(globalPos.y(), globalPos.z()));
      lookupTable.at(id.sensitive()).atLocalBins(bin).push_back(sl);
    }

    m_lookupTables = lookupTable;
  };

  // Get the source link grid for a given geometry id
  eGrid getSourceLinkTable(
      const Acts::GeometryIdentifier& geoId) const override {
    return m_lookupTables.at(geoId.sensitive());
  }
};

// A simple path width provider to set
// the grid lookup boundaries around the
// intersection point
class PathWidthProvider {
 public:
  ActsScalar m_pathWidth = 0.1;

  std::pair<ActsScalar, ActsScalar> operator()(
      const GeometryContext& /*geoCtx*/,
      const GeometryIdentifier& /*geoId*/) const {
    return {m_pathWidth, m_pathWidth};
  }
};

// Calibrator to transform the source links
// to global coordinates
class SourceLinkCalibrator {
 public:
  SourceLinkSurfaceAccessor m_surfaceAccessor;

  Acts::Vector3 operator()(const GeometryContext& geoCtx,
                           const Acts::SourceLink& sourceLink) const {
    auto ssl = sourceLink.get<Acts::detail::Test::TestSourceLink>();
    auto res =
        m_surfaceAccessor(sourceLink)
            ->localToGlobal(geoCtx, ssl.parameters, Acts::Vector3{0, 1, 0});
    return res;
  }
};

// Estimator of the particle's energy,
// vertex, momentum direction at the IP
// and the direction at the first hit
class TrackEstimator {
 public:
  Vector3 ip;

  std::tuple<ActsScalar, Vector3, Vector3, Vector3> operator()(
      const GeometryContext& /*geoCtx*/, const Vector3& pivot) const {
    Vector3 direction = (pivot - ip).normalized();
    return {1._GeV, ip, direction, direction};
  };
};

// Construct a simple telescope detector
std::shared_ptr<Experimental::Detector> constructTelescopeDetector() {
  RotationMatrix3 rotation;
  double angle = 90_degree;
  Vector3 xPos(cos(angle), 0., sin(angle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-sin(angle), 0., cos(angle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Create a bunch of surface bounds
  auto surf0bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf1bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf2bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  auto surf3bounds = std::make_unique<RectangleBounds>(halfY, halfZ);

  // Create a bunch of surfaces
  auto transform0 = Transform3::Identity();
  auto surf0 = Surface::makeShared<PlaneSurface>(
      transform0 * Transform3(rotation), std::move(surf0bounds));
  auto geoId0 = GeometryIdentifier();
  geoId0.setSensitive(1);

  auto transform1 =
      Transform3::Identity() * Translation3(Vector3(2 * deltaX, 0, 0));
  auto surf1 = Surface::makeShared<PlaneSurface>(
      transform1 * Transform3(rotation), std::move(surf1bounds));
  auto geoId1 = GeometryIdentifier();
  geoId1.setSensitive(2);

  auto transform2 =
      Transform3::Identity() * Translation3(Vector3(4 * deltaX, 0, 0));
  auto surf2 = Surface::makeShared<PlaneSurface>(
      transform2 * Transform3(rotation), std::move(surf2bounds));
  auto geoId2 = GeometryIdentifier();
  geoId2.setSensitive(3);

  auto transform3 =
      Transform3::Identity() * Translation3(Vector3(6 * deltaX, 0, 0));
  auto surf3 = Surface::makeShared<PlaneSurface>(
      transform3 * Transform3(rotation), std::move(surf3bounds));
  auto geoId3 = GeometryIdentifier();
  geoId3.setSensitive(4);

  // Create a bunch of volume bounds
  auto vol0bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol1bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol2bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  auto vol3bounds = std::make_unique<CuboidVolumeBounds>(
      deltaX, halfY + deltaYZ, halfZ + deltaYZ);

  // Create a bunch of volumes
  auto vol0 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol0",
      transform0, std::move(vol0bounds), {surf0}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol1 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol1",
      transform1, std::move(vol1bounds), {surf1}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol2 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol2",
      transform2, std::move(vol2bounds), {surf2}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  auto vol3 = Experimental::DetectorVolumeFactory::construct(
      Experimental::defaultPortalAndSubPortalGenerator(), gctx, "vol3",
      transform3, std::move(vol3bounds), {surf3}, {},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes = {
      vol0, vol1, vol2, vol3};

  // Connect the volumes
  auto portalContainer = Experimental::detail::CuboidalDetectorHelper::connect(
      gctx, volumes, BinningValue::binX, {}, Logging::VERBOSE);

  // Make sure that the geometry ids are
  // independent of the potential Id generation
  // changes
  int id = 1;

  // Volume ids
  for (auto& volume : volumes) {
    volume->assignGeometryId(id);
    id++;
  }
  // Intervolume portal ids
  for (auto& volume : volumes) {
    for (auto& port : volume->portalPtrs()) {
      if (port->surface().geometryId() == 0) {
        port->surface().assignGeometryId(id);
        id++;
      }
    }
  }
  // Surface ids
  for (auto& surf : {surf0, surf1, surf2, surf3}) {
    auto geoId = GeometryIdentifier();
    geoId.setSensitive(id);
    surf->assignGeometryId(geoId);
    id++;
  }

  auto detector = Experimental::Detector::makeShared(
      "TelescopeDetector", volumes, Experimental::tryRootVolumes());

  return detector;
}

std::vector<SourceLink> createSourceLinks(
    const GeometryContext& geoCtx, const Experimental::Detector& detector) {
  NoFieldIntersectionFinder intersectionFinder;

  for (auto volume : detector.volumes()) {
    for (auto surface : volume->surfaces()) {
      intersectionFinder.m_surfaces.push_back(surface);
    }
  }

  Vector3 vertex(-5., 0., 0.);
  std::vector<ActsScalar> phis = {-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15};

  std::vector<SourceLink> sourceLinks;
  for (ActsScalar phi : phis) {
    Vector3 direction(cos(phi), sin(phi), 0.);

    auto intersections = intersectionFinder(geoCtx, vertex, direction);

    SquareMatrix2 cov = Acts::SquareMatrix2::Identity();

    for (auto& [id, refPoint] : intersections) {
      auto surf = *detector.sensitiveHierarchyMap().find(id);
      Acts::Vector2 val =
          surf->globalToLocal(geoCtx, refPoint, direction).value();

      Acts::detail::Test::TestSourceLink sourceLink(eBoundLoc0, eBoundLoc1, val,
                                                    cov, id, id.value());

      Acts::SourceLink sl{sourceLink};
      sourceLinks.push_back(sl);
    }
  }

  return sourceLinks;
}

BOOST_AUTO_TEST_CASE(PathSeederZeroField) {
  // Create detector
  auto detector = constructTelescopeDetector();

  // Create source links
  auto sourceLinks = createSourceLinks(gctx, *detector);

  // Prepare the PathSeeder
  auto pathSeederCfg = Acts::PathSeeder<Axis>::Config();

  // Grid to bin the source links
  Acts::detail::Test::TestSourceLink::SurfaceAccessor surfaceAccessor;
  surfaceAccessor.detector = detector.get();
  auto sourceLinkGrid = std::make_shared<SourceLinkGrid>();
  sourceLinkGrid->m_surfaceAccessor.connect<
      &Acts::detail::Test::TestSourceLink::SurfaceAccessor::operator()>(
      &surfaceAccessor);
  pathSeederCfg.sourceLinkGrid = sourceLinkGrid;

  // Estimator of the IP and first hit
  // parameters of the track
  TrackEstimator trackEstimator;
  trackEstimator.ip = Vector3(-5., 0., 0.);
  pathSeederCfg.trackEstimator.connect<&TrackEstimator::operator()>(
      &trackEstimator);

  // Transforms the source links to global coordinates
  SourceLinkCalibrator sourceLinkCalibrator;
  sourceLinkCalibrator.m_surfaceAccessor.connect<
      &Acts::detail::Test::TestSourceLink::SurfaceAccessor::operator()>(
      &surfaceAccessor);
  pathSeederCfg.sourceLinkCalibrator.connect<&SourceLinkCalibrator::operator()>(
      &sourceLinkCalibrator);

  // Intersection finder
  auto intersectionFinder = NoFieldIntersectionFinder();
  for (auto volume : detector->volumes()) {
    for (auto surface : volume->surfaces()) {
      intersectionFinder.m_surfaces.push_back(surface);
    }
  }
  pathSeederCfg.intersectionFinder
      .connect<&NoFieldIntersectionFinder::operator()>(&intersectionFinder);

  // Path width provider
  auto pathWidthProvider = PathWidthProvider();
  pathSeederCfg.pathWidthProvider.connect<&PathWidthProvider::operator()>(
      &pathWidthProvider);

  // First tracking layer
  Extent firstLayerExtent;
  firstLayerExtent.set(Acts::binX, -0.1, 0.1);
  firstLayerExtent.set(Acts::binY, -halfY - deltaYZ, halfY + deltaYZ);
  firstLayerExtent.set(Acts::binZ, -halfZ - deltaYZ, halfZ + deltaYZ);

  pathSeederCfg.firstLayerExtent = firstLayerExtent;

  // Create the PathSeeder
  Acts::PathSeeder<Axis> pathSeeder(pathSeederCfg);

  // Get the seeds
  auto seeds = pathSeeder.getSeeds(gctx, sourceLinks);

  // Check the seeds
  BOOST_CHECK_EQUAL(seeds.size(), 7);
  for (auto& seed : seeds) {
    BOOST_CHECK_EQUAL(seed.sourceLinks.size(), 4);
    BOOST_CHECK_EQUAL(seed.ipVertex, Vector3(-5., 0., 0.));
    BOOST_CHECK_EQUAL(seed.ipP, 1._GeV);
  }
}

BOOST_AUTO_TEST_SUITE_END()
