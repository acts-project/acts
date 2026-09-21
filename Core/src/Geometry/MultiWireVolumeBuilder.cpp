// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/MultiWireVolumeBuilder.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <map>

// A struct to avoid having surfaces that belong to
// the same layer into different bins along Y because of floating point
// comparison
namespace {

struct PrecisionDelimiter {
  PrecisionDelimiter(const double limit = 1e-3) : m_limit{limit} {}
  bool operator()(const double a, const double b) const {
    return a + m_limit < b;
  }

 private:
  double m_limit{1e-3};
};

}  // namespace

namespace Acts {

MultiWireVolumeBuilder::MultiWireVolumeBuilder(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_config(config), m_logger(std::move(logger)) {
  if (m_config.mlSurfaces.empty()) {
    throw std::invalid_argument(
        "MultiWireStructureBuilder: No surfaces are given");
  }

  // check that straw tube surfaces hav been passed tot he builder
  const bool allStraw =
      std::ranges::all_of(m_config.mlSurfaces, [](const auto& s) {
        return s->bounds().type() == SurfaceBounds::BoundsType::eLine;
      });
  if (!allStraw) {
    throw std::invalid_argument(
        "MultiWireVolumeBuilder: all surfaces must have LineBounds");
  }
}

std::unique_ptr<TrackingVolume> MultiWireVolumeBuilder::buildVolume() const {
  // Create the tracking volume

  ACTS_VERBOSE("Building a tracking volume with name "
               << m_config.name << " ,translation"
               << toString(m_config.transform.translation())
               << " and number of surfaces " << m_config.mlSurfaces.size());

  auto boundsType = m_config.bounds ? m_config.bounds->type()
                                    : VolumeBounds::BoundsType::eOther;
  if (!(boundsType == VolumeBounds::BoundsType::eTrapezoid ||
        boundsType == VolumeBounds::BoundsType::eCuboid ||
        boundsType == VolumeBounds::BoundsType::eDiamond)) {
    throw std::invalid_argument(
        "MultiWireStructureBuilder: Only trapezoid cuboid or diamond bounds "
        "are "
        "supported");
  }

  std::unique_ptr<TrackingVolume> trackingVolume{};

  if (m_config.alignablePlacement == nullptr) {
    trackingVolume = std::make_unique<TrackingVolume>(
        m_config.transform, m_config.bounds, m_config.name);
  } else {
    trackingVolume = std::make_unique<TrackingVolume>(
        *m_config.alignablePlacement, m_config.bounds, m_config.name);
  }

  // Add the surfaces to the tracking volume
  for (auto& surface : m_config.mlSurfaces) {
    trackingVolume->addSurface(surface);
  }

  return trackingVolume;
}

std::tuple<AxisSpec, AxisSpec, std::vector<double>>
MultiWireVolumeBuilder::deriveGridParameters(
    const GeometryContext& gctx) const {
  // The global to local transformation to get to the volume's frame
  const Transform3 globalToLoc =
      (m_config.alignablePlacement == nullptr)
          ? m_config.transform.inverse()
          : m_config.alignablePlacement->globalToLocalTransform(gctx);

  if (m_config.binning.size() != 2u) {
    throw std::invalid_argument(
        "MultiWireVolumeBuilder: exactly two binning directions required");
  }
  const AxisDirection dirA = std::get<0>(m_config.binning.at(0));
  const AxisDirection dirB = std::get<0>(m_config.binning.at(1));
  const AxisDirection shiftDir = m_config.shiftDirection;

  if (shiftDir != dirA && shiftDir != dirB) {
    throw std::invalid_argument(
        "MultiWireVolumeBuilder: shiftDirection must be one of the binning "
        "directions");
  }

  // The layer direction
  const AxisDirection layerDir = (shiftDir == dirA) ? dirB : dirA;

  // project every tube onto both directions, grouped by layer
  std::map<double, std::set<double>, PrecisionDelimiter> coordsPerLayer;
  for (const auto& surf : m_config.mlSurfaces) {
    const Vector3 cLocal = globalToLoc * surf->center(gctx);
    const long long layerKey = VectorHelpers::cast(cLocal, layerDir);
    const long long shiftKey = VectorHelpers::cast(cLocal, shiftDir);
    coordsPerLayer[layerKey].insert(shiftKey);
  }

  // min gap between adjacent tubes , find the pitch
  const auto minGap = [&](const auto& sortedSet) -> double {
    double g = std::numeric_limits<double>::max();
    for (auto it = std::next(sortedSet.begin()); it != sortedSet.end(); ++it) {
      g = std::min(g, *it - *std::prev(it));
    }
    return g;
  };

  // shift axis: from the combined tube positions along shiftDir
  // Use the 1st plane as the reference plane
  const auto& firstPlane = coordsPerLayer.begin()->second;
  const double shiftPitch = minGap(firstPlane);
  const double shiftLow = *firstPlane.begin() - 0.5 * shiftPitch;
  const double shiftHigh = *firstPlane.rbegin() + 0.5 * shiftPitch;
  const auto nShiftBins = static_cast<std::size_t>(
      std::lround((shiftHigh - shiftLow) / shiftPitch));

  // layerr axis: from the layer keys (one distinct value per layer)
  //     assumes uniform layer spacing (see check below).
  double layerPitch = std::numeric_limits<double>::max();
  for (auto it = std::next(coordsPerLayer.begin()); it != coordsPerLayer.end();
       ++it) {
    layerPitch = std::min(layerPitch, (it->first - std::prev(it)->first));
  }
  const double layerLow = coordsPerLayer.begin()->first - 0.5 * layerPitch;
  const double layerHigh = coordsPerLayer.rbegin()->first + 0.5 * layerPitch;
  const auto nLayerBins = static_cast<std::size_t>(
      std::lround((layerHigh - layerLow) / layerPitch));

  // per-layer staggering correction, indexed by layer-axis bin (if enabled)
  // Use as reference the first tube of the first layer and then calculate the
  // offsets along the shift direction
  std::vector<double> layerShifts(nLayerBins, 0.0);
  if (m_config.correctOffsets) {
    const double refPos = *coordsPerLayer.begin()->second.begin();
    for (const auto& [lk, coords] : coordsPerLayer) {
      const double firstCoord = *coords.begin();
      // how far this layer is staggered w.r.t layer 0
      double d = refPos - firstCoord;
      // calculate the index of the vector which corresponds to the layer index
      const auto index =
          static_cast<std::size_t>(std::floor((lk - layerLow) / layerPitch));
      layerShifts.at(index) = d;
    }
  }
  // assemble the axis specs
  AxisSpec shiftAxis = AxisSpec::Equidistant(nShiftBins, shiftLow, shiftHigh,
                                             AxisBoundaryType::Bound, shiftDir);
  AxisSpec layerAxis = AxisSpec::Equidistant(nLayerBins, layerLow, layerHigh,
                                             AxisBoundaryType::Bound, layerDir);

  return {std::move(shiftAxis), std::move(layerAxis), std::move(layerShifts)};
}

std::unique_ptr<Acts::NavigationPolicyFactory>
MultiWireVolumeBuilder::createNavigationPolicyFactory(
    const GeometryContext& gctx) const {
  if (m_config.binning.size() != 2u) {
    throw std::invalid_argument(
        "MultiWireStructureBuilder: Invalid binning provided");
  }
  auto [axisDirectionA, expansionA] = m_config.binning.at(0);
  auto [axisDirectionB, expansionB] = m_config.binning.at(1);

  if (axisDirectionA == axisDirectionB) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: The axis directions need to be different for "
        "a two-dimensional grid");
  }

  // check if the direction along which the tubes are shifted is also consistent
  if (m_config.shiftDirection != axisDirectionA &&
      m_config.shiftDirection != axisDirectionB) {
    throw std::invalid_argument(
        "MultiWireVolumeBuilder: shiftDirection must be one of the two "
        "binning axis directions");
  }

  const auto [shiftAxisSpec, layerAxisSpec, layerShifts] =
      deriveGridParameters(gctx);

  if (shiftAxisSpec.isDeferred() || layerAxisSpec.isDeferred()) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need a fully specified range");
  }

  if (!shiftAxisSpec.isEquidistant() || !layerAxisSpec.isEquidistant()) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need to be equidistant");
  }

  const auto& shiftParams = shiftAxisSpec.asEquidistant();
  const auto& layerParams = layerAxisSpec.asEquidistant();

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisShift(
      *shiftParams.min, *shiftParams.max, shiftParams.nBins);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisLayer(
      *layerParams.min, *layerParams.max, layerParams.nBins);

  Grid<std::vector<std::size_t>, decltype(axisShift), decltype(axisLayer)> grid(
      axisShift, axisLayer);
  ACTS_VERBOSE(
      "MultiWireVolumeBuilder: Assign Multi-layer Navigation Policy with Grid "
      "axis: "
      << axisShift << "," << axisLayer);
  // The indexed grid to be filled from the navigation policy
  // The first axis direction corresponds to the shift direction (the direction
  // the tubes are aligned) The second axis direction corresponds to the
  // direction from one layer to another
  const auto* placement = m_config.alignablePlacement;
  auto indexedGrid =
      placement == nullptr
          ? IndexGrid<decltype(grid)>{std::move(grid),
                                      {*shiftAxisSpec.direction(),
                                       *layerAxisSpec.direction()},
                                      m_config.transform.inverse()}
          : IndexGrid<decltype(grid)>{
                std::move(grid),
                {*shiftAxisSpec.direction(), *layerAxisSpec.direction()},
                [placement](const GeometryContext& gctx2) -> const Transform3& {
                  return placement->globalToLocalTransform(gctx2);
                }};

  TryAllNavigationPolicy::Config tryAllConfig;
  tryAllConfig.portals = true;
  tryAllConfig.sensitives = false;

  const std::size_t shiftExp =
      (m_config.shiftDirection == axisDirectionA) ? expansionA : expansionB;
  const std::size_t layerExp =
      (m_config.shiftDirection == axisDirectionA) ? expansionB : expansionA;

  // build grid, IndexGrid with casts {shiftDir, layerDir}
  MultiLayerNavigationPolicy::Config navConfig;
  navConfig.binExpansion = {shiftExp, layerExp};
  navConfig.layerOffsets = layerShifts;

  // Create the navigation policy factory
  std::unique_ptr<NavigationPolicyFactory> factory =
      NavigationPolicyFactory{}
          .add<TryAllNavigationPolicy>(tryAllConfig)
          .add<MultiLayerNavigationPolicy>(navConfig, indexedGrid)
          .asUniquePtr();

  return factory;
}

}  // namespace Acts
