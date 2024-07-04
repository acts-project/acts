// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/container/flat_set.hpp>

template <typename SpacePoint>
Acts::PlanarSpacePointGrid<SpacePoint>
Acts::PlanarSpacePointGridCreator::createGrid(
    const Acts::PlanarSpacePointGridConfig& config,
    const Acts::PlanarSpacePointGridOptions& options) {
  if (!config.isInInternalUnits) {
    throw std::runtime_error(
        "PlanarSpacePointGridConfig not in ACTS internal units in "
        "PlanarSpacePointGridCreator::createGrid");
  }
  using AxisScalar = Acts::Vector3::Scalar;
  using namespace Acts::UnitLiterals;

  // vector that will store the edges of the bins of z
  std::vector<AxisScalar> zValues;

  if (config.zBinEdges.empty()) {
    throw std::runtime_error(
        "Empty zBinEdges ");
  } else {
    // Use the zBinEdges defined in the config
    zValues.reserve(config.zBinEdges.size());
    for (float bin : config.zBinEdges) {
      zValues.push_back(bin);
    }
  }

  // vector that will store the edges of the bins of x
  std::vector<AxisScalar> xValues;

  if (config.xBinEdges.empty()) {
    throw std::runtime_error(
        "Empty xBinEdges ");
  } else {
    // Use the xBinEdges defined in the config
    xValues.reserve(config.xBinEdges.size());
    for (float bin : config.xBinEdges) {
      xValues.push_back(bin);
    }
  }

  Axis<AxisType::Variable, AxisBoundaryType::Bound>
      xAxis(std::move(xValues));

  Axis<AxisType::Variable, AxisBoundaryType::Bound>
      zAxis(std::move(zValues));

  return Acts::PlanarSpacePointGrid<SpacePoint>(
      std::make_tuple(std::move(xAxis), std::move(zAxis)));
}

template <typename external_spacepoint_t,
          typename external_spacepoint_iterator_t, typename callable_t>
void Acts::PlanarSpacePointGridCreator::fillGrid(
    const Acts::SeedFinderConfigNA60<external_spacepoint_t>& config,
    const Acts::SeedFinderOptionsNA60& options,
    Acts::PlanarSpacePointGrid<external_spacepoint_t>& grid,
    external_spacepoint_iterator_t spBegin,
    external_spacepoint_iterator_t spEnd, callable_t&& toGlobal,
    Acts::Extent& rRangeSPExtent) {
  using iterated_value_t =
      typename std::iterator_traits<external_spacepoint_iterator_t>::value_type;
  using iterated_t = typename std::remove_const<
      typename std::remove_pointer<iterated_value_t>::type>::type;
  static_assert(std::is_pointer<iterated_value_t>::value,
                "Iterator must contain pointers to space points");
  static_assert(std::is_same<iterated_t, external_spacepoint_t>::value,
                "Iterator does not contain type this class was templated with");

  if (!config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderConfigNA60 not in ACTS internal units in BinnedSPGroup");
  }
  if (config.seedFilter == nullptr) {
    throw std::runtime_error("SeedFinderConfigNA60 has a null SeedFilter object");
  }
  if (!options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptionsNA60 not in ACTS internal units in BinnedSPGroup");
  }

  // get region of interest (or full detector if configured accordingly)
  float xMin = config.xMin;
  float xMax = config.xMax;
  float zMin = config.zMin;
  float zMax = config.zMax;

  // keep track of changed bins while sorting
  boost::container::flat_set<std::size_t> yBinsIndex;

  std::size_t counter = 0ul;
  for (external_spacepoint_iterator_t it = spBegin; it != spEnd;
       it++, ++counter) {
    if (*it == nullptr) {
      continue;
    }
    const external_spacepoint_t& sp = **it;
    const auto& [spPosition, variance, spTime] =
        toGlobal(sp, config.zAlign, config.rAlign, config.sigmaError);

    float spX = spPosition[0];
    float spY = spPosition[1];
    float spZ = spPosition[2];

    // store x,y,z values in extent
    rRangeSPExtent.extend({spX, spY, spZ});

    // remove SPs outside z and phi region
    if (spZ > zMax || spZ < zMin) {
      continue;
    }
    if (spX > xMax || spX < xMin) {
      continue;
    }

    auto isp = std::make_unique<InternalSpacePoint<external_spacepoint_t>>(
        counter, sp, spPosition, options.beamPos, variance, spTime);

    // fill rbins into grid
    Acts::Vector2 spLocation(isp->x(), isp->z());
    std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>&
        rbin = grid.atPosition(spLocation);
    rbin.push_back(std::move(isp));

    // keep track of the bins we modify so that we can later sort the SPs in
    // those bins only
    if (rbin.size() > 1) {
      yBinsIndex.insert(grid.globalBinFromPosition(spLocation));
    }
  }

  /// sort SPs in R for each filled bin
  for (auto& binIndex : yBinsIndex) {
    auto& rbin = grid.atPosition(binIndex);
    std::sort(rbin.begin(), rbin.end(), [](const auto& a, const auto& b) {
      return a->y() < b->y();
    });
  }
}
