// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/container/flat_set.hpp>

template <typename SpacePoint>
Acts::CylindricalSpacePointGrid<SpacePoint>
Acts::CylindricalSpacePointGridCreator::createGrid(
    const Acts::CylindricalSpacePointGridConfig& config,
    const Acts::CylindricalSpacePointGridOptions& options) {
  if (!config.isInInternalUnits) {
    throw std::runtime_error(
        "CylindricalSpacePointGridConfig not in ACTS internal units in "
        "CylindricalSpacePointGridCreator::createGrid");
  }
  if (!options.isInInternalUnits) {
    throw std::runtime_error(
        "CylindricalSpacePointGridOptions not in ACTS internal units in "
        "CylindricalSpacePointGridCreator::createGrid");
  }
  using AxisScalar = Acts::Vector3::Scalar;
  using namespace Acts::UnitLiterals;

  int phiBins = 0;
  // for no magnetic field, create 100 phi-bins
  if (options.bFieldInZ == 0) {
    phiBins = 100;
  } else {
    // calculate circle intersections of helix and max detector radius
    float minHelixRadius =
        config.minPt /
        (1_T * 1e6 *
         options.bFieldInZ);  // in mm -> R[mm] =pT[GeV] / (3·10−4×B[T])
                              // = pT[MeV] / (300 *Bz[kT])

    // sanity check: if yOuter takes the square root of a negative number
    if (minHelixRadius < config.rMax / 2) {
      throw std::domain_error(
          "The value of minHelixRadius cannot be smaller than rMax / 2. Please "
          "check the configuration of bFieldInZ and minPt");
    }

    float maxR2 = config.rMax * config.rMax;
    float xOuter = maxR2 / (2 * minHelixRadius);
    float yOuter = std::sqrt(maxR2 - xOuter * xOuter);
    float outerAngle = std::atan(xOuter / yOuter);
    // intersection of helix and max detector radius minus maximum R distance
    // from middle SP to top SP
    float innerAngle = 0;
    float rMin = config.rMax;
    if (config.rMax > config.deltaRMax) {
      rMin = config.rMax - config.deltaRMax;
      float innerCircleR2 =
          (config.rMax - config.deltaRMax) * (config.rMax - config.deltaRMax);
      float xInner = innerCircleR2 / (2 * minHelixRadius);
      float yInner = std::sqrt(innerCircleR2 - xInner * xInner);
      innerAngle = std::atan(xInner / yInner);
    }

    // evaluating the azimutal deflection including the maximum impact parameter
    float deltaAngleWithMaxD0 =
        std::abs(std::asin(config.impactMax / (rMin)) -
                 std::asin(config.impactMax / config.rMax));

    // evaluating delta Phi based on the inner and outer angle, and the azimutal
    // deflection including the maximum impact parameter
    // Divide by config.phiBinDeflectionCoverage since we combine
    // config.phiBinDeflectionCoverage number of consecutive phi bins in the
    // seed making step. So each individual bin should cover
    // 1/config.phiBinDeflectionCoverage of the maximum expected azimutal
    // deflection
    float deltaPhi = (outerAngle - innerAngle + deltaAngleWithMaxD0) /
                     config.phiBinDeflectionCoverage;

    // sanity check: if the delta phi is equal to or less than zero, we'll be
    // creating an infinite or a negative number of bins, which would be bad!
    if (deltaPhi <= 0.f) {
      throw std::domain_error(
          "Delta phi value is equal to or less than zero, leading to an "
          "impossible number of bins (negative or infinite)");
    }

    // divide 2pi by angle delta to get number of phi-bins
    // size is always 2pi even for regions of interest
    phiBins = static_cast<int>(std::ceil(2 * M_PI / deltaPhi));
    // need to scale the number of phi bins accordingly to the number of
    // consecutive phi bins in the seed making step.
    // Each individual bin should be approximately a fraction (depending on this
    // number) of the maximum expected azimutal deflection.

    // set protection for large number of bins, by default it is large
    if (phiBins > config.maxPhiBins) {
      phiBins = config.maxPhiBins;
    }
  }

  Acts::detail::Axis<detail::AxisType::Equidistant,
                     detail::AxisBoundaryType::Closed>
      phiAxis(config.phiMin, config.phiMax, phiBins);

  // vector that will store the edges of the bins of z
  std::vector<AxisScalar> zValues;

  // If zBinEdges is not defined, calculate the edges as zMin + bin * zBinSize
  if (config.zBinEdges.empty()) {
    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering
    float zBinSize = config.cotThetaMax * config.deltaRMax;
    float zBins =
        std::max(1.f, std::floor((config.zMax - config.zMin) / zBinSize));

    zValues.reserve(static_cast<int>(zBins));
    for (int bin = 0; bin <= static_cast<int>(zBins); bin++) {
      AxisScalar edge =
          config.zMin + bin * ((config.zMax - config.zMin) / zBins);
      zValues.push_back(edge);
    }

  } else {
    // Use the zBinEdges defined in the config
    zValues.reserve(config.zBinEdges.size());
    for (float bin : config.zBinEdges) {
      zValues.push_back(bin);
    }
  }

  detail::Axis<detail::AxisType::Variable, detail::AxisBoundaryType::Bound>
      zAxis(std::move(zValues));
  return Acts::CylindricalSpacePointGrid<SpacePoint>(
      std::make_tuple(std::move(phiAxis), std::move(zAxis)));
}

template <typename external_spacepoint_t,
          typename external_spacepoint_iterator_t, typename callable_t>
void Acts::CylindricalSpacePointGridCreator::fillGrid(
    const Acts::SeedFinderConfig<external_spacepoint_t>& config,
    const Acts::SeedFinderOptions& options,
    Acts::CylindricalSpacePointGrid<external_spacepoint_t>& grid,
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
        "SeedFinderConfig not in ACTS internal units in BinnedSPGroup");
  }
  if (config.seedFilter == nullptr) {
    throw std::runtime_error("SeedFinderConfig has a null SeedFilter object");
  }
  if (!options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptions not in ACTS internal units in BinnedSPGroup");
  }

  // get region of interest (or full detector if configured accordingly)
  float phiMin = config.phiMin;
  float phiMax = config.phiMax;
  float zMin = config.zMin;
  float zMax = config.zMax;

  // sort by radius
  // add magnitude of beamPos to rMax to avoid excluding measurements
  // create number of bins equal to number of millimeters rMax
  // (worst case minR: configured minR + 1mm)
  // binSizeR allows to increase or reduce numRBins if needed
  std::size_t numRBins = static_cast<std::size_t>(
      (config.rMax + options.beamPos.norm()) / config.binSizeR);

  // keep track of changed bins while sorting
  boost::container::flat_set<std::size_t> rBinsIndex;

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
    float spPhi = std::atan2(spY, spX);
    if (spPhi > phiMax || spPhi < phiMin) {
      continue;
    }

    auto isp = std::make_unique<InternalSpacePoint<external_spacepoint_t>>(
        counter, sp, spPosition, options.beamPos, variance, spTime);
    // calculate r-Bin index and protect against overflow (underflow not
    // possible)
    std::size_t rIndex =
        static_cast<std::size_t>(isp->radius() / config.binSizeR);
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins) {
      continue;
    }

    // fill rbins into grid
    Acts::Vector2 spLocation(isp->phi(), isp->z());
    std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>&
        rbin = grid.atPosition(spLocation);
    rbin.push_back(std::move(isp));

    // keep track of the bins we modify so that we can later sort the SPs in
    // those bins only
    if (rbin.size() > 1) {
      rBinsIndex.insert(grid.globalBinFromPosition(spLocation));
    }
  }

  /// sort SPs in R for each filled bin
  for (auto& binIndex : rBinsIndex) {
    auto& rbin = grid.atPosition(binIndex);
    std::sort(rbin.begin(), rbin.end(), [](const auto& a, const auto& b) {
      return a->radius() < b->radius();
    });
  }
}
