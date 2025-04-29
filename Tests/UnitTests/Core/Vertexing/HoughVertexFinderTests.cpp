// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HoughVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace Acts::Test {

/// @brief SpacePoint definition to be used for the unit tests. Implements all the relevant methods.
struct SpacePoint4HVFT {
  SpacePoint4HVFT(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}
  double m_x;
  double m_y;
  double m_z;
  double x() const { return m_x; }
  double y() const { return m_y; }
  double z() const { return m_z; }
  double r() const { return std::sqrt(m_x * m_x + m_y * m_y); }
};

/// @brief Provides random double number between $from and $to
/// @param gen random number generator
/// @param from lower threshold
/// @param to upper threshold
/// @return random number in [from,to)
double getRndDouble(std::mt19937& gen, double from, double to) {
  std::uniform_real_distribution<double> dist(from, to);
  return dist(gen);
}

/// @brief Provides random integer number between $from and $to
/// @param gen random number generator
/// @param from lower threshold
/// @param to upper threshold
/// @return random number in [from,to]
int getRndInt(std::mt19937& gen, int from, int to) {
  std::uniform_int_distribution<int> dist(from, to);
  return dist(gen);
}

/// @brief Unit test for HoughVertexFinder. Compare the result to the easy-to-calculate expected result
BOOST_AUTO_TEST_CASE(hough_vertex_finder_small_test) {
  Acts::HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = 1000;
  houghVtxCfg.minAbsEta = 0.3;
  houghVtxCfg.maxAbsEta = 3.0;
  houghVtxCfg.minHits = 3;
  houghVtxCfg.fillNeighbours = 0;
  houghVtxCfg.absEtaRanges = std::vector<double>({3.0});
  houghVtxCfg.absEtaFractions = std::vector<double>({1.0});
  houghVtxCfg.rangeIterZ =
      std::vector<double>({100.05 * Acts::UnitConstants::mm});
  houghVtxCfg.nBinsZIterZ = std::vector<unsigned int>({2001});
  houghVtxCfg.nBinsCotThetaIterZ = std::vector<unsigned int>({1000});
  houghVtxCfg.binsCotThetaDecrease = 1.0;
  houghVtxCfg.peakWidth = 1;
  houghVtxCfg.defVtxPosition[0] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[1] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[2] = 0. * Acts::UnitConstants::mm;

  Acts::HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(
      std::move(houghVtxCfg));

  double vtxX = 0., vtxY = 0., vtxZ = 20.;

  std::vector<std::vector<double>> positions = {
      {10., 0., 25.},   {20., 0., 30.},   {30., 0., 35.},      // track 1
      {0., 5., 19.},    {0., 10., 18.},   {0, 15., 17.},       // track 2
      {-6., -4., 22.5}, {-12., -8., 25.}, {-18., -12., 27.5},  // track 3
      {-8., 2., 23.5},  {-16., 4., 27.},  {-24., 6., 30.5}};   // track 4

  std::vector<SpacePoint4HVFT> inputSpacepoints;
  for (auto pos : positions) {
    inputSpacepoints.emplace_back(pos[0], pos[1], pos[2]);
  }

  auto vtx = houghVertexFinder.find(inputSpacepoints);

  bool vtxFound = false;
  if (vtx.ok()) {
    // check if the found vertex has a compatible position
    if (std::abs(vtxX - vtx.value()[0]) < 1e-3 &&
        std::abs(vtxY - vtx.value()[1]) < 1e-3 &&
        std::abs(vtxZ - vtx.value()[2]) < 1e-3) {
      vtxFound = true;
    }
  }

  BOOST_CHECK(vtxFound);
}

/// @brief Unit test for HoughVertexFinder. Generates real-looking sets of the spacepoints, then finds a vertex, and then verifies the reconstructed vertex is actually near the original one
BOOST_AUTO_TEST_CASE(hough_vertex_finder_full_test) {
  HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = 1000;
  houghVtxCfg.minAbsEta = 0.3;
  houghVtxCfg.maxAbsEta = 3.0;
  houghVtxCfg.minHits = 3;
  houghVtxCfg.fillNeighbours = 0;
  houghVtxCfg.absEtaRanges = std::vector<double>({3.0});
  houghVtxCfg.absEtaFractions = std::vector<double>({1.0});
  houghVtxCfg.rangeIterZ =
      std::vector<double>({100.05 * Acts::UnitConstants::mm});
  houghVtxCfg.nBinsZIterZ = std::vector<unsigned int>({2001});
  houghVtxCfg.nBinsCotThetaIterZ = std::vector<unsigned int>({1000});
  houghVtxCfg.binsCotThetaDecrease = 1.35;
  houghVtxCfg.peakWidth = 3;
  houghVtxCfg.defVtxPosition[0] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[1] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[2] = 0. * Acts::UnitConstants::mm;

  HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(houghVtxCfg);

  std::mt19937 gen(299792458);

  int vtxFound = 0;
  int nEvents = 20;
  for (int event = 0; event < nEvents; event++) {
    double vtxX = getRndDouble(gen, -0.1, 0.1);
    double vtxY = getRndDouble(gen, -0.1, 0.1);
    double vtxZ = getRndDouble(gen, -50., 50.);

    std::vector<SpacePoint4HVFT> inputSpacepoints;

    // make straight lines originating from the given vertex
    int nTracks = getRndInt(gen, 200, 1000);
    for (int track = 0; track < nTracks; ++track) {
      // direction of the track
      double dirX = getRndDouble(gen, -1., 1.);
      double dirY = getRndDouble(gen, -1., 1.);
      double dirZ = getRndDouble(gen, -1., 1.);
      // use upper or lower intersection?
      int part = getRndInt(gen, 0, 1) * 2 - 1;

      for (int rIndx = 1; rIndx <= 3; rIndx += 1) {
        double sgn = std::copysign(1., dirY);
        double dirR2 = dirX * dirX + dirY * dirY;
        double D = vtxX * (vtxY + dirY) - vtxY * (vtxX + dirX);
        // add some smearing to the layers
        // layers are (9-11), (19-21), and (29-31)
        double r = rIndx * 10 + getRndDouble(gen, -1., 1.);
        // intersection of the layer and the straigh line
        double x1 =
            (D * dirY + part * sgn * dirX * std::sqrt(r * r * dirR2 - D * D)) /
            dirR2;
        double y1 = (-D * dirX +
                     part * std::abs(dirY) * std::sqrt(r * r * dirR2 - D * D)) /
                    dirR2;
        // how many units from the vertex to the intersection
        double zDist = std::abs((x1 - vtxX) / dirX);
        // use the same amount of units for distance in Z
        inputSpacepoints.emplace_back(x1, y1, zDist * dirZ + vtxZ);
      }
    }

    auto vtx = houghVertexFinder.find(inputSpacepoints);

    if (vtx.ok()) {
      // check if the found vertex has a compatible position
      if (std::abs(vtxX - vtx.value()[0]) < 0.2 &&
          std::abs(vtxY - vtx.value()[1]) < 0.2 &&
          std::abs(vtxZ - vtx.value()[2]) < 0.2) {
        ++vtxFound;
      }
    }
  }

  BOOST_CHECK_EQUAL(vtxFound, nEvents);
}

/// @brief Unit test for HoughVertexFinder. Provides no input spacepoints
BOOST_AUTO_TEST_CASE(hough_vertex_finder_empty_test) {
  Acts::HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  Acts::HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(
      std::move(houghVtxCfg));

  // no input spacepoints
  std::vector<SpacePoint4HVFT> inputSpacepoints;

  auto vtx = houghVertexFinder.find(inputSpacepoints);

  bool vtxFound = false;
  if (vtx.ok()) {
    vtxFound = true;
  }

  BOOST_CHECK(!vtxFound);
}

/// @brief Unit test for HoughVertexFinder. Does not provides enough spacepoints
BOOST_AUTO_TEST_CASE(hough_vertex_finder_insufficient_test) {
  Acts::HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = 1000;
  houghVtxCfg.minAbsEta = 0.3;
  houghVtxCfg.maxAbsEta = 3.0;
  houghVtxCfg.minHits = 3;  // requires 3 spacepoints per track
  houghVtxCfg.fillNeighbours = 0;
  houghVtxCfg.absEtaRanges = std::vector<double>({3.0});
  houghVtxCfg.absEtaFractions = std::vector<double>({1.0});
  houghVtxCfg.rangeIterZ =
      std::vector<double>({100.05 * Acts::UnitConstants::mm});
  houghVtxCfg.nBinsZIterZ = std::vector<unsigned int>({2001});
  houghVtxCfg.nBinsCotThetaIterZ = std::vector<unsigned int>({1000});
  houghVtxCfg.binsCotThetaDecrease = 1.0;
  houghVtxCfg.peakWidth = 1;
  houghVtxCfg.defVtxPosition[0] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[1] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[2] = 0. * Acts::UnitConstants::mm;

  Acts::HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(
      std::move(houghVtxCfg));

  // only 2 spacepoints per track provided
  std::vector<std::vector<double>> positions = {
      {10., 0., 25.},   {20., 0., 30.},     // track 1
      {0., 5., 19.},    {0., 10., 18.},     // track 2
      {-6., -4., 22.5}, {-12., -8., 25.}};  // track 3

  std::vector<SpacePoint4HVFT> inputSpacepoints;
  for (auto pos : positions) {
    inputSpacepoints.emplace_back(pos[0], pos[1], pos[2]);
  }

  auto vtx = houghVertexFinder.find(inputSpacepoints);

  bool vtxFound = false;
  if (vtx.ok()) {
    vtxFound = true;
  }

  BOOST_CHECK(!vtxFound);
}

/// @brief Unit test for HoughVertexFinder. Misconfigured #1
BOOST_AUTO_TEST_CASE(hough_vertex_finder_misconfig1_test) {
  Acts::HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = 1000;
  houghVtxCfg.minAbsEta = 0.3;
  houghVtxCfg.maxAbsEta = 3.0;
  houghVtxCfg.minHits = 3;
  houghVtxCfg.fillNeighbours = 0;
  houghVtxCfg.absEtaRanges = std::vector<double>({3.0, 4.0});  // misconfigured
  houghVtxCfg.absEtaFractions = std::vector<double>({1.0});
  houghVtxCfg.rangeIterZ =
      std::vector<double>({100.05 * Acts::UnitConstants::mm});
  houghVtxCfg.nBinsZIterZ = std::vector<unsigned int>({2001});
  houghVtxCfg.nBinsCotThetaIterZ = std::vector<unsigned int>({1000});
  houghVtxCfg.binsCotThetaDecrease = 1.0;
  houghVtxCfg.peakWidth = 1;
  houghVtxCfg.defVtxPosition[0] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[1] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[2] = 0. * Acts::UnitConstants::mm;

  BOOST_CHECK_THROW(Acts::HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(
                        std::move(houghVtxCfg)),
                    std::invalid_argument);
}

/// @brief Unit test for HoughVertexFinder. Misconfigured #2
BOOST_AUTO_TEST_CASE(hough_vertex_finder_misconfig2_test) {
  Acts::HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = 1000;
  houghVtxCfg.minAbsEta = 0.3;
  houghVtxCfg.maxAbsEta = 3.0;
  houghVtxCfg.minHits = 3;
  houghVtxCfg.fillNeighbours = 0;
  houghVtxCfg.absEtaRanges = std::vector<double>({3.0});
  houghVtxCfg.absEtaFractions = std::vector<double>({1.0});
  houghVtxCfg.rangeIterZ = std::vector<double>(
      {100.05 * Acts::UnitConstants::mm, 100.05 * Acts::UnitConstants::mm});
  houghVtxCfg.nBinsZIterZ = std::vector<unsigned int>({2001});  // misconfigured
  houghVtxCfg.nBinsCotThetaIterZ = std::vector<unsigned int>({1000, 1000});
  houghVtxCfg.binsCotThetaDecrease = 1.0;
  houghVtxCfg.peakWidth = 1;
  houghVtxCfg.defVtxPosition[0] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[1] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[2] = 0. * Acts::UnitConstants::mm;

  BOOST_CHECK_THROW(Acts::HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(
                        std::move(houghVtxCfg)),
                    std::invalid_argument);
}

/// @brief Unit test for HoughVertexFinder. Misconfigured #3
BOOST_AUTO_TEST_CASE(hough_vertex_finder_misconfig3_test) {
  Acts::HoughVertexFinder<SpacePoint4HVFT>::Config houghVtxCfg;
  houghVtxCfg.targetSPs = 1000;
  houghVtxCfg.minAbsEta = 0.3;
  houghVtxCfg.maxAbsEta = 3.0;
  houghVtxCfg.minHits = 3;
  houghVtxCfg.fillNeighbours = 0;
  houghVtxCfg.absEtaRanges = std::vector<double>({3.0});
  houghVtxCfg.absEtaFractions = std::vector<double>({1.0});
  houghVtxCfg.rangeIterZ = std::vector<double>(
      {100.05 * Acts::UnitConstants::mm, 100.05 * Acts::UnitConstants::mm});
  houghVtxCfg.nBinsZIterZ = std::vector<unsigned int>({2001, 2001});
  houghVtxCfg.nBinsCotThetaIterZ =
      std::vector<unsigned int>({1000});  // misconfigured
  houghVtxCfg.binsCotThetaDecrease = 1.0;
  houghVtxCfg.peakWidth = 1;
  houghVtxCfg.defVtxPosition[0] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[1] = 0. * Acts::UnitConstants::mm;
  houghVtxCfg.defVtxPosition[2] = 0. * Acts::UnitConstants::mm;

  BOOST_CHECK_THROW(Acts::HoughVertexFinder<SpacePoint4HVFT> houghVertexFinder(
                        std::move(houghVtxCfg)),
                    std::invalid_argument);
}

}  // namespace Acts::Test
