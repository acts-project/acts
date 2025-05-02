// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/Segmentation.hpp"

#include <iostream>

ActsFatras::CartesianSegmentation::CartesianSegmentation(
    const Acts::ProtoAxis& xAxis, const Acts::ProtoAxis& yAxis)
    : ISegmentation(), m_xAxis(xAxis), m_yAxis(yAxis) {}

std::array<std::size_t, 2> ActsFatras::CartesianSegmentation::bin(
    const Acts::Vector2& pos) const {
  return {m_xAxis.getAxis().getBin(pos.x()) - 1u,
          m_yAxis.getAxis().getBin(pos.y()) - 1u};
}

Acts::Vector2 ActsFatras::CartesianSegmentation::position(
    const std::array<std::size_t, 2>& bin) const {
  return {m_xAxis.getAxis().getBinCenter(bin[0] + 1u),
          m_yAxis.getAxis().getBinCenter(bin[1] + 1u)};
}

std::vector<ActsFatras::ChannelStep>
ActsFatras::CartesianSegmentation::channelSteps(
    const Acts::Vector2& start, const Acts::Vector2& end) const {
  // Prepare the channel steps
  std::vector<ActsFatras::ChannelStep> cSteps;

  auto bStart = bin(start);
  auto bEnd = bin(end);

  // Fast single channel exit
  if (bStart == bEnd) {
    return cSteps;
  }

  Acts::Vector2 segment2d = (end - start);

  // The're change in the x direction
  if (bStart[0] != bEnd[0]) {
    double k = segment2d.y() / segment2d.x();
    double d = start.y() - k * start.x();

    const auto& xBoundaries = m_xAxis.getAxis().getBinEdges();
    std::vector<double> xbBounds = {
        xBoundaries.begin() + std::min(bStart[0], bEnd[0]) + 1,
        xBoundaries.begin() + std::max(bStart[0], bEnd[0]) + 1};
    for (const auto& x : xbBounds) {
      cSteps.emplace_back(ChannelStep{
          {(bStart[0] < bEnd[0] ? 1 : -1), 0}, {x, k * x + d}, start});
    }
  }

  // The lines channel segment lines along y
  if (bStart[1] != bEnd[1]) {
    double k = segment2d.x() / segment2d.y();
    double d = start.x() - k * start.y();
    const auto& yBoundaries = m_yAxis.getAxis().getBinEdges();
    std::vector<double> ybBounds = {
        yBoundaries.begin() + std::min(bStart[1], bEnd[1]) + 1,
        yBoundaries.begin() + std::max(bStart[1], bEnd[1]) + 1};
    for (const auto& y : ybBounds) {
      cSteps.emplace_back(ChannelStep{
          {0, (bStart[1] < bEnd[1] ? 1 : -1)}, {k * y + d, y}, start});
    }
  }

  return cSteps;
}
