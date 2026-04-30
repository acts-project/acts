// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/HoughAccumulatorSection.hpp"

namespace Acts {
HoughAccumulatorSection::HoughAccumulatorSection(
    float xs, float ys, float xBegin, float yBegin, int div,
    const std::vector<unsigned> &indices, const std::vector<float> &history)
    : m_xSize(xs),
      m_ySize(ys),
      m_xBegin(xBegin),
      m_yBegin(yBegin),
      m_divisionLevel(div),
      m_indices(indices),
      m_history(history) {}

void HoughAccumulatorSection::updateDimensions(float xs, float ys, float xBegin,
                                               float yBegin) {
  m_xSize = xs;
  m_ySize = ys;
  m_xBegin = xBegin;
  m_yBegin = yBegin;
}

void HoughAccumulatorSection::expand(float xs, float ys) {
  m_xBegin = m_xBegin + 0.5f * m_xSize - m_xSize * xs * 0.5f;
  m_yBegin = m_yBegin + 0.5f * m_ySize - m_ySize * ys * 0.5f;
  m_xSize *= xs;
  m_ySize *= ys;
}

HoughAccumulatorSection HoughAccumulatorSection::bottomLeft(
    bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize * 0.5, m_ySize * 0.5, m_xBegin, m_yBegin, m_divisionLevel + 1,
      (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::topLeft(
    bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize * 0.5, m_ySize * 0.5, m_xBegin,
      m_yBegin + m_ySize - m_ySize * 0.5, m_divisionLevel + 1,
      (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::topRight(
    bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize * 0.5, m_ySize * 0.5, m_xBegin + m_xSize - m_xSize * 0.5,
      m_yBegin + m_ySize - m_ySize * 0.5, m_divisionLevel + 1,
      (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::bottomRight(
    bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize * 0.5, m_ySize * 0.5, m_xBegin + m_xSize - m_xSize * 0.5,
      m_yBegin, m_divisionLevel + 1,
      (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::bottom(
    bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize, m_ySize * 0.5, m_xBegin, m_yBegin, m_divisionLevel + 1,
      (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::top(bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize, m_ySize * 0.5, m_xBegin, m_yBegin + m_ySize * 0.5,
      m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
      m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::left(bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize * 0.5, m_ySize, m_xBegin, m_yBegin, m_divisionLevel + 1,
      (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

HoughAccumulatorSection HoughAccumulatorSection::right(bool copyIndices) const {
  return HoughAccumulatorSection(
      m_xSize * 0.5, m_ySize, m_xBegin + m_xSize * 0.5, m_yBegin,
      m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
      m_history);
}

}  // namespace Acts