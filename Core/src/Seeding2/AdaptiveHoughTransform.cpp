#include "Acts/Seeding2/AdaptiveHoughTransform.hpp"

namespace Acts {
AccumulatorSection::AccumulatorSection(float xs, float ys, float xBegin,
                                       float yBegin, int div,
                                       const std::vector<unsigned> &indices,
                                       const std::vector<float> &history)
    : m_xSize(xs),
      m_ySize(ys),
      m_xBegin(xBegin),
      m_yBegin(yBegin),
      m_divisionLevel(div),
      m_indices(indices),
      m_history(history) {}

void AccumulatorSection::updateDimensions(float xs, float ys, float xBegin,
                                          float yBegin) {
  m_xSize = xs;
  m_ySize = ys;
  m_xBegin = xBegin;
  m_yBegin = yBegin;
}

void AccumulatorSection::expand(float xs, float ys) {
  m_xBegin = m_xBegin + 0.5f * m_xSize - m_xSize * xs * 0.5f;
  m_yBegin = m_yBegin + 0.5f * m_ySize - m_ySize * ys * 0.5f;
  m_xSize *= xs;
  m_ySize *= ys;
}

AccumulatorSection AccumulatorSection::bottomLeft(bool copyIndices) const {
  return AccumulatorSection(m_xSize * 0.5, m_ySize * 0.5, m_xBegin,
                            m_yBegin, m_divisionLevel + 1,
                            (copyIndices ? m_indices : std::vector<unsigned>()), m_history);
}

AccumulatorSection AccumulatorSection::topLeft(bool copyIndices) const {
  return AccumulatorSection(m_xSize * 0.5, m_ySize * 0.5, m_xBegin,
                            m_yBegin + m_ySize - m_ySize * 0.5,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);
}

AccumulatorSection AccumulatorSection::topRight(bool copyIndices) const {
  return AccumulatorSection(m_xSize * 0.5, m_ySize * 0.5,
                            m_xBegin + m_xSize - m_xSize * 0.5,
                            m_yBegin + m_ySize - m_ySize * 0.5,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);
}

AccumulatorSection AccumulatorSection::bottomRight(bool copyIndices) const {
  return AccumulatorSection(m_xSize * 0.5, m_ySize * 0.5,
                            m_xBegin + m_xSize - m_xSize * 0.5, m_yBegin,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);
}

AccumulatorSection AccumulatorSection::bottom(bool copyIndices) const {
  return AccumulatorSection(m_xSize, m_ySize * 0.5,
                            m_xBegin, m_yBegin,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);

}

AccumulatorSection AccumulatorSection::top(bool copyIndices) const {
  return AccumulatorSection(m_xSize, m_ySize * 0.5, m_xBegin,
                            m_yBegin + m_ySize * 0.5,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);

}

AccumulatorSection AccumulatorSection::left(bool copyIndices) const {
  return AccumulatorSection(m_xSize * 0.5, m_ySize, m_xBegin,
                            m_yBegin,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);

}

AccumulatorSection AccumulatorSection::right(bool copyIndices) const {
    return AccumulatorSection(m_xSize * 0.5, m_ySize,
                            m_xBegin + m_xSize * 0.5,
                            m_yBegin,
                            m_divisionLevel + 1, (copyIndices ? m_indices : std::vector<unsigned>()),
                            m_history);

}

}  // namespace Acts