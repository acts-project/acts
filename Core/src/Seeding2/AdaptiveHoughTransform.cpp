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

AccumulatorSection AccumulatorSection::bottomLeft(float xFraction,
                                                  float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction, m_xBegin,
                            m_yBegin, m_divisionLevel + 1, m_indices,
                            m_history);
}
AccumulatorSection AccumulatorSection::topLeft(float xFraction,
                                               float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction, m_xBegin,
                            m_yBegin + m_ySize - m_ySize * yFraction,
                            m_divisionLevel + 1, m_indices, m_history);
}
AccumulatorSection AccumulatorSection::topRight(float xFraction,
                                                float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction,
                            m_xBegin + m_xSize - m_xSize * xFraction,
                            m_yBegin + m_ySize - m_ySize * yFraction,
                            m_divisionLevel + 1, m_indices, m_history);
}
AccumulatorSection AccumulatorSection::bottomRight(float xFraction,
                                                   float yFraction) const {
  return AccumulatorSection(m_xSize * xFraction, m_ySize * yFraction,
                            m_xBegin + m_xSize - m_xSize * xFraction, m_yBegin,
                            m_divisionLevel + 1, m_indices, m_history);
}
AccumulatorSection AccumulatorSection::bottom(float yFraction) const {
  return bottomLeft(1.0, yFraction);
}
AccumulatorSection AccumulatorSection::top(float yFraction) const {
  return topLeft(1.0, yFraction);
}
AccumulatorSection AccumulatorSection::left(float xFraction) const {
  return bottomLeft(xFraction, 1.0);
}
AccumulatorSection AccumulatorSection::right(float xFraction) const {
  return bottomRight(xFraction, 1.0);
}

}