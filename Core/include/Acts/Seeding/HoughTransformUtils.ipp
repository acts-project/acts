
template <class PointType>
void Acts::HoughTransformUtils::HoughPlane::fill(
    const PointType& measurement, const houghAxisRanges& axisRanges,
    lineParametrisation<PointType> linePar,
    lineParametrisation<PointType> widthPar, idType identifier, unsigned layer,
    yieldType weight) {
  for (int xBin = 0; xBin < m_cfg.nBinsX; xBin++) {
    auto x = binCenter(axisRanges.xMin, axisRanges.xMax, m_cfg.nBinsX, xBin);
    // now evaluate the line equation provided by the user
    coordType y = linePar(x, measurement);
    coordType dy = widthPar(x, measurement);
    auto yBinDown =
        binIndex(axisRanges.yMin, axisRanges.yMax, m_cfg.nBinsY, y - dy);
    auto yBinUp =
        binIndex(axisRanges.yMin, axisRanges.yMax, m_cfg.nBinsY, y + dy);
    for (int yBin = yBinDown; yBin <= yBinUp; ++yBin) {
      if (yBin >= m_cfg.nBinsY || yBin < 0)
        continue;  // out of bounds
      fillBin(xBin, yBin, identifier, layer, weight);
    }
  }
}
