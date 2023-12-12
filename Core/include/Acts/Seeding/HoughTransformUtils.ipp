
template <class PointType>
void Acts::HoughTransformUtils::HoughPlane::fill(
    const PointType& measurement, lineParametrisation<PointType> linePar,
    lineParametrisation<PointType> widthPar, idGetter<PointType> idGetter,
    layerGetter<PointType> layerGetter) {
  auto id = idGetter(measurement);
  auto layer = layerGetter(measurement);
  assert(layer < m_layerHists.size());

  for (int xBin = 0; xBin < m_cfg.nBinsY; xBin++) {
        auto x = m_binCentersX.at(xBin);
        // now evaluate the line equation provided by the user 
        double y = linePar(x,measurement);
        auto yBin = quant(m_cfg.yMin, m_cfg.yMax, m_cfg.nBinsY, y); 
        if (yBin >= m_cfg.nBinsY || yBin  < 0) continue; // out of bounds 
        // also get the width in the x coordinate 
        double dy = widthPar(x,measurement); 
        // fill the nominal bin
        m_layerHists[layer](xBin,yBin).first++; 
        m_layerHists[layer](xBin,yBin).second.insert(id); 
        // fill extra bins - to the left 
        for (int by = yBin; by > 0; --by){
            if (std::abs(m_binCentersY.at(by) - y) <= dy){
                m_layerHists[layer](xBin,by).first++; 
                m_layerHists[layer](xBin,by).second.insert(id); 
            }
        }
        // ... and to the right 
        for (int by = yBin; by < m_cfg.nBinsY ; ++by){
            if (std::abs(m_binCentersY.at(by) - y) <= dy){
                m_layerHists[layer](xBin,by).first++; 
                m_layerHists[layer](xBin,by).second.insert(id); 
            }
        }
  }
  m_mergedLayers = false; 
}
