

void Acts::ATL_CovarianceTool::setCovariance(SPForSeed& sp, float zAlign, float rAlign){
  const InDet::SiCluster*           c  = static_cast<const InDet::SiCluster*>(sp->spacepoint->clusterList().first);
  const InDetDD::SiDetectorElement* de = c ->detectorElement();

  if( de->isPixel() ) {
    
    const Amg::MatrixX& v =  c->localCovariance();
    float f22 = float(v(1,1) );
    float wid = float(c->width().z());
    float cov = wid*wid*.08333; if(cov < f22) cov = f22;
    if(de->isBarrel()) {m_covz = 9.*cov; m_covr = .06;}
    else               {m_covr = 9.*cov; m_covz = .06;}
  }
  else {
    const Amg::MatrixX& v = sp->spacepoint->localCovariance();
    float f22 = float(v(1,1));
    if(de->isBarrel()) {m_covz = 8.*f22; m_covr = .1;} 
    else               {m_covr = 8.*f22; m_covz = .1;} 
  }
}
