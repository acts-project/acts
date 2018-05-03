#include "TrkSpacePoint/SpacePoint.h"

#include <boost/any.hpp>

void Acts::ATL_CovarianceTool::setCovariance(SPForSeed& sp, float zAlign, float rAlign){

  Trk::SpacePoint* atlasSp = boost::any_cast(sp->spacepoint);
  
  const InDet::SiCluster*           c  = static_cast<const InDet::SiCluster*>(atlasSp->clusterList().first);
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
