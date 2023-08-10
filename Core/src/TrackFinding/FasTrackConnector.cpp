#include "Acts/TrackFinding/FasTrackConnector.hpp"
#include <iostream>
#include <cstring>
#include <fstream>

namespace Acts {
  
FasTrackConnection::FasTrackConnection(unsigned int s, unsigned int d) : m_src(s), m_dst(d) { 

}

FasTrackConnector::FasTrackConnector(std::ifstream& inFile) {
  
  m_connMap.clear();

  int nLinks;



  inFile >> nLinks >> m_etaBin;


  for(int l=0;l<nLinks;l++) {

    unsigned int stage, lIdx, src, dst, nEntries;
    int height, width;

    inFile >> lIdx >> stage >> src >> dst >> height >> width >> nEntries;
    
    FasTrackConnection* pC = new FasTrackConnection(src, dst);
    
    int dummy;

    for(int i=0;i<height;i++) {
      for(int j=0;j<width;j++) inFile >> dummy;//pC->m_binTable[j+i*width];
    }

    int vol_id = src / 1000;

    if(vol_id == 13 || vol_id == 12 || vol_id == 14) {
      delete pC;
      continue;
    }

    vol_id = dst / 1000;
    
    if(vol_id == 13 || vol_id == 12 || vol_id == 14) {
      delete pC;
      continue;
    }

    std::map<int, std::vector<FasTrackConnection*> >::iterator it = m_connMap.find(stage);
    
    if(it == m_connMap.end()) {
      std::vector<FasTrackConnection*> v(1, pC);
      m_connMap.insert(std::make_pair(stage, v));
    } else (*it).second.push_back(pC);
  }
}

FasTrackConnector::~FasTrackConnector() {

  for(std::map<int, std::vector<FasTrackConnection*> >::iterator it = m_connMap.begin();it!=m_connMap.end();++it) {
    for(std::vector<FasTrackConnection*>::iterator cIt=(*it).second.begin();cIt!=(*it).second.end();++cIt) {
      delete (*cIt);
    }
    (*it).second.clear();
  }
  m_connMap.clear();

}
} 