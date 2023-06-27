// #include "GNN_Geometry.h"
#include "GNN_DataStorage.pp"

#include<cmath>
#include<cstring>
#include<algorithm>

namespace Acts {
template <typename space_point_t>  
TrigFTF_GNN_Node<space_point_t>::TrigFTF_GNN_Node(const std::vector<space_point_t>& p, float minT = -100.0, float maxT = 100.0)  : m_sp(p), m_minCutOnTau(minT), m_maxCutOnTau(maxT) {
  m_in.clear();
  m_out.clear();
}
template <typename space_point_t>  
TrigFTF_GNN_Node<space_point_t>::~TrigFTF_GNN_Node() {

}

template <typename space_point_t>  
TrigFTF_GNN_EtaBin<space_point_t>::TrigFTF_GNN_EtaBin() {
  m_vn.clear();
}

template <typename space_point_t>  
TrigFTF_GNN_EtaBin<space_point_t>::~TrigFTF_GNN_EtaBin() {
  for(std::vector<TrigFTF_GNN_Node*>::iterator it = m_vn.begin();it!=m_vn.end();++it) {
    delete (*it);
  }
  m_vn.clear();
}

template <typename space_point_t> 
void TrigFTF_GNN_EtaBin<space_point_t>::sortByPhi() {
  std::sort(m_vn.begin(), m_vn.end(), TrigFTF_GNN_Node::CompareByPhi());
}

template <typename space_point_t> 
void TrigFTF_GNN_EtaBin<space_point_t>::generatePhiIndexing(float dphi) {

  for(unsigned int nIdx=0;nIdx<m_vn.size();nIdx++) {

    TrigFTF_GNN_Node<space_point_t>* pN = m_vn.at(nIdx);
    float phi = pN->m_sp.phi();
    if(phi <= M_PI-dphi) continue;
    
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi - 2*M_PI, nIdx));
    
  }

  for(unsigned int nIdx=0;nIdx<m_vn.size();nIdx++) {
    TrigFTF_GNN_Node<space_point_t>* pN = m_vn.at(nIdx);
    float phi = pN->m_sp.phi();
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi, nIdx));
  }

  for(unsigned int nIdx=0;nIdx<m_vn.size();nIdx++) {
    TrigFTF_GNN_Node<space_point_t>* pN = m_vn.at(nIdx);
    float phi = pN->m_sp.phi();
    if(phi >= -M_PI + dphi) break;
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi + 2*M_PI, nIdx));
  }
}

template <typename space_point_t>  
TrigFTF_GNN_DataStorage<space_point_t>::TrigFTF_GNN_DataStorage(const TrigFTF_GNN_Geometry<space_point_t>& g) : m_geo(g) {

  for(int k=0;k<g.num_bins();k++) {
    m_etaBins.push_back(TrigFTF_GNN_EtaBin());
  }
}

template <typename space_point_t>  
TrigFTF_GNN_DataStorage<space_point_t>::~TrigFTF_GNN_DataStorage() {

}
//input will be FTF SP so need to call .SP for other funcitons 
template <typename space_point_t>  
int TrigFTF_GNN_DataStorage<space_point_t>::addSpacePoint(const FTF_SP<space_point_t> sp, bool useML = false) { 

  const TrigFTF_GNN_Layer<space_point_t>* pL = m_geo.getTrigFTF_GNN_LayerByIndex<space_point_t>(sp.FTF_ID); //want ftf layer 

  if(pL==nullptr) return -1;

  int binIndex = pL->getEtaBin(sp.SP.z(), sp.SP.r());

  if(binIndex == -1) {
    return -2;
  }
  
  bool isBarrel = (pL->m_layer.m_type == 0);

  if(isBarrel) {
    float min_tau = -100.0;
    float max_tau =  100.0;
  // //cant do this bit yet as dont have cluster width, these loops just change the values on tau 
  //   // if (useML) {
  //   //   const Trk::SpacePoint* osp = sp.offlineSpacePoint();
  //   //   const InDet::PixelCluster* pCL = dynamic_cast<const InDet::PixelCluster*>(osp->clusterList().first);
  //   //   float cluster_width = pCL->width().widthPhiRZ().y();
  //   //   min_tau = 6.7*(cluster_width - 0.2);
  //   //   max_tau = 1.6 + 0.15/(cluster_width + 0.2) + 6.1*(cluster_width - 0.2);
  //   // }
    m_etaBins.at(binIndex).m_vn.push_back(new TrigFTF_GNN_Node(sp.SP, min_tau, max_tau));
  }
  else {
  //   // if (useML) {
  //   //   const Trk::SpacePoint* osp = sp.offlineSpacePoint();
  //   //   const InDet::PixelCluster* pCL = dynamic_cast<const InDet::PixelCluster*>(osp->clusterList().first);
  //   //   float cluster_width = pCL->width().widthPhiRZ().y();
  //   //   if(cluster_width > 0.2) return -3;
  //   // }
    m_etaBins.at(binIndex).m_vn.push_back(new TrigFTF_GNN_Node(sp.SP)); //add template here? 
  }

  return 0;
}



template <typename space_point_t> 
unsigned int TrigFTF_GNN_DataStorage<space_point_t>::numberOfNodes() const {

  unsigned int n=0;
  
  for(auto& b : m_etaBins) {
    n += b.m_vn.size();
  }
  return n;
}

template <typename space_point_t> 
void TrigFTF_GNN_DataStorage<space_point_t>::sortByPhi() {
  for(auto& b : m_etaBins) b.sortByPhi();
}

void TrigFTF_GNN_DataStorage<space_point_t>::generatePhiIndexing(float dphi) {
  for(auto& b : m_etaBins) b.generatePhiIndexing(dphi);
}


void TrigFTF_GNN_DataStorage::getConnectingNodes(std::vector<const TrigFTF_GNN_Node*>& vn) {
  
  vn.clear();
  vn.reserve(numberOfNodes());
  
  for(const auto& b : m_etaBins) {
    for(std::vector<TrigFTF_GNN_Node*>::const_iterator nIt = b.m_vn.begin();nIt!=b.m_vn.end();++nIt) {
      if((*nIt)->m_in.empty()) continue;
      if((*nIt)->m_out.empty()) continue;
      vn.push_back(*nIt);
    }
  }
}
} //end of ACTS namespace 