//5 classes geometry node, eta bin, data storage and edge 


#include<vector>
#include<map>
#include<algorithm>

#define MAX_SEG_PER_NODE 1000 //was 30
#define N_SEG_CONNS  6 //was 6

#include "Acts/Seeding/GNN_Geometry.h"


template <typename space_point_t>  
class TrigFTF_GNN_Node {
public:

  struct CompareByPhi {

    bool operator()(const TrigFTF_GNN_Node<space_point_t>*  n1, const TrigFTF_GNN_Node<space_point_t>*  n2) {
      return n1->m_sp.phi() < n2->m_sp.phi();
    }

  };
  //want constructor to take simspace point 
  TrigFTF_GNN_Node(const space_point_t&, float, float);  
  ~TrigFTF_GNN_Node();
  

 inline void addIn(int i) {
//     if(m_in.size()<MAX_SEG_PER_NODE) {
//       m_in.push_back(i);
//     }
  }

  inline void addOut(int i) {
//     if(m_out.size()<MAX_SEG_PER_NODE) {
//       m_out.push_back(i);
//     }
  }
  
  inline bool isConnector() const {
//     if(m_in.empty() || m_out.empty()) return false;
//     return true;
  }

  inline bool isFull() const {
//     if(m_in.size()==MAX_SEG_PER_NODE && m_out.size()==MAX_SEG_PER_NODE) return true;
//     else return false;
  }

  const std::vector<space_point_t>& m_sp;
  
  std::vector<unsigned int> m_in;//indices of the edges in the edge storage
  std::vector<unsigned int> m_out;
  float m_minCutOnTau, m_maxCutOnTau;

};

template <typename space_point_t>  
class TrigFTF_GNN_EtaBin {
public:
  TrigFTF_GNN_EtaBin();
  ~TrigFTF_GNN_EtaBin();

  void 
  hi();

  bool empty() const {
  //   return m_vn.empty();
  }
  
  void generatePhiIndexing(float);
  
  std::vector<TrigFTF_GNN_Node<space_point_t>*> m_vn;
  std::vector<std::pair<float, unsigned int> > m_vPhiNodes;

};

template <typename space_point_t>  
class TrigFTF_GNN_DataStorage {
public:
  TrigFTF_GNN_DataStorage(const TrigFTF_GNN_Geometry<space_point_t>& );
  ~TrigFTF_GNN_DataStorage();
 
  int addSpacePoint(const space_point_t&, bool); 

  unsigned int numberOfNodes() const;
  void getConnectingNodes(std::vector<const TrigFTF_GNN_Node<space_point_t>*>&);
  void sortByPhi();
  void generatePhiIndexing(float);


  const TrigFTF_GNN_EtaBin<space_point_t>& getEtaBin(int idx) const {
  //   if(idx >= static_cast<int>(m_etaBins.size())) idx = idx-1;
  //   return m_etaBins.at(idx);
  }

protected:

  const TrigFTF_GNN_Geometry<space_point_t>&  m_geo;

  std::vector<TrigFTF_GNN_EtaBin<space_point_t>> m_etaBins; 

};

template <typename space_point_t>  
class TrigFTF_GNN_Edge {
public:

  struct CompareLevel {
  // public:
  //   bool operator()(const TrigFTF_GNN_Edge* pS1, const TrigFTF_GNN_Edge* pS2) {
  //     return pS1->m_level > pS2->m_level;
  //   }
  };

 TrigFTF_GNN_Edge() : m_n1(nullptr), m_n2(nullptr), m_level(-1), m_next(-1), m_nNei(0) {};

 TrigFTF_GNN_Edge(const TrigFTF_GNN_Edge<space_point_t>&  e) : m_n1(e.m_n1), m_n2(e.m_n2) {};

  inline void initialize(TrigFTF_GNN_Node<space_point_t>* n1, TrigFTF_GNN_Node<space_point_t>* n2) {
    m_n1 = n1; 
    m_n2 = n2;
    m_level = 1;
    m_next = 1;
    m_nNei = 0;
  }


  TrigFTF_GNN_Node<space_point_t>* m_n1{nullptr};
  TrigFTF_GNN_Node<space_point_t>* m_n2{nullptr};
  
  signed char m_level{-1}, m_next{-1};

  unsigned char m_nNei{0};
  float m_p[4]{};
  
  unsigned int m_vNei[N_SEG_CONNS]{};//global indices of the connected edges

};

