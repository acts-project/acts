//TODO: update to C++17 style 
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace Acts {

struct FasTrackConnection {
public:
  FasTrackConnection(unsigned int, unsigned int);
  ~FasTrackConnection(){};

  unsigned int m_src, m_dst;
  std::vector<int> m_binTable;
};

class FasTrackConnector {
public:

  struct LayerGroup {
  LayerGroup(unsigned int l1Key, const std::vector<const Acts::FasTrackConnection*>& v) : m_dst(l1Key), m_sources(v) {};

    unsigned int m_dst;//the target layer of the group
    std::vector<const Acts::FasTrackConnection*> m_sources;//the source layers of the group
  };

  FasTrackConnector(std::ifstream &);

  ~FasTrackConnector();

  float m_etaBin;

  std::map<int, std::vector<struct LayerGroup> > m_layerGroups;
  std::map<int, std::vector<Acts::FasTrackConnection *>> m_connMap;
};

} // namespace Acts