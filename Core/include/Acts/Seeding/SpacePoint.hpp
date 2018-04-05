#ifndef SpacePoint_h
#define SpacePoint_h

#include <utility>

struct SpacePoint{
  public:
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  std::pair<int,int> m_clusterList = std::pair<int,int>(1,1);
  void setClusterList(int first, int second) {
    m_clusterList = std::pair<int,int>(first,second);
  }
  const std::pair<int,int> clusterList() const {
    return m_clusterList;
  }
  int surface;

  float& x(){
    return m_x;
  }
  float& y(){
    return m_y;
  }
  float& z(){
    return m_z;
  }
  float& r(){
    return m_r;
  }

};
#endif
