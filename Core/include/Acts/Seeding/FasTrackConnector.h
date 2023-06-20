
#include<fstream>
#include<vector>
#include<map>

struct FasTrackConnection {
public:
  FasTrackConnection(unsigned int, unsigned int);
  ~FasTrackConnection() {};

  unsigned int m_src, m_dst;
  std::vector<int> m_binTable;

};


class FasTrackConnector {
 public:
  FasTrackConnector(std::ifstream&);
  ~FasTrackConnector();

  float m_etaBin;

  std::map<int, std::vector<FasTrackConnection*> > m_connMap;

} ;


