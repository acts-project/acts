#pragma once

#include <vector>

namespace Acts {
template <typename SpacePoint>
class Seed
{

  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

public:
  Seed();
  Seed(const SpacePoint*, const SpacePoint*, const SpacePoint*, float);
  Seed(const Seed&) = default;
  Seed&
  operator=(const Seed&);
  virtual ~Seed();
  const std::vector<const SpacePoint*>&
  sp() const
  {
    return m_spacepoints;
  }
  double
  z() const
  {
    return m_zvertex;
  }

private:
  std::vector<const SpacePoint*> m_spacepoints;
  float                          m_zvertex;
};

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
Seed<SpacePoint>::Seed()
{
}

template <typename SpacePoint>
Seed<SpacePoint>::Seed(const SpacePoint* b,
                       const SpacePoint* m,
                       const SpacePoint* u,
                       float             vertex)
{
  m_zvertex = vertex;
  m_spacepoints.push_back(b);
  m_spacepoints.push_back(m);
  m_spacepoints.push_back(u);
}

template <typename SpacePoint>
Seed<SpacePoint>::~Seed()
{
}
}  // end of Acts namespace
