#ifndef Seed_h
#define Seed_h
#include <list>
#include "ACTS/Seeding/SpacePoint.hpp"

namespace Acts {
namespace Seeding{
  
  class Seed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    Seed();
    Seed
      (const SpacePoint*&,const SpacePoint*&,const SpacePoint*&,
       const double&);
    Seed(const Seed&) = default;
    Seed& operator = (const Seed&);
    virtual ~Seed();
    void                                     erase();
    void                                     add(const SpacePoint*&);
    void                                     setZVertex(const double&);
    const std::list<const SpacePoint*>& spacePoints() const;
    const double&                            zVertex    () const;
    
    /////////////////////////////////////////////////////////////////////////////////
    // Protected data members
    /////////////////////////////////////////////////////////////////////////////////
    
  protected:
    
    std::list<const SpacePoint*> m_spacepoints;
    double                            m_zvertex    ;  
  };

  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////
 

  
  inline const std::list<const SpacePoint*>& Seed::spacePoints() const 
    {
      return this->m_spacepoints;
    }

  inline void Seed::erase() 
    {
      m_spacepoints.erase(m_spacepoints.begin(),m_spacepoints.end());
    }

  inline void Seed::add(const SpacePoint*& p) 
    {
      m_spacepoints.push_back(p);
    }
  
  inline void Seed::setZVertex(const double& z) 
    {
      m_zvertex = z;
    }

  inline const double& Seed::zVertex() const 
    {
      return m_zvertex;
    }
  
  ///////////////////////////////////////////////////////////////////////////////
  // Constructors
  ///////////////////////////////////////////////////////////////////////////////

  
  Seed::Seed(){}

  
  Seed::Seed(const SpacePoint*& b, const SpacePoint*& m, const SpacePoint*& u,
         const double& vertex){
    m_zvertex = vertex;
    m_spacepoints.push_back(b);
    m_spacepoints.push_back(m);
    m_spacepoints.push_back(u);

  }

  
  Seed::~Seed(){
  }
} // end of Seeding namespace
} // end of Acts namespace

#endif  // Seed_h
