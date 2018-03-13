#ifndef Seed_h
#define Seed_h
#include <list>

namespace Acts {

  template <typename SpacePoint>
  class Seed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    Seed();
    Seed(const SpacePoint*&,const SpacePoint*&, const double&);
    Seed
      (const SpacePoint*&,const SpacePoint*&,const SpacePoint*&,
       const double&);
    Seed(const Seed&);
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
 

  template <typename SpacePoint>
  inline const std::list<const SpacePoint*>& Seed<SpacePoint>::spacePoints() const 
    {
      return this->m_spacepoints;
    }

  template<typename SpacePoint>
  inline void Seed<SpacePoint>::erase() 
    {
      m_spacepoints.erase(m_spacepoints.begin(),m_spacepoints.end());
    }

  template<typename SpacePoint>
  inline void Seed<SpacePoint>::add(const SpacePoint*& p) 
    {
      m_spacepoints.push_back(p);
    }
  
  template<typename SpacePoint>
  inline void Seed<SpacePoint>::setZVertex(const double& z) 
    {
      m_zvertex = z;
    }

  template<typename SpacePoint>
  inline const double& Seed<SpacePoint>::zVertex() const 
    {
      return m_zvertex;
    }
  
  ///////////////////////////////////////////////////////////////////////////////
  // Constructors
  ///////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  Seed<SpacePoint>::Seed(){}

  template <typename SpacePoint>
  Seed<SpacePoint>::Seed(const SpacePoint*& b,const SpacePoint*& u, const double& vertex){
    m_zvertex = vertex;
    m_spacepoints.push_back(b);
    m_spacepoints.push_back(u);

  }

  template <typename SpacePoint>
  Seed<SpacePoint>::Seed(const SpacePoint*& b, const SpacePoint*& m, const SpacePoint*& u,
         const double& vertex){
    m_zvertex = vertex;
    m_spacepoints.push_back(b);
    m_spacepoints.push_back(m);
    m_spacepoints.push_back(u);

  }

  template <typename SpacePoint>
  Seed<SpacePoint>::~Seed(){
  }
} // end of name space

#endif  // Seed_h
