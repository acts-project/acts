#ifndef InternalSeed_h
#define InternalSeed_h

#include "ACTS/Seeding/SPForSeed.hpp"
#include "ACTS/Seeding/Seed.hpp"

//TODO: still unchanged from ATLAS
namespace Acts {
namespace Seeding {
  class InternalSeed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    InternalSeed();
    InternalSeed(SPForSeed*&,SPForSeed*&,SPForSeed*&,float);
    InternalSeed(const InternalSeed&);
    virtual ~InternalSeed();
    InternalSeed& operator  = (const InternalSeed&);

    SPForSeed* spacepoint0() {return m_s0;}
    SPForSeed* spacepoint1() {return m_s1;}
    SPForSeed* spacepoint2() {return m_s2;}
    const float&             z() const {return m_z ;}
    const float&       quality() const {return m_q ;}

    void set(SPForSeed*&,SPForSeed*&,SPForSeed*&,float);

    void setQuality(float, bool);

    bool set3(Seed&);

    void set2(Seed&);

  protected:
    
    SPForSeed* m_s0;
    SPForSeed* m_s1;
    SPForSeed* m_s2;
    float                m_z ;
    float                m_q ;
  };


  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////

  inline InternalSeed::InternalSeed ()
    {
      m_s0 = 0 ;
      m_s1 = 0 ;
      m_s2 = 0 ;
      m_z  = 0.;
      m_q  = 0.;
    }

  inline InternalSeed& InternalSeed::operator = 
    (const InternalSeed& sp) 
    {
      if(&sp!=this) {

    m_z   = sp.m_z ;
    m_q   = sp.m_q ;
    m_s0  = sp.m_s0;
    m_s1  = sp.m_s1;
    m_s2  = sp.m_s2;
      }
      return(*this);
    }
 
  inline InternalSeed::InternalSeed
    (SPForSeed*& s0,SPForSeed*& s1,SPForSeed*& s2,float z)
    {
      set(s0,s1,s2,z); m_q = 0.;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Copy constructor
  /////////////////////////////////////////////////////////////////////////////////

  inline InternalSeed::InternalSeed (const InternalSeed& sp): m_s0(sp.m_s0),m_s1(sp.m_s1),m_s2(sp.m_s2)
    {
      *this = sp;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Destructor
  /////////////////////////////////////////////////////////////////////////////////

  inline InternalSeed::~InternalSeed() 
  {
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Set 
  /////////////////////////////////////////////////////////////////////////////////

  inline void InternalSeed::set
    (SPForSeed*& s0,SPForSeed*& s1,SPForSeed*& s2,float z)
    {
      m_z   = z ;
      m_s0  = s0;
      m_s1  = s1;
      m_s2  = s2;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Set two space points seed
  /////////////////////////////////////////////////////////////////////////////////

  inline void InternalSeed::set2(Seed& s)
    {
      s.erase();
      s.add(m_s0->spacepoint);
      s.add(m_s1->spacepoint);
      s.setZVertex(double(m_z));
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Set three space points seed
  /////////////////////////////////////////////////////////////////////////////////

  inline bool InternalSeed::set3(Seed& s)
    {
      
      bool pixb = !m_s0->spacepoint->clusterList().second;
      bool pixt = !m_s2->spacepoint->clusterList().second;
      
      if(pixb!=pixt) {
    if(m_q > m_s0->quality() && m_q > m_s1->quality() && m_q > m_s2->quality()) return false;
      }
     
      m_s0->setQuality(m_q);
      m_s1->setQuality(m_q);
      m_s2->setQuality(m_q);
      
      s.erase();
      s.add(m_s0->spacepoint);
      s.add(m_s1->spacepoint);
      s.add(m_s2->spacepoint);
      s.setZVertex(double(m_z)); 
      return true;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Set quality in pro seed
  /////////////////////////////////////////////////////////////////////////////////

  inline void InternalSeed::setQuality(float q, bool setSpQuality)
    {
      m_q = q;
      if(setSpQuality) {
        m_s0->setQuality(q);
        m_s1->setQuality(q);
        m_s2->setQuality(q);
      }
    }

} // end of Seeding namespace
} // end of Acts namespace

#endif  // InternalSeed_h
