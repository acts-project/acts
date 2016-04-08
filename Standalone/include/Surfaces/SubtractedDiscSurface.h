///////////////////////////////////////////////////////////////////
// SubtractedDiscSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SUBTRACTEDDISCSURFACE_H
#define ACTS_SURFACES_SUBTRACTEDDISCSURFACE_H 1

// Geometry module
#include "Surfaces/DiscSurface.h"
#include "GeometryUtils/AreaExcluder.h"
// EventData module
#include "CoreUtils/ParameterDefinitions.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"

class MsgStream;

namespace Acts {
  
  /**
   @class SubtractedDiscSurface
   Class for a planar subtracted/shared surface in the ATLAS detector.
   It owns its surface bounds and subtracted volume. 
    
   @author Sarka.Todorova@cern.ch
   */

  class SubtractedDiscSurface : public DiscSurface {
    public:
      /** Default Constructor - needed for persistency*/
      SubtractedDiscSurface();  
      
      /** Copy Constructor*/
      SubtractedDiscSurface(const SubtractedDiscSurface& psf);
      
      /** Copy Constructor*/
      SubtractedDiscSurface(const SubtractedDiscSurface& psf, const Transform3D& shift);
      
      /** Constructor */      
      SubtractedDiscSurface(const DiscSurface& ps , AreaExcluder* vol, bool shared);
      
      /**Destructor*/
      virtual ~SubtractedDiscSurface();  
      
      /**Assignment operator*/
      SubtractedDiscSurface& operator=(const SubtractedDiscSurface& psf);
      
      /**Equality operator*/
      bool operator==(const Surface& sf) const;

      /** This method indicates the subtraction mode */
      bool shared() const;

      /**This method calls the inside() method of the Bounds*/
      bool insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk=true) const;       

      /**This method allows access to the subtracted part*/
      std::shared_ptr<AreaExcluder> subtractedVolume() const; 
            
      /** Return properly formatted class name for screen output */
      std::string name() const { return "Acts::SubtractedDiscSurface"; }
      
    protected:
      std::shared_ptr<AreaExcluder>     m_subtrVol; 
      bool                              m_shared;  
    };

  inline bool SubtractedDiscSurface::insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const
  {
    // no subtracted Volume exists
    if (!m_subtrVol.get()) return (this->bounds().inside(locpos,bchk)); 
    // subtracted Volume exists, needs to be checked
    double rPos   = locpos[Acts::eLOC_R];
    double phiPos = locpos[Acts::eLOC_PHI];
    Vector3D gp(rPos*cos(phiPos),rPos*sin(phiPos),0.);
    if (m_shared) return (this->bounds().inside(locpos,bchk) && m_subtrVol->inside(gp,0.) );
    bool in(this->bounds().inside(locpos,bchk) && !m_subtrVol->inside(gp,0.)) ;
    return in;
  }

  inline bool SubtractedDiscSurface::shared() const { return m_shared;} 

  inline std::shared_ptr<AreaExcluder> SubtractedDiscSurface::subtractedVolume() const 
  {
    return m_subtrVol;
  }
    
} // end of namespace

#endif // ACTS_SURFACES_SUBTRACTEDDISCSURFACE_H
