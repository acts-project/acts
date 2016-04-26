///////////////////////////////////////////////////////////////////
// SubtractedCylinderSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SUBTRACTEDCYLINDERSURFACE_H
#define ACTS_SURFACES_SUBTRACTEDCYLINDERSURFACE_H 1

#include "ACTS/Surfaces/CylinderSurface.h"
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Utilities/AreaExcluder.h"
#include "ACTS/Utilities/ParameterDefinitions.h"

class Identifier;

namespace Acts {

  /**
   @class SubtractedCylinderSurface
   Class for a cylinder subtracted/shared surface in the ATLAS detector.
   It owns its surface bounds and subtracted volume.

   @author Sarka.Todorova@cern.ch
   */

  class SubtractedCylinderSurface : public CylinderSurface {
    public:
      /** Default Constructor - needed for persistency*/
      SubtractedCylinderSurface();

      /** Copy Constructor*/
      SubtractedCylinderSurface(const SubtractedCylinderSurface& psf);

      /** Copy Constructor with shift*/
      SubtractedCylinderSurface(const SubtractedCylinderSurface& psf, const Transform3D& transf);

      /** Constructor */
      SubtractedCylinderSurface(const CylinderSurface& cs, AreaExcluder* vol, bool shared);

      /**Destructor*/
      virtual ~SubtractedCylinderSurface();

      /**Assignment operator*/
      SubtractedCylinderSurface& operator=(const SubtractedCylinderSurface& psf);

      /**Equality operator*/
      bool operator==(const Surface& sf) const;

      /** This method indicates the subtraction mode */
      bool shared() const;

      /**This method calls the inside() method of the Bounds*/
      bool insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk=true) const;

      /**This method allows access to the subtracted part*/
      std::shared_ptr<AreaExcluder> subtractedVolume() const;

      /** Return properly formatted class name for screen output */
      std::string name() const { return "Acts::SubtractedCylinderSurface"; }

    protected:
      std::shared_ptr<AreaExcluder>           m_subtrVol;
      bool                                    m_shared;
    };

  inline bool SubtractedCylinderSurface::insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const
  {
    // no subtracted volume exists
    if (!m_subtrVol.get()) return (bounds().inside(locpos,bchk));
    // subtracted volume exists, needs to be checked
    double rCyl    = bounds().r();
    double phiPos  = locpos[Acts::eLOC_RPHI]/rCyl;
    const Vector3D gp(rCyl*cos(phiPos),rCyl*sin(phiPos),locpos[Acts::eLOC_Z]);

    bool inside_shared(this->bounds().inside(locpos,bchk) && m_subtrVol->inside(gp,0.) );
    bool inside(this->bounds().inside(locpos,bchk) && !m_subtrVol->inside(gp,0.) );

    if (m_shared) return inside_shared;
    return inside;
  }

  inline bool SubtractedCylinderSurface::shared() const { return m_shared;}

  inline std::shared_ptr< AreaExcluder> SubtractedCylinderSurface::subtractedVolume() const
  {
    return m_subtrVol;
  }

} // end of namespace

#endif // ACTS_SURFACES_SUBTRACTEDCYLINDERSURFACE_H


