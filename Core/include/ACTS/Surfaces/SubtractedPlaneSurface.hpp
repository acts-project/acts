///////////////////////////////////////////////////////////////////
// SubtractedPlaneSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SUBTRACTEDPLANESURFACE_H
#define ACTS_SURFACES_SUBTRACTEDPLANESURFACE_H

#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Utilities/AreaExcluder.hpp"
#include "ACTS/Utilities/Definitions.hpp"

class Identifier;

namespace Acts {

  /**
   @class SubtractedPlaneSurface

   Class for a planar subtracted/shared surface in the ATLAS detector.
   It owns its surface bounds and subtracted volume.

   @author Sarka.Todorova@cern.ch
   */

  class SubtractedPlaneSurface : public PlaneSurface {
    public:
      /** Default Constructor - needed for persistency*/
      SubtractedPlaneSurface();

      /** Copy Constructor*/
      SubtractedPlaneSurface(const SubtractedPlaneSurface& psf);

      /** Copy Constructor with shift*/
      SubtractedPlaneSurface(const SubtractedPlaneSurface& psf, const Transform3D& transf);

      /** Constructor */
      SubtractedPlaneSurface(const PlaneSurface& ps , AreaExcluder* vol, bool shared);

      /**Destructor*/
      virtual ~SubtractedPlaneSurface();

      /**Assignment operator*/
      SubtractedPlaneSurface& operator=(const SubtractedPlaneSurface& psf);

      /**Equality operator*/
      bool operator==(const Surface& sf) const;

      /** This method indicates the subtraction mode */
      bool shared() const;

      /**This method calls the inside() method of the Bounds*/
      bool insideBounds(const Vector2D& locpos,  const BoundaryCheck& bchk=true) const;

      /**This method allows access to the subtracted part*/
      std::shared_ptr<AreaExcluder>  subtractedVolume() const;

      /** Return properly formatted class name for screen output */
      std::string name() const { return "Acts::SubtractedPlaneSurface"; }

    protected:
      std::shared_ptr<AreaExcluder>          m_subtrVol;
      bool                                   m_shared;
    };

  inline bool SubtractedPlaneSurface::insideBounds(const Vector2D& locpos,  const BoundaryCheck& bchk) const
  {
    // no subtracted volume exists
    if (!m_subtrVol.get()) return (this->bounds().inside(locpos,bchk));
    // subtracted volume exists, needs to be checked
    Vector3D gp(locpos.x(),locpos.y(),0.);
    if (m_shared) return (this->bounds().inside(locpos,bchk) && m_subtrVol->inside(gp,0.) );
    bool in(this->bounds().inside(locpos,bchk) && !m_subtrVol->inside(gp,0.)) ;
    return in;
  }

  inline bool SubtractedPlaneSurface::shared() const { return m_shared;}

  inline std::shared_ptr<AreaExcluder> SubtractedPlaneSurface::subtractedVolume() const
  {
    return m_subtrVol;
  }

} // end of namespace

#endif // TRKGEOMETRYSURFACES_SUBTRACTEDPLANESURFACE_H


