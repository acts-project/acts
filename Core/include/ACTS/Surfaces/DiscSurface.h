///////////////////////////////////////////////////////////////////
// DiscSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACE_SDISCSURFACE_H
#define ACTS_SURFACE_SDISCSURFACE_H 1

#include "ACTS/Surfaces/Surface.h"
#include "ACTS/Surfaces/DiscBounds.h"
#include "ACTS/Surfaces/NoBounds.h"
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Utilities/Identifier.h"

namespace Acts {

  class RadialBounds;
  class DiscTrapezoidalBounds;
  class DetectorElementBase;

  /**
   @class DiscSurface

   Class for a DiscSurface in the TrackingGEometry.
   It inherits from Surface.

   @author Andreas.Salzburger@cern.ch
   */

  class DiscSurface : public Surface {

    public:
      /** Default Constructor*/
      DiscSurface();

      /** Constructor for Discs from Transform3D, \f$ r_{min}, r_{max} \f$ */
      DiscSurface(std::shared_ptr<Transform3D> htrans, double rmin, double rmax);

      /** Constructor for Discs from Transform3D, \f$ r_{min}, r_{max}, \phi_{hsec} \f$ */
      DiscSurface(std::shared_ptr<Transform3D> htrans, double rmin, double rmax, double hphisec);

      /** Constructor for Discs from Transform3D, \f$ r_{min}, r_{max}, hx_{min}, hx_{max} \f$
      	 In this case you have DiscTrapezoidalBounds*/
      DiscSurface(std::shared_ptr<Transform3D> htrans, double minhalfx, double maxhalfx, double maxR, double minR, double avephi, double stereo = 0.);

      /** Constructor for Discs from Transform3D and DiscBounds - ownership of bounds is passed */
      DiscSurface(std::shared_ptr<Transform3D> htrans, const RadialBounds* rbounds);

      /** Constructor for Discs from Transform3D and DiscTrapezoidalBounds - ownership of bounds is passed */
      DiscSurface(std::shared_ptr<Transform3D> htrans, const DiscTrapezoidalBounds* dtbounds);

      /** Constructor for Discs from Transform3D and shared DiscBounds */
      DiscSurface(std::shared_ptr<Transform3D> htrans, std::shared_ptr<const DiscBounds> dbounds);

      /** Constructor for Discs from Transform3D by unique_ptr
       - bounds is not set */
      DiscSurface(std::unique_ptr<Transform3D> htrans);

      /** Constructor from DetectorElementBase*/
      DiscSurface(const DetectorElementBase& detelement, const Identifier& identifier = Identifier());

      /** Copy Constructor*/
      DiscSurface(const DiscSurface& psf);

      /** Copy Constructor with shift*/
      DiscSurface(const DiscSurface& psf, const Transform3D& transf);

      /** Destructor*/
      virtual ~DiscSurface();

      /** Assignement operator*/
      DiscSurface& operator=(const DiscSurface&dsf);

      /** Equality operator*/
      virtual bool operator==(const Surface& sf) const override;

      /** Virtual constructor - shift can be given optionally */
      virtual DiscSurface* clone(const Transform3D* shift = nullptr) const override;

      /** Return the surface type */
      virtual SurfaceType type() const override { return Surface::Disc; }

      /**This method returns the bounds by reference*/
      const SurfaceBounds& bounds() const  override;

      /** This method returns true if the GlobalPosition is on the Surface for both, within
        or without check of whether the local position is inside boundaries or not */
      virtual bool isOnSurface(const Vector3D& glopo, const BoundaryCheck& bchk=true) const  override;

      /** Specialized for DiscSurface: LocalToGlobal method without dynamic memory allocation */
      virtual void localToGlobal(const Vector2D& locp, const Vector3D& mom, Vector3D& glob) const  override;

      /** Specialized for DiscSurface: GlobalToLocal method without dynamic memory allocation - boolean checks if on surface */
      virtual bool globalToLocal(const Vector3D& glob, const Vector3D& mom, Vector2D& loc) const  override;

      /**  Special method for DiscSurface : local<->local transformations polar <-> cartesian */
      const Vector2D localPolarToCartesian(const Vector2D& locpol) const;

      /**  Special method for Disc surface : local<->local transformations polar <-> cartesian */
      const Vector2D localCartesianToPolar(const Vector2D& loccart) const;

      /**  Special method for DiscSurface : local<->local transformations polar <-> cartesian */
      const Vector2D localPolarToLocalCartesian(const Vector2D& locpol) const;

      /** Special method for DiscSurface :  local<->global transformation when provided cartesian coordinates */
      const Vector3D localCartesianToGlobal(const Vector2D& locpos) const;

      /** Special method for DiscSurface : global<->local from cartesian coordinates */
      const Vector2D globalToLocalCartesian(const Vector3D& glopos, double tol=0.) const;

      /** fast straight line intersection schema - standard: provides closest intersection and (signed) path length
          forceDir is to provide the closest forward solution

          <b>mathematical motivation:</b>

          the equation of the plane is given by: <br>
          @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
          where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of the plane,
          @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on the plane and @f$ \vec x = (x,y,z) @f$ all possible points
          on the plane.<br>
          Given a line with:<br>
          @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
          the solution for @f$ u @f$ can be written:
          @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
          If the denominator is 0 then the line lies:
          - either in the plane
          - perpenticular to the normal of the plane
       */
      virtual Intersection intersectionEstimate(const Vector3D& pos,
                                                    const Vector3D& dir,
                                                    bool forceDir = false,
                                                    const BoundaryCheck& bchk = false) const  override;

      /** Return properly formatted class name for screen output */
      virtual std::string name() const override { return "Acts::DiscSurface"; }

    protected: //!< data members
      mutable std::shared_ptr<const DiscBounds>     m_bounds;          //!< bounds (shared)
      mutable Vector3D*                             m_referencePoint;  //!< a reference point on the Surface
      static NoBounds                               s_boundless;       //!< static member for boundless approach
 };

  inline DiscSurface* DiscSurface::clone(const Transform3D* shift) const
  {
    if (shift) return new DiscSurface(*this,*shift);
    return new DiscSurface(*this);
  }

  inline const SurfaceBounds& DiscSurface::bounds() const
  {
    if (m_bounds) return (*(m_bounds.get()));
    if (Surface::m_associatedDetElement && Surface::m_associatedDetElementId.is_valid()){
     return m_associatedDetElement->bounds(Surface::m_associatedDetElementId);
    }
    if (Surface::m_associatedDetElement) return m_associatedDetElement->bounds();
    return s_boundless;
  }

  inline const Vector2D DiscSurface::localPolarToCartesian(const Vector2D& locpol) const
  { return(Vector2D(locpol[Acts::eLOC_R]*cos(locpol[Acts::eLOC_PHI]),locpol[Acts::eLOC_R]*sin(locpol[Acts::eLOC_PHI]))); }

  inline const Vector2D DiscSurface::localCartesianToPolar(const Vector2D& loccart) const
  {
    return (Vector2D(sqrt(loccart[Acts::eLOC_X]*loccart[Acts::eLOC_X]+loccart[Acts::eLOC_Y]*loccart[Acts::eLOC_Y]),
                          atan2(loccart[Acts::eLOC_Y], loccart[Acts::eLOC_X])));
  }

  inline Intersection DiscSurface::intersectionEstimate(const Vector3D& pos,
                                                            const Vector3D& dir,
                                                            bool forceDir,
                                                            const BoundaryCheck& bchk) const
  {
      double denom = dir.dot(normal());
      if (denom){
        double u = (normal().dot((center()-pos)))/(denom);
        Vector3D intersectPoint(pos + u * dir);
        // evaluate the intersection in terms of direction
        bool isValid = forceDir ?  ( u > 0.) : true;
        // evaluate (if necessary in terms of boundaries)
        isValid = bchk ? (isValid && isOnSurface(intersectPoint,bchk) ) : isValid;
        // return the result
        return Intersection(intersectPoint,u,isValid);
      }
      return Intersection(pos,0.,false);
  }

} // end of namespace

#endif // ACTS_SURFACES_DISCSURFACE_H
