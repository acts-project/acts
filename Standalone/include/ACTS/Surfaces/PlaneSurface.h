///////////////////////////////////////////////////////////////////
// PlaneSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PLANESURFACE_H
#define ACTS_SURFACES_PLANESURFACE_H 1

// Geometry module
#include "ACTS/Surface.h"
#include "ACTS/Surfaces/PlanarBounds.h"
#include "ACTS/Surfaces/NoBounds.h"
// Core Module
#include "ACTS/Utilities/AlgebraDefinitions.h"
#include "ACTS/Utilities/Identifier.h"

namespace Acts {

  class DetectorElementBase;

  /**
   @class PlaneSurface
   Class for a planaer rectangular or trapezoidal surface in the TrackingGeometry.
   It inherits from Surface.

   The Acts::PlaneSurface extends the Surface class with the possibility to convert
   in addition to local to global positions, also local to global direction (vice versa).

   @image html PlaneSurface.gif

   @author Andreas.Salzburger@cern.ch
   */

  class PlaneSurface : public Surface {
    public:

      /** Default Constructor - needed for persistency*/
      PlaneSurface();

      /** Copy Constructor*/
      PlaneSurface(const PlaneSurface& psf);

      /** Copy Constructor with shift*/
      PlaneSurface(const PlaneSurface& psf, const Transform3D& transf);

      /** Dedicated Constructor with normal vector */
      PlaneSurface(const Vector3D& position, const Vector3D& normal);

      /** Constructor from DetectorElementBase - potentially with identifier */
      PlaneSurface(const DetectorElementBase& detelement, const Identifier& identifier = Identifier());

      /** Constructor for planar Surface without Bounds */
      PlaneSurface(std::shared_ptr<Transform3D> htrans);

      /** Constructor for planar Surface from unique_ptr without Bounds */
      PlaneSurface(std::unique_ptr<Transform3D> htrans);

      /** Constructor for Rectangular Planes*/
      PlaneSurface(std::shared_ptr<Transform3D> htrans, double halephi, double haleta);

      /** Constructor for Trapezoidal Planes*/
      PlaneSurface(std::shared_ptr<Transform3D> htrans, double minhalephi, double maxhalephi, double haleta);

      /** Constructor for Planes with a pointer - passing ownership */
      PlaneSurface(std::shared_ptr<Transform3D> htrans, const PlanarBounds* pbounds);

      /** Constructor for Planes with shared bounds object */
      PlaneSurface(std::shared_ptr<Transform3D> htrans, std::shared_ptr<const PlanarBounds> pbounds);

      /** Destructor*/
      virtual ~PlaneSurface();

      /** Assignment operator*/
      PlaneSurface& operator=(const PlaneSurface& psf);

      /** Equality operator*/
      virtual bool operator==(const Surface& sf) const override;

      /** Virtual constructor - shift is optionally */
      virtual PlaneSurface* clone(const Transform3D* shift = nullptr) const override;

      /** Return the surface type */
      virtual SurfaceType type() const override { return Surface::Plane; }

      /**This method returns the bounds by reference, static NoBounds in case of no boundaries*/
      virtual const SurfaceBounds& bounds() const override;

      /**This method calls the inside() method of the Bounds*/
      virtual bool insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const override;

      /** This method returns true if the GlobalPosition is on the Surface for both, within
        or without check of whether the local position is inside boundaries or not */
      virtual bool isOnSurface(const Vector3D& glopo, const BoundaryCheck& bchk=true) const override;

      /** Specified for PlaneSurface: LocalToGlobal method without dynamic memory allocation */
      virtual void localToGlobal(const Vector2D& locp, const Vector3D& mom, Vector3D& glob) const override;

      /** Specified for PlaneSurface: GlobalToLocal method without dynamic memory allocation - boolean checks if on surface */
      virtual bool globalToLocal(const Vector3D& glob, const Vector3D& mom, Vector2D& loc) const override;

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
                                                    bool forceDir,
                                                    const BoundaryCheck& bchk = true) const override;

      /** Return properly formatted class name for screen output */
      virtual std::string name() const  override { return "Acts::PlaneSurface"; }

    protected: //!< data members
      mutable std::shared_ptr<const PlanarBounds>   m_bounds;      //!< bounds (shared)
      static NoBounds                               s_boundless;   //!< NoBounds as return object when no bounds are declared

    };


  inline PlaneSurface* PlaneSurface::clone(const Transform3D* shift) const
  {
    if (shift) new PlaneSurface(*this,*shift);
    return new PlaneSurface(*this);
  }

  inline bool PlaneSurface::insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const
  { return (bounds().inside(locpos,bchk));}

  inline const SurfaceBounds& PlaneSurface::bounds() const
  {
    if (m_bounds.get()) return (*m_bounds.get());
    if (Surface::m_associatedDetElement && Surface::m_associatedDetElementId.is_valid()){
     return m_associatedDetElement->bounds(Surface::m_associatedDetElementId);
    }
    if (Surface::m_associatedDetElement) return m_associatedDetElement->bounds();
    return s_boundless;
  }


  inline Intersection PlaneSurface::intersectionEstimate(const Vector3D& pos,
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

#endif // ACTS_SURFACES_PLANESURFACE_H
