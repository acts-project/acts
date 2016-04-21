///////////////////////////////////////////////////////////////////
// ConeSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_CONESURFACE_H
#define ACTS_SURFACES_CONESURFACE_H 1

// Geometry module
#include "ACTS/Utilities/ParameterDefinitions.h"
#include "ACTS/Surfaces/Surface.h"
#include "ACTS/Surfaces/ConeBounds.h"
// EventData module
#include "ACTS/Utilities/AlgebraDefinitions.h"

namespace Acts {

  /**
   @class ConeSurface

   Class for a conical surface in the Tracking geometry.
   It inherits from Surface.

   The ConeSurface is special since no corresponding
   TrackParameters exist since they're numerical instable
   at the tip of the cone.
   Propagations to a cone surface will be returned in
   curvilinear coordinates.

   @author Ian.Watson@cern.ch, Andreas.Salzburger@cern.ch
   */

  class ConeSurface : public Surface {

    public:
      /** Default Constructor*/
      ConeSurface();

      /** Constructor form HepTransform and an opening angle */
      ConeSurface(std::shared_ptr<Transform3D> htrans, double alpha, bool symmetric=false);

      /** Constructor form HepTransform, radius halfphi, and halflenght*/
      ConeSurface(std::shared_ptr<Transform3D> htrans, double alpha, double locZmin, double locZmax, double halfPhi=M_PI);

      /** Constructor from HepTransform and CylinderBounds
        - ownership of the bounds is passed */
      ConeSurface(std::shared_ptr<Transform3D> htrans, std::shared_ptr<const ConeBounds> cbounds);

      /** Constructor from HepTransform by unique_ptr.
         - bounds is not set. */
      ConeSurface(std::unique_ptr<Transform3D> htrans);

      /** Copy constructor */
      ConeSurface(const ConeSurface& csf);

      /** Copy constructor with shift */
      ConeSurface(const ConeSurface& csf, const Transform3D& transf);

      /** Destructor*/
      virtual ~ConeSurface();

      /** Assignment operator*/
      ConeSurface& operator=(const ConeSurface& csf);

      /** Equality operator*/
      virtual bool operator==(const Surface& sf) const override;

      /** Implicit Constructor*/
      virtual ConeSurface* clone(const Transform3D* shift = nullptr) const override;

      /** The binning position method - is overloaded for r-type binning */
      virtual Vector3D binningPosition(BinningValue bValue) const override;

      /** Return the surface type */
      virtual SurfaceType type() const override { return Surface::Cone; }

      /** Return the measurement frame - this is needed for alignment, in particular for StraightLine and Perigee Surface
        - the default implementation is the the RotationMatrix3D of the transform */
      virtual const RotationMatrix3D measurementFrame(const Vector3D& glopos, const Vector3D& glomom) const override;

      /**Return method for surface normal information
         at a given local point, overwrites the normal() from base class.*/
      virtual const Vector3D& normal() const override;

      /**Return method for surface normal information
         at a given local point, overwrites the normal() from base class.*/
      virtual const Vector3D normal(const Vector2D& locpo) const override;

      /**Return method for surface normal information
         at a given local point, overwrites the normal() from base class.*/
      virtual const Vector3D normal(const Vector3D& global) const override;

      /**Return method for the rotational symmetry axis - the z-Axis of the HepTransform */
      virtual const Vector3D& rotSymmetryAxis() const;

      /**This method returns the ConeBounds by reference
       (NoBounds is not possible for cone)*/
      virtual const ConeBounds& bounds() const override;

      /** Specialized for ConeSurface : LocalToGlobal method without dynamic memory allocation */
      virtual void localToGlobal(const Vector2D& locp, const Vector3D& mom, Vector3D& glob) const override;

      /** Specialized for ConeSurface : GlobalToLocal method without dynamic memory allocation - boolean checks if on surface */
      virtual bool globalToLocal(const Vector3D& glob, const Vector3D& mom, Vector2D& loc) const override;

      /** fast straight line intersection schema - provides closest intersection and (signed) path length

      <b>mathematical motivation:</b>

        The calculation will be done in the 3-dim frame of the cone,
        i.e. the symmetry axis of the cone is the z-axis, x- and y-axis are perpenticular
        to the the z-axis. In this frame the cone is centered around the origin.
        Therefore the two points describing the line have to be first recalculated into the new frame.
        Suppose, this is done, the points of intersection can be
        obtained as follows:<br>

        The cone is described by the implicit equation
        @f$x^2 + y^2 = z^2 \tan \alpha@f$
        where @f$\alpha@f$ is opening half-angle of the cone  the and
        the line by the parameter equation (with @f$t@f$ the
        parameter and @f$x_1@f$ and @f$x_2@f$ are points on the line)
        @f$(x,y,z) = \vec x_1 + (\vec x_2 - \vec x_2) t @f$.
        The intersection is the given to the value of @f$t@f$ where
        the @f$(x,y,z)@f$ coordinates of the line satisfy the implicit
        equation of the cone. Inserting the expression for the points
        on the line into the equation of the cone and rearranging to
        the form of a  gives (letting @f$ \vec x_d = \frac{\vec x_2 - \vec
        x_1}{\abs{\vec x_2 - \vec x_1}} @f$):
        @f$t^2 (x_d^2 + y_d^2 - z_d^2 \tan^2 \alpha) + 2 t (x_1 x_d +
        y_1 y_d - z_1 z_d \tan^2 \alpha) + (x_1^2 + y_1^2 - z_1^2
        \tan^2 \alpha) = 0 @f$
        Solving the above for @f$t@f$ and putting the values into the
        equation of the line gives the points of intersection. @f$t@f$
        is also the length of the path, since we normalized @f$x_d@f$
        to be unit length.
      */
      virtual Intersection intersectionEstimate(const Vector3D& pos,const Vector3D& dir, bool forceDir=false, const BoundaryCheck& bchk = false) const override;

      /** the pathCorrection for derived classes with thickness */
      virtual double pathCorrection(const Vector3D&, const Vector3D&) const override;

      /** Return properly formatted class name for screen output */
      virtual std::string name() const override { return "Acts::ConeSurface"; }

    protected: //!< data members
      std::shared_ptr<const ConeBounds>   m_bounds;                //!< bounds (shared)
      mutable Vector3D*                   m_rotSymmetryAxis;       //!< The rotational symmetry axis
  };

  inline ConeSurface* ConeSurface::clone(const Transform3D* shift) const
  {
    if (shift) new ConeSurface(*this,*shift);
    return new ConeSurface(*this);
  }

  inline const Vector3D& ConeSurface::normal() const
  { return Surface::normal(); }

  inline const Vector3D ConeSurface::normal(const Vector2D& lp) const
  {
    // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
    double phi = lp[Acts::eLOC_RPHI]/(bounds().r(lp[Acts::eLOC_Z])),
           sgn = lp[Acts::eLOC_Z] > 0 ? -1. : +1.;
    Vector3D localNormal(cos(phi) * bounds().cosAlpha(),
				sin(phi) * bounds().cosAlpha(),
				sgn*bounds().sinAlpha());
    return Vector3D(transform().rotation()*localNormal);
  }

  inline const Vector3D ConeSurface::normal(const Vector3D& gp) const {
      Vector2D local(0.,0.);
      globalToLocal(gp,Vector3D::UnitX(),local);
      return normal(local);
  }

  inline const ConeBounds& ConeSurface::bounds() const {
      return (*m_bounds.get());
  }

} // end of namespace

#endif // ACTS_SURFACESCONESURFACE_H


