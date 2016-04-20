/////////////////////////////////////////////////////////////////
// PerigeeSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PERIGEESURFACE_H
#define ACTS_SURFACES_PERIGEESURFACE_H 1

// Geometry  module
#include "ACTS/Surfaces/Surface.h"
#include "ACTS/Surfaces/NoBounds.h"
// EventData module
#include "ACTS/GeometryUtils/GeometryStatics.h"
// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"

namespace Acts {

  /**
   @class PerigeeSurface

   Class describing the Line to which the Perigee refers to.
   The Surface axis is fixed to be the z-axis of the Tracking frame.
   It inherits from Surface.

   @author Andreas.Salzburger@cern.ch
   */

  class PerigeeSurface : public Surface {

    public:
      /** Default Constructor - needed for persistency*/
      PerigeeSurface();

      /** Constructor from GlobalPosition*/
      PerigeeSurface(const Vector3D& gp);

      /** Constructor with a Transform - needed for tilt */
      PerigeeSurface(std::shared_ptr<Transform3D> tTransform);

      /** Constructor with a Transform by unique_ptr - needed for tilt */
      PerigeeSurface(std::unique_ptr<Transform3D> tTransform);

      /** Copy constructor*/
      PerigeeSurface(const PerigeeSurface& pesf);

      /** Copy constructor with shift*/
      PerigeeSurface(const PerigeeSurface& pesf, const Transform3D& transf);

      /** Destructor*/
      virtual ~PerigeeSurface();

      /** Virtual constructor*/
      virtual PerigeeSurface* clone(const Transform3D* shift = nullptr) const override;

      /** Assignment operator*/
      PerigeeSurface& operator=(const PerigeeSurface& slsf );

      /** Equality operator*/
      virtual bool operator==(const Surface& sf) const override;

      /** Return the surface type */
      virtual SurfaceType type() const override { return Surface::Perigee; }

      /**Return method for transfromation, overwrites the transform() form base class*/
      virtual const Transform3D& transform() const override;

      /**Return method for surface center infromation, overwrites the center() form base class*/
      virtual const Vector3D& center() const override;

      /** Return the measurement frame - this is needed for alignment, in particular for StraightLine and Perigee Surface
          - the default implementation is the the RotationMatrix3D of the transform */
      virtual const RotationMatrix3D measurementFrame(const Vector3D& glopos, const Vector3D& glomom) const override;

      /** LocalToGlobal method without dynamic memory allocation */
      virtual void localToGlobal(const Vector2D& locp, const Vector3D& mom, Vector3D& glob) const override;

      /** GlobalToLocal method without dynamic memory allocation - boolean checks if on surface
         \image html SignOfDriftCircleD0.gif */
      virtual bool globalToLocal(const Vector3D& glob, const Vector3D& mom, Vector2D& loc) const override;

      /** fast straight line intersection schema - standard: provides closest intersection and (signed) path length
          forceDir is to provide the closest forward solution

          b>mathematical motivation:</b>
          Given two lines in parameteric form:<br>
          - @f$ \vec l_{a}(\lambda) = \vec m_a + \lambda \cdot \vec e_{a} @f$ <br>
          - @f$ \vec l_{b}(\mu) = \vec m_b + \mu \cdot \vec e_{b} @f$ <br>
          the vector between any two points on the two lines is given by:
          - @f$ \vec s(\lambda, \mu) = \vec l_{b} - l_{a} = \vec m_{ab} + \mu \cdot \vec e_{b} - \lambda \cdot \vec e_{a} @f$, <br>
          when @f$ \vec m_{ab} = \vec m_{b} - \vec m_{a} @f$.<br>
          @f$ \vec s(\lambda_0, \mu_0) @f$  denotes the vector between the two closest points <br>
          @f$ \vec l_{a,0} = l_{a}(\lambda_0) @f$ and @f$ \vec l_{b,0} = l_{b}(\mu_0) @f$ <br>
          and is perpenticular to both, @f$ \vec e_{a} @f$ and @f$ \vec e_{b} @f$.

          This results in a system of two linear equations:<br>
          - (i) @f$ 0 = \vec s(\lambda_0, \mu_0) \cdot \vec e_a = \vec m_ab \cdot \vec e_a + \mu_0 \vec e_a \cdot \vec e_b - \lambda_0 @f$ <br>
          - (ii) @f$ 0 = \vec s(\lambda_0, \mu_0) \cdot \vec e_b = \vec m_ab \cdot \vec e_b + \mu_0  - \lambda_0 \vec e_b \cdot \vec e_a @f$ <br>

          Solving (i), (ii) for @f$ \lambda_0 @f$ and @f$ \mu_0 @f$ yields:
          - @f$ \lambda_0 = \frac{(\vec m_ab \cdot \vec e_a)-(\vec m_ab \cdot \vec e_b)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
          - @f$ \mu_0 = - \frac{(\vec m_ab \cdot \vec e_b)-(\vec m_ab \cdot \vec e_a)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
       */
      virtual Intersection intersectionEstimate(const Vector3D& pos,
                                                   const Vector3D& dir,
                                                   bool forceDir = false,
                                                   const BoundaryCheck& bchk = false) const override;


      /** the pathCorrection for derived classes with thickness */
      virtual double pathCorrection(const Vector3D&, const Vector3D&) const override { return 1.; }

        /**This method checks if a globalPosition in on the Surface or not*/
      virtual bool isOnSurface(const Vector3D& glopo, const BoundaryCheck& bchk=true) const override;

      /**This surface calls the iside method of the bounds*/
      virtual bool insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const override;

      /** Special method for StraightLineSurface - provides the Line direction from cache: speedup */
      const Vector3D& lineDirection() const;

      /** Return bounds() method */
      virtual const NoBounds& bounds() const override;

      /** Return properly formatted class name for screen output */
      virtual std::string name() const override { return "Acts::PerigeeSurface"; }

      /** Output Method for std::ostream*/
      virtual std::ostream& dump(std::ostream& sl) const override;

    protected: //!< data members

      mutable Vector3D*                      m_lineDirection;  //!< cache of the line direction (speeds up)
      static NoBounds                             s_perigeeBounds;

  };


  inline PerigeeSurface* PerigeeSurface::clone(const Transform3D* shift) const
  {
    if (shift) return new PerigeeSurface(*this,*shift);
    return new PerigeeSurface(*this);
  }

  inline const Transform3D& PerigeeSurface::transform() const
  {
    if (!Surface::m_transform) return(s_idTransform);
    return(*Surface::m_transform);
  }

  inline const Vector3D& PerigeeSurface::center() const
  {
    if (!Surface::m_center && !Surface::m_transform) return(s_origin);
    else if (!Surface::m_center) m_center = new Vector3D(m_transform->translation());
    return(*Surface::m_center);
  }

  inline bool PerigeeSurface::insideBounds(const Vector2D&, const BoundaryCheck& ) const { return true;}

  inline bool PerigeeSurface::isOnSurface(const Vector3D&, const BoundaryCheck&) const {return true;}

  inline const NoBounds& PerigeeSurface::bounds() const { return s_perigeeBounds; }

  inline Intersection PerigeeSurface::intersectionEstimate(const Vector3D& pos,
                                                               const Vector3D& dir,
                                                               bool forceDir,
                                                               const BoundaryCheck&) const
  {
       // following nominclature found in header file and doxygen documentation
       // line one is the straight track
       const Vector3D   ma  = pos;
       const Vector3D&  ea  = dir;
       // line two is the line surface
       const Vector3D& mb = center();
       const Vector3D& eb = lineDirection();
       // now go ahead
       Vector3D  mab(mb - ma);
       double eaTeb = ea.dot(eb);
       double denom = 1 - eaTeb*eaTeb;
       if (fabs(denom)>10e-7){
          double lambda0 = (mab.dot(ea) - mab.dot(eb)*eaTeb)/denom;
          // evaluate the direction, bounds are always true for Perigee
          bool isValid = forceDir ? (lambda0 > 0.) : true;
          return Intersection( (ma+lambda0*ea),lambda0, isValid );
       }
      return Intersection(pos,0.,false);
  }

  inline const Vector3D& PerigeeSurface::lineDirection() const {
       if (m_lineDirection)
           return (*m_lineDirection);
       if (!m_lineDirection && Surface::m_transform) {
           m_lineDirection = new Vector3D(transform().rotation().col(2));
           return (*m_lineDirection);
       }
       return Acts::s_zAxis;

   }

} // end of namespace

#endif // ACTS_SURFACESPERIGEESURFACE_H
