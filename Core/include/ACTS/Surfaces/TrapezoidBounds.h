///////////////////////////////////////////////////////////////////
// TrapezoidBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_TRAPEZOIDBOUNDS_H
#define ACTS_SURFACES_TRAPEZOIDBOUNDS_H 1

#include "ACTS/Surfaces/PlanarBounds.h"
#include "ACTS/Utilities/PrecisionDefinition.h"
#include "ACTS/Utilities/ParameterDefinitions.h"
#include "ACTS/Utilities/AlgebraDefinitions.h"
#include <math.h>

namespace Acts {

  /**
   @class TrapezoidBounds
   Bounds for a trapezoidal, planar Surface.

   @image html TrapezoidalBounds.gif
      <br>
   @image html SCT_Trapezoid.gif

   @todo can be speed optimized by calculating kappa/delta and caching it

   @author Andreas.Salzburger@cern.ch
   */

  class TrapezoidBounds : public PlanarBounds {

    public:
      /** @enum BoundValues - for readability */
      enum BoundValues {
          bv_minHalfX = 0,
          bv_maxHalfX = 1,
          bv_halfY    = 2,
          bv_length   = 3
      };

      /** Default Constructor, needed for persistency*/
      TrapezoidBounds();

      /** Constructor for symmetric Trapezoid*/
      TrapezoidBounds(double minhalex, double maxhalex,  double haley);

      /** Constructor for arbitrary Trapezoid*/
      TrapezoidBounds(double minhalex, double haley, double alpha, double beta);

      /** Copy constructor*/
      TrapezoidBounds(const TrapezoidBounds& trabo);

      /** Destructor*/
      virtual ~TrapezoidBounds();

      /** Virtual constructor*/
      virtual TrapezoidBounds* clone() const override;

      /** Return the type of the bounds for persistency */
      virtual BoundsType type() const override { return SurfaceBounds::Trapezoid; }

      /** Assignment operator*/
      TrapezoidBounds& operator=(const TrapezoidBounds& sbo);

      /** Equality operator*/
      virtual bool operator==(const SurfaceBounds& trabo) const override;

      /** This method returns the minimal halflength in phi (first coordinate of local surface frame)*/
      double minHalflengthPhi() const;

      /** This method returns the maximal halflength in phi (first coordinate of local surface frame)*/
      double maxHalflengthPhi() const;

      /** This method returns the halflength in eta (second coordinate of local surface frame)*/
      double halflengthEta() const;

      /** This method returns the minimal halflength in X (first coordinate of local surface frame)*/
      double minHalflengthX() const;

      /** This method returns the maximal halflength in X (first coordinate of local surface frame)*/
      double maxHalflengthX() const;

      /** This method returns the halflength in Y (second coordinate of local surface frame)*/
      double halflengthY() const;

      /** This method returns the maximal extension on the local plane*/
      virtual double r() const override;

      /** This method returns the opening angle alpha in point A (negative local phi) */
      double alpha() const;

      /** This method returns the opening angle beta in point B (positive local phi) */
      double beta() const;

      /** The orientation of the Trapezoid is according to the figure above,
       in words: the shorter of the two parallel sides of the trapezoid intersects
       with the negative @f$ y @f$ - axis of the local frame.

       <br>
       The cases are:<br>
       (0) @f$ y @f$ or @f$ x @f$ bounds are 0 || 0<br>
       (1) the local position is outside @f$ y @f$ bounds <br>
       (2) the local position is inside @f$ y @f$ bounds, but outside maximum @f$ x @f$ bounds  <br>
       (3) the local position is inside @f$ y @f$ bounds AND inside minimum @f$ x @f$ bounds <br>
       (4) the local position is inside @f$ y @f$ bounds AND inside maximum @f$ x @f$ bounds, so that
       it depends on the @f$ eta @f$ coordinate
       (5) the local position fails test of (4) <br>

       The inside check is done using single equations of straight lines and one has to take care if a point
       lies on the positive @f$ x @f$ half area(I) or the negative one(II). Denoting @f$ |x_{min}| @f$ and
       @f$ | x_{max} | @f$ as \c minHalfX respectively \c maxHalfX, such as @f$ | y_{H} | @f$ as \c halfY,
       the equations for the straing lines in (I) and (II) can be written as:<br>
        <br>
       - (I):  @f$ y = \kappa_{I} x + \delta_{I} @f$ <br>
       - (II): @f$ y = \kappa_{II} x + \delta_{II} @f$ ,<br>
        <br>
       where @f$  \kappa_{I} = - \kappa_{II} = 2 \frac{y_{H}}{x_{max} - x_{min}} @f$ <br>
       and   @f$  \delta_{I} = \delta_{II} = - \frac{1}{2}\kappa_{I}(x_{max} + x_{min}) @f$  */
      virtual bool inside(const Vector2D& locpo, double tol1=0., double tol2=0.) const override;
      virtual bool inside(const Vector2D& locpo, const BoundaryCheck& bchk) const override;

      /** This method checks inside bounds in loc1
      - loc1/loc2 correspond to the natural coordinates of the surface
      - As loc1/loc2 are correlated the single check doesn't make sense :
         -> check is done on enclosing Rectangle ! */

      virtual bool insideLoc1(const Vector2D& locpo, double tol1=0.) const override;

      /** This method checks inside bounds in loc2
      - loc1/loc2 correspond to the natural coordinates of the surface
      - As loc1/loc2 are correlated the single check doesn't make sense :
         -> check is done on enclosing Rectangle !  */
      virtual bool insideLoc2(const Vector2D& locpo, double tol2=0.) const override;


      /** Minimal distance to boundary ( > 0 if outside and <=0 if inside) */
      virtual double minDistance(const Vector2D& pos) const override;
      
      /** Return the vertices - or, the points of the extremas */
      virtual const std::vector< Vector2D > vertices() const override;

      /** Output Method for std::ostream */
      virtual std::ostream& dump(std::ostream& sl) const override;

   private:
      /** inside() method for a full symmetric trapezoid */
      bool insideFull(const Vector2D& locpo, double tol1=0., double tol2=0.) const;

      /** inside() method for the triangular exclude area for an arbitrary trapezoid */
      bool insideExclude(const Vector2D& locpo, double tol1=0., double tol2=0.) const;

      /** isAbove() method for checking whether a point lies above or under a straight line */
      bool isAbove(const Vector2D& locpo,
                   double tol1,
                   double tol2,
                   double k,
                   double d) const;

      std::vector<TDD_real_t>   m_boundValues;
      TDD_real_t                m_alpha;
      TDD_real_t                m_beta;

  };

  inline TrapezoidBounds* TrapezoidBounds::clone() const { return new TrapezoidBounds(*this); }

  inline double TrapezoidBounds::minHalflengthX() const { return m_boundValues[TrapezoidBounds::bv_minHalfX]; }

  inline double TrapezoidBounds::maxHalflengthX() const { return m_boundValues[TrapezoidBounds::bv_maxHalfX]; }

  inline double TrapezoidBounds::halflengthY() const    { return m_boundValues[TrapezoidBounds::bv_halfY]; }

  inline double TrapezoidBounds::minHalflengthPhi() const { return minHalflengthX(); }

  inline double TrapezoidBounds::maxHalflengthPhi() const { return maxHalflengthX(); }

  inline double TrapezoidBounds::halflengthEta() const { return halflengthY(); }

  inline double TrapezoidBounds::alpha() const { return m_alpha; }

  inline double TrapezoidBounds::beta() const { return m_beta; }

  inline double TrapezoidBounds::r() const
   { return sqrt(m_boundValues[TrapezoidBounds::bv_maxHalfX]*m_boundValues[TrapezoidBounds::bv_maxHalfX]
               + m_boundValues[TrapezoidBounds::bv_halfY]*m_boundValues[TrapezoidBounds::bv_halfY]); }

  inline bool TrapezoidBounds::inside(const Vector2D& locpo, const BoundaryCheck& bchk) const
  {
	if(bchk.bcType==0)	return TrapezoidBounds::inside(locpo, bchk.toleranceLoc1, bchk.toleranceLoc2);

	// a fast FALSE
	double fabsY = fabs(locpo[Acts::eLOC_Y]);
	double max_ell = bchk.lCovariance(0,0) > bchk.lCovariance(1,1) ? bchk.lCovariance(0,0) :bchk.lCovariance(1,1);
	double limit = bchk.nSigmas*sqrt(max_ell);
	if (fabsY > ( m_boundValues[TrapezoidBounds::bv_halfY] + limit)) return false;
	// a fast FALSE
	double fabsX = fabs(locpo[Acts::eLOC_X]);
	if (fabsX > (m_boundValues[TrapezoidBounds::bv_maxHalfX] + limit)) return false;
	// a fast TRUE
	double min_ell = bchk.lCovariance(0,0) < bchk.lCovariance(1,1) ? bchk.lCovariance(0,0) : bchk.lCovariance(1,1);
	limit = bchk.nSigmas*sqrt(min_ell);
	if (fabsX < (m_boundValues[TrapezoidBounds::bv_minHalfX] + limit) && fabsY < (m_boundValues[TrapezoidBounds::bv_halfY] + limit)) return true;

	// compute KDOP and axes for surface polygon
    std::vector<KDOP> elementKDOP(3);
    std::vector<Vector2D> elementP(4);
    float theta = (bchk.lCovariance(1,0) != 0 && (bchk.lCovariance(1,1)-bchk.lCovariance(0,0))!=0 ) ? .5*bchk.FastArcTan( 2*bchk.lCovariance(1,0)/(bchk.lCovariance(1,1)-bchk.lCovariance(0,0)) ) : 0.;
    sincosCache scResult = bchk.FastSinCos(theta);
    ActsMatrixD<2,2> rotMatrix ;
    rotMatrix << scResult.cosC, scResult.sinC,
                -scResult.sinC, scResult.cosC;
	ActsMatrixD<2,2> normal ;
    normal    << 0, -1,
                 1,  0;
	// ellipse is always at (0,0), surface is moved to ellipse position and then rotated
    Vector2D p;
    p << m_boundValues[TrapezoidBounds::bv_minHalfX],-m_boundValues[TrapezoidBounds::bv_halfY];
    elementP[0] =( rotMatrix * (p - locpo) );
    p << -m_boundValues[TrapezoidBounds::bv_minHalfX],-m_boundValues[TrapezoidBounds::bv_halfY];
    elementP[1] =( rotMatrix * (p - locpo) );
	scResult = bchk.FastSinCos(m_beta);
    p << m_boundValues[TrapezoidBounds::bv_minHalfX] + (2.*m_boundValues[TrapezoidBounds::bv_halfY])*(scResult.sinC/scResult.cosC),m_boundValues[TrapezoidBounds::bv_halfY];
    elementP[2] =( rotMatrix * (p - locpo) );
	scResult = bchk.FastSinCos(m_alpha);
    p << -(m_boundValues[TrapezoidBounds::bv_minHalfX] + (2.*m_boundValues[TrapezoidBounds::bv_halfY])*(scResult.sinC/scResult.cosC)),m_boundValues[TrapezoidBounds::bv_halfY];
    elementP[3] =( rotMatrix * (p - locpo) );
    std::vector<Vector2D> axis = {normal*(elementP[1]-elementP[0]), normal*(elementP[3]-elementP[1]), normal*(elementP[2]-elementP[0])};
    bchk.ComputeKDOP(elementP, axis, elementKDOP);
	// compute KDOP for error ellipse
    std::vector<KDOP> errelipseKDOP(3);
	bchk.ComputeKDOP(bchk.EllipseToPoly(3), axis, errelipseKDOP);
	// check if KDOPs overlap and return result
	return bchk.TestKDOPKDOP(elementKDOP, errelipseKDOP);
  }

  inline bool TrapezoidBounds::insideLoc1(const Vector2D &locpo, double tol1) const
    { return (fabs(locpo[Acts::eLOC_X]) < m_boundValues[TrapezoidBounds::bv_maxHalfX] + tol1); }

  inline bool TrapezoidBounds::insideLoc2(const Vector2D &locpo, double tol2) const
    { return (fabs(locpo[Acts::eLOC_Y]) < m_boundValues[TrapezoidBounds::bv_halfY] + tol2); }

  inline const std::vector< Vector2D > TrapezoidBounds::vertices() const {
      // create the return vector
      std::vector< Vector2D > vertices;
      // fill the vertices
      vertices.reserve(4);
      vertices.push_back(Vector2D(m_boundValues[TrapezoidBounds::bv_minHalfX],-m_boundValues[TrapezoidBounds::bv_halfY]));   // [0]
      vertices.push_back(Vector2D(m_boundValues[TrapezoidBounds::bv_maxHalfX], m_boundValues[TrapezoidBounds::bv_halfY]));   // [1]
      vertices.push_back(Vector2D(-m_boundValues[TrapezoidBounds::bv_maxHalfX],m_boundValues[TrapezoidBounds::bv_halfY]));   // [1]
      vertices.push_back(Vector2D(-m_boundValues[TrapezoidBounds::bv_minHalfX],-m_boundValues[TrapezoidBounds::bv_halfY]));  // [3]
      return vertices;
      
  }

} // end of namespace

#endif // ACTS_SURFACES_TRAPEZOIDBOUNDS_H

