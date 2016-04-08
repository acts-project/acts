///////////////////////////////////////////////////////////////////
// RectangleBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_RECTANGLEBOUNDS_H
#define ACTS_SURFACES_RECTANGLEBOUNDS_H 1

// Geometry module
#include "Surfaces/PlanarBounds.h"
#include "GeometryUtils/PrecisionDefinition.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"

class MsgStream;

namespace Acts {

   /**
    @class RectangleBounds 

    Bounds for a rectangular, planar surface.
    The two local coordinates Acts::eLOC_X, Acts::eLOC_Y are for legacy reasons
    also called @f$ phi @f$ respectively @f$ eta @f$. The orientation 
    with respect to the local surface framce can be seen in the attached
    illustration.

    @image html RectangularBounds.gif
    
    @author Andreas.Salzburger@cern.ch */

  class RectangleBounds : public PlanarBounds {

    public:
        
      /** @enum BoundValues for readability */
      enum BoundValues {
          bv_halfX  = 0,
          bv_halfY  = 1,
          bv_length = 2
      };
        
      /** Default Constructor - needed for persistency*/
      RectangleBounds();  

      /** Constructor with halflength in phi and halflength in eta*/
      RectangleBounds(double halephi, double haleta);
      
      /** Copy constructor*/
      RectangleBounds(const RectangleBounds& recbo);
      
      /** Destructor*/
      virtual ~RectangleBounds();
      
      /** Assignment Operator*/
      RectangleBounds& operator=(const RectangleBounds& recbo);
      
      /** Equality operator*/
      virtual bool operator==(const SurfaceBounds& sbo) const override;
      
      /** Virtual constructor*/
      virtual RectangleBounds* clone() const override;
    
      /** Return the type of the bounds for persistency */
      virtual BoundsType type() const override { return SurfaceBounds::Rectangle; }
    
      /** This method checks if the provided local coordinates are inside the surface bounds*/
      virtual bool inside(const Vector2D &locpo, double tol1=0., double tol2=0.) const override;
      
      /** This method checks if the provided local coordinates are inside the surface bounds*/
      virtual bool inside(const Vector2D& locpo, const BoundaryCheck& bchk) const override;

      /** This method checks inside bounds in loc1
        - loc1/loc2 correspond to the natural coordinates of the surface */
      virtual bool insideLoc1(const Vector2D& locpo, double tol1=0.) const override;

      /** This method checks inside bounds in loc2 
        - loc1/loc2 correspond to the natural coordinates of the surface */
      virtual bool insideLoc2(const Vector2D& locpo, double tol2=0.) const override;

      /** Minimal distance to boundary ( > 0 if outside and <=0 if inside) */
      virtual double minDistance(const Vector2D& pos) const override;

      /** This method returns the halflength in phi (first coordinate of local surface frame)*/
      double halflengthPhi() const;
      
      /** This method returns the halflength in Eta (second coordinate of local surface frame)*/
      double halflengthEta() const;
      
      /* *for consistant naming*/
      double halflengthX() const;
      
      /** for consitant naming*/
      double halflengthY() const;    
      
      /** This method returns the maximal extension on the local plane, i.e. @f$s\sqrt{h_{\phi}^2 + h_{\eta}^2}\f$*/
      virtual double r() const override;
      
      /** Return the vertices - or, the points of the extremas */
      virtual const std::vector< Vector2D > vertices() const override;
    
      /** Output Method for MsgStream*/
      virtual MsgStream& dump(MsgStream& sl) const override;
      
      /** Output Method for std::ostream */
      virtual std::ostream& dump(std::ostream& sl) const override;
  
    private:
      /** The internal version of the bounds can be float/double*/
      std::vector<TDD_real_t>   m_boundValues;

  };

  inline RectangleBounds* RectangleBounds::clone() const
    { return new RectangleBounds(*this); }

  inline bool RectangleBounds::inside(const Vector2D &locpo, double tol1, double tol2) const
    { return ((fabs(locpo[Acts::eLOC_X]) < m_boundValues[RectangleBounds::bv_halfX] + tol1) && (fabs(locpo[Acts::eLOC_Y]) < m_boundValues[RectangleBounds::bv_halfY] + tol2)  ); }
  
  inline bool RectangleBounds::inside(const Vector2D& locpo, const BoundaryCheck& bchk) const
  { 
	if(bchk.bcType==0)	return RectangleBounds::inside(locpo, bchk.toleranceLoc1, bchk.toleranceLoc2);
	
	// a fast FALSE
	double max_ell = bchk.lCovariance(0,0) > bchk.lCovariance(1,1) ? bchk.lCovariance(0,0) : bchk.lCovariance(1,1);
	double limit = bchk.nSigmas*sqrt(max_ell);
    if (!RectangleBounds::inside(locpo, limit, limit)) return false;
	// a fast TRUE
	double min_ell = bchk.lCovariance(0,0) < bchk.lCovariance(1,1) ? bchk.lCovariance(0,0) : bchk.lCovariance(1,1);
	limit = bchk.nSigmas*sqrt(min_ell);
    if (RectangleBounds::inside(locpo, limit, limit)) return true;

	// compute KDOP and axes for surface polygon
    std::vector<KDOP> elementKDOP(4);
    std::vector<Vector2D> elementP(4);
    float theta = (bchk.lCovariance(1,0) != 0 && (bchk.lCovariance(1,1)-bchk.lCovariance(0,0))!=0 ) ? .5*bchk.FastArcTan( 2*bchk.lCovariance(1,0)/(bchk.lCovariance(1,1)-bchk.lCovariance(0,0)) ) : 0.;
    sincosCache scResult = bchk.FastSinCos(theta);
    ActsMatrixD<2,2> rotMatrix ;
    rotMatrix << scResult.cosC, scResult.sinC,
                -scResult.sinC, scResult.cosC;   
	// ellipse is always at (0,0), surface is moved to ellipse position and then rotated
    Vector2D p;
    p << m_boundValues[RectangleBounds::bv_halfX],m_boundValues[RectangleBounds::bv_halfY];
    elementP[0] = ( rotMatrix * (p - locpo) );
    p << m_boundValues[RectangleBounds::bv_halfX],-m_boundValues[RectangleBounds::bv_halfY];
    elementP[1] =( rotMatrix * (p - locpo) );
    p << -m_boundValues[RectangleBounds::bv_halfX],m_boundValues[RectangleBounds::bv_halfY];
    elementP[2] =( rotMatrix * (p - locpo) );
    p << -m_boundValues[RectangleBounds::bv_halfX],-m_boundValues[RectangleBounds::bv_halfY];
    elementP[3] =( rotMatrix * (p - locpo) );
    std::vector<Vector2D> axis = {elementP[0]-elementP[1], elementP[0]-elementP[2], elementP[0]-elementP[3], elementP[1]-elementP[2]};
    bchk.ComputeKDOP(elementP, axis, elementKDOP);
	// compute KDOP for error ellipse
    std::vector<KDOP> errelipseKDOP(4);
	bchk.ComputeKDOP(bchk.EllipseToPoly(3), axis, errelipseKDOP);
	// check if KDOPs overlap and return result
	return bchk.TestKDOPKDOP(elementKDOP, errelipseKDOP);
  }

  inline bool RectangleBounds::insideLoc1(const Vector2D &locpo, double tol1) const
    { return (fabs(locpo[Acts::eLOC_X]) < m_boundValues[RectangleBounds::bv_halfX] + tol1); }

  inline bool RectangleBounds::insideLoc2(const Vector2D &locpo, double tol2) const
    { return (fabs(locpo[Acts::eLOC_Y]) < m_boundValues[RectangleBounds::bv_halfY] + tol2); }

  inline double RectangleBounds::halflengthPhi() const { return this->halflengthX(); }
  
  inline double RectangleBounds::halflengthEta() const { return this->halflengthY(); }

  inline double RectangleBounds::halflengthX() const { return m_boundValues[RectangleBounds::bv_halfX]; }
  
  inline double RectangleBounds::halflengthY() const { return m_boundValues[RectangleBounds::bv_halfY]; }
 
  inline double RectangleBounds::r() const { 
        return sqrt(m_boundValues[RectangleBounds::bv_halfX]*m_boundValues[RectangleBounds::bv_halfX] 
                  + m_boundValues[RectangleBounds::bv_halfY]*m_boundValues[RectangleBounds::bv_halfY]); 
  }
  
  inline const std::vector< Vector2D > RectangleBounds::vertices() const {
      // create the return vector
      std::vector< Vector2D > vertices;
      // fill the vertices
      vertices.reserve(4);
      vertices.push_back(Vector2D(m_boundValues[RectangleBounds::bv_halfX],-m_boundValues[RectangleBounds::bv_halfY]));   // [0]
      vertices.push_back(Vector2D(m_boundValues[RectangleBounds::bv_halfX], m_boundValues[RectangleBounds::bv_halfY]));   // [1]
      vertices.push_back(Vector2D(-m_boundValues[RectangleBounds::bv_halfX],m_boundValues[RectangleBounds::bv_halfY]));   // [1]
      vertices.push_back(Vector2D(-m_boundValues[RectangleBounds::bv_halfX],-m_boundValues[RectangleBounds::bv_halfY]));  // [3]
      return vertices;
      
  }

} // end of namespace

#endif // ACTS_SURFACES_RECTANGLEBOUNDS_H
