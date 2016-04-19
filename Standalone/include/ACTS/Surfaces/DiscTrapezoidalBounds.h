///////////////////////////////////////////////////////////////////
// DiscTrapezoidalBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_DISCTRAPEZOIDALBOUNDS_H
#define ACTS_SURFACES_DISCTRAPEZOIDALBOUNDS_H

// Geometry module
#include "ACTS/Surfaces/DiscBounds.h"
#include "ACTS/GeometryUtils/PrecisionDefinition.h"
// EventData module
#include "ACTS/Utilities/AlgebraDefinitions.h"
// STD/STL
#include <math.h>
#include <cmath>
#include "ACTS/Utilities/ParameterDefinitions.h"

namespace Acts {

  /**
     @class DiscTrapezoidalBounds

     Class to describe the bounds for a planar DiscSurface.
     By providing an argument for hphisec, the bounds can
     be restricted to a phirange around the center position.

     @author Noemi.Calace@cern.ch , Andreas.Salzburger@cern.ch
  */

  class DiscTrapezoidalBounds : public DiscBounds {

  public:
    /** enumeration for readability */
    enum BoundValues {
      bv_rMin          = 0,
      bv_rMax          = 1,
      bv_minHalfX      = 2,
      bv_maxHalfX      = 3,
      bv_halfY         = 4,
      bv_halfPhiSector = 5,
      bv_averagePhi    = 6,
      bv_rCenter       = 7,
      bv_stereo        = 8,
      bv_length        = 9
    };

    /** Default Constructor*/
    DiscTrapezoidalBounds();

    /** Constructor for a symmetric Trapezoid giving min X lenght, max X lenght, Rmin and R max */
    DiscTrapezoidalBounds(double minhalfx, double maxhalfx, double rMin, double rMax, double avephi, double stereo = 0.);

    /** Copy constructor*/
    DiscTrapezoidalBounds(const DiscTrapezoidalBounds& disctrbo);

    /** Destructor*/
    virtual ~DiscTrapezoidalBounds();

    /** Assignment operator*/
    DiscTrapezoidalBounds& operator=(const DiscTrapezoidalBounds& disctrbo);

    /** Equality operator*/
    bool operator==(const SurfaceBounds& sbo) const;

    /** Virtual constructor*/
    virtual DiscTrapezoidalBounds* clone() const override;

    /** Return the type - mainly for persistency */
    virtual SurfaceBounds::BoundsType type() const override { return SurfaceBounds::DiscTrapezoidal; }

    /** This method cheks if the radius given in the LocalPosition is inside [rMin,rMax]
       if only tol1 is given and additional in the phi sector is tol2 is given */
    virtual bool inside(const Vector2D &locpo, double tol1=0., double tol2=0.) const override;
    virtual bool inside(const Vector2D& locpo, const BoundaryCheck& bchk = true) const override;

    /** This method checks inside bounds in loc1
	- loc1/loc2 correspond to the natural coordinates of the surface */
    virtual bool insideLoc1(const Vector2D& locpo, double tol1=0.) const override;

    /** This method checks inside bounds in loc2
	- loc1/loc2 correspond to the natural coordinates of the surface */
    virtual bool insideLoc2(const Vector2D& locpo, double tol2=0.) const override;

    /** Minimal distance to boundary ( > 0 if outside and <=0 if inside) */
    virtual double minDistance(const Vector2D& pos) const override;

    /** This method returns inner radius*/
    double rMin() const;

    /** This method returns outer radius*/
    double rMax() const;

    /** This method returns the maximum expansion on the plane (=rMax)*/
    virtual double r() const override;

    /** This method returns the average phi*/
    double averagePhi() const;

    /** This method returns the center radius*/
    double rCenter() const;

    /** This method returns the stereo angle*/
    double stereo() const;

    /** This method returns the halfPhiSector which is covered by the disc*/
    double halfPhiSector() const;

    /** This method returns the minimal halflength in X */
    double minHalflengthX() const;

    /** This method returns the maximal halflength in X */
    double maxHalflengthX() const;

    /** This method returns the halflength in Y (this is Rmax -Rmin)*/
    double halflengthY() const;

    /** Output Method for std::ostream */
    virtual std::ostream& dump(std::ostream& sl) const override;

  private:
    /** Internal members of the bounds (float/double)*/
    std::vector<TDD_real_t>                    m_boundValues;

  };

  inline DiscTrapezoidalBounds* DiscTrapezoidalBounds::clone() const
  { return new DiscTrapezoidalBounds(*this); }

  inline bool DiscTrapezoidalBounds::inside(const Vector2D &locpo, double , double ) const
  {
    double rMedium = m_boundValues[DiscTrapezoidalBounds::bv_rCenter];
    double phi     = m_boundValues[DiscTrapezoidalBounds::bv_averagePhi];

    Vector2D cartCenter(rMedium*cos(phi),
    			     rMedium*sin(phi));

    Vector2D cartPos(locpo[Acts::eLOC_R]*cos(locpo[Acts::eLOC_PHI]),
    			  locpo[Acts::eLOC_R]*sin(locpo[Acts::eLOC_PHI]));

    Vector2D Pos = cartPos - cartCenter;

    Vector2D DeltaPos(Pos[Acts::eLOC_X]*sin(phi) - Pos[Acts::eLOC_Y]*cos(phi),
    			   Pos[Acts::eLOC_Y]*sin(phi) + Pos[Acts::eLOC_X]*cos(phi));

    bool insideX = (fabs(DeltaPos[Acts::eLOC_X])<m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]);
    bool insideY = (fabs(DeltaPos[Acts::eLOC_Y])<m_boundValues[DiscTrapezoidalBounds::bv_halfY]);

    if (!insideX || !insideY) return false;

    if (fabs(DeltaPos[Acts::eLOC_X]) < m_boundValues[DiscTrapezoidalBounds::bv_minHalfX]) return true;

    double m = 2.*m_boundValues[DiscTrapezoidalBounds::bv_halfY]/(m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX] - m_boundValues[DiscTrapezoidalBounds::bv_minHalfX]);

    double q = m_boundValues[DiscTrapezoidalBounds::bv_halfY]*(m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX] + m_boundValues[DiscTrapezoidalBounds::bv_minHalfX])/(m_boundValues[DiscTrapezoidalBounds::bv_minHalfX] - m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]);

    bool inside = (DeltaPos[Acts::eLOC_Y] > (m*fabs(DeltaPos[Acts::eLOC_X])+q));

    return inside;

  }

  inline bool DiscTrapezoidalBounds::inside(const Vector2D& locpo, const BoundaryCheck& bchk) const
    {
      if(bchk.bcType==0 || bchk.nSigmas==0 || m_boundValues[DiscTrapezoidalBounds::bv_rMin]!=0)	return DiscTrapezoidalBounds::inside(locpo, bchk.toleranceLoc1, bchk.toleranceLoc2);
      double alpha = fabs(locpo[Acts::eLOC_PHI]-m_boundValues[DiscTrapezoidalBounds::bv_averagePhi]);
      if ( alpha>M_PI) alpha = 2*M_PI - alpha;
      // a fast FALSE
      sincosCache scResult = bchk.FastSinCos(locpo(1,0));
      double dx = bchk.nSigmas*sqrt(bchk.lCovariance(0,0));
      double dy = bchk.nSigmas*sqrt(scResult.sinC*scResult.sinC*bchk.lCovariance(0,0) + locpo(0,0)*locpo(0,0)*scResult.cosC*scResult.cosC*bchk.lCovariance(1,1)+2*scResult.cosC*scResult.sinC*locpo(0,0)*bchk.lCovariance(0,1));
      double max_ell = dx > dy ? dx : dy;
      if (locpo(0,0) > ( m_boundValues[DiscTrapezoidalBounds::bv_rMax]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha) + max_ell))
	return false;
      // a fast TRUE
      double min_ell = dx < dy ? dx : dy;
      if (locpo(0,0) < ( m_boundValues[DiscTrapezoidalBounds::bv_rMax]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha) + min_ell))
	return true;

      // we are not using the KDOP approach here but rather a highly optimized one
      class EllipseCollisionTest {
      private:
	int m_maxIterations;
	bool iterate(double x, double y, double c0x, double c0y, double c2x, double c2y, double rr) const {
	  std::vector<double> innerPolygonCoef(m_maxIterations+1);
	  std::vector<double> outerPolygonCoef(m_maxIterations+1);
		/*
		   t2______t4
               --_     	\
                    --_	 \	              				    /¨¨¨ ¨¨\
                    	      t1 = (0, 0)   	           (     t   )
                     		  | \          					    \__ _ /
                        	  |   \
                  	          |    t3
                    	      |   /
                     	      | /
                           t0
		*/
	  for (int t = 1; t <= m_maxIterations; t++) {
	    int numNodes = 4 << t;
	    innerPolygonCoef[t] = 0.5/cos(4*acos(0.0)/numNodes);
	    double c1x = (c0x + c2x)*innerPolygonCoef[t];
	    double c1y = (c0y + c2y)*innerPolygonCoef[t];
	    double tx = x - c1x; // t indicates a translated coordinate
	    double ty = y - c1y;
	    if (tx*tx + ty*ty <= rr) {
	      return true;	// collision with t1
	    }
	    double t2x = c2x - c1x;
	    double t2y = c2y - c1y;
	    if (tx*t2x + ty*t2y >= 0 && tx*t2x + ty*t2y <= t2x*t2x + t2y*t2y &&
		(ty*t2x - tx*t2y >= 0 || rr*(t2x*t2x + t2y*t2y) >= (ty*t2x - tx*t2y)*(ty*t2x - tx*t2y))) {
	      return true;	// collision with t1---t2
		  }
	    double t0x = c0x - c1x;
	    double t0y = c0y - c1y;
	    if (tx*t0x + ty*t0y >= 0 && tx*t0x + ty*t0y <= t0x*t0x + t0y*t0y &&
		(ty*t0x - tx*t0y <= 0 || rr*(t0x*t0x + t0y*t0y) >= (ty*t0x - tx*t0y)*(ty*t0x - tx*t0y))) {
	      return true;	// collision with t1---t0
	    }
	    outerPolygonCoef[t] = 0.5/(cos(2*acos(0.0)/numNodes)*cos(2*acos(0.0)/numNodes));
	    double c3x = (c0x + c1x)*outerPolygonCoef[t];
	    double c3y = (c0y + c1y)*outerPolygonCoef[t];
	    if ((c3x-x)*(c3x-x) + (c3y-y)*(c3y-y) < rr) {
	      c2x = c1x;
	      c2y = c1y;
	      continue;	// t3 is inside circle
	    }
	    double c4x = c1x - c3x + c1x;
	    double c4y = c1y - c3y + c1y;
	    if ((c4x-x)*(c4x-x) + (c4y-y)*(c4y-y) < rr) {
	      c0x = c1x;
	      c0y = c1y;
	      continue;	// t4 is inside circle
	    }
	    double t3x = c3x - c1x;
	    double t3y = c3y - c1y;
	    if (ty*t3x - tx*t3y <= 0 || rr*(t3x*t3x + t3y*t3y) > (ty*t3x - tx*t3y)*(ty*t3x - tx*t3y)) {
	      if (tx*t3x + ty*t3y > 0) {
		if (std::abs(tx*t3x + ty*t3y) <= t3x*t3x + t3y*t3y || (x-c3x)*(c0x-c3x) + (y-c3y)*(c0y-c3y) >= 0) {
		  c2x = c1x;
		  c2y = c1y;
		  continue;	// circle center is inside t0---t1---t3
		}
	      } else if (-(tx*t3x + ty*t3y) <= t3x*t3x + t3y*t3y || (x-c4x)*(c2x-c4x) + (y-c4y)*(c2y-c4y) >= 0) {
		c0x = c1x;
		c0y = c1y;
		continue;	// circle center is inside t1---t2---t4
	      }
	    }
	    return false; // no collision possible
	  }
	  return false; // out of iterations so it is unsure if there was a collision. But have to return something.
	}
      public:
	// test for collision between an ellipse of horizontal radius w and vertical radius h at (x0, y0) and a circle of radius r at (x1, y1)
	bool collide(double x0, double y0, double w, double h, double x1, double y1, double r) const {
	  double x = fabs(x1 - x0);
	  double y = fabs(y1 - y0);
	  if (x*x + (h - y)*(h - y) <= r*r || (w - x)*(w - x) + y*y <= r*r || x*h + y*w <= w*h	       		 // collision with (0, h)
	      || ((x*h + y*w - w*h)*(x*h + y*w - w*h) <= r*r*(w*w + h*h) && x*w - y*h >= -h*h && x*w - y*h <= w*w)) {   // collision with (0, h)---(w, 0)
	    return true;
	  } else {
	    if ((x-w)*(x-w) + (y-h)*(y-h) <= r*r || (x <= w && y - r <= h) || (y <= h && x - r <= w)) {
	      return iterate(x, y, w, 0, 0, h, r*r);	 // collision within triangle (0, h) (w, h) (0, 0) is possible
	    }
	    return false;
	  }
	}
	EllipseCollisionTest(int maxIterations) {
	  this->m_maxIterations = maxIterations;
	}
      };

      EllipseCollisionTest test(4);
      // convert to cartesian coordinates
      ActsMatrixD<2,2> covRotMatrix ;
      covRotMatrix << scResult.cosC, -locpo(0,0)*scResult.sinC,
	scResult.sinC, locpo(0,0)*scResult.cosC;
      ActsMatrixD<2,2> lCovarianceCar = covRotMatrix*bchk.lCovariance*covRotMatrix.transpose();
      Vector2D locpoCar(covRotMatrix(1,1), -covRotMatrix(0,1));

      // ellipse is always at (0,0), surface is moved to ellipse position and then rotated
      double w = bchk.nSigmas*sqrt( lCovarianceCar(0,0));
      double h = bchk.nSigmas*sqrt( lCovarianceCar(1,1));
      double x0 = 0;
      double y0 = 0;
      float theta = (lCovarianceCar(1,0) != 0 && (lCovarianceCar(1,1)-lCovarianceCar(0,0))!=0 ) ? .5*bchk.FastArcTan( 2*lCovarianceCar(1,0)/(lCovarianceCar(1,1)-lCovarianceCar(0,0)) ) : 0.;
      scResult = bchk.FastSinCos(theta);
      ActsMatrixD<2,2> rotMatrix ;
      rotMatrix << scResult.cosC, scResult.sinC,
	-scResult.sinC, scResult.cosC;
      Vector2D tmp = rotMatrix * (- locpoCar) ;
      double x1 = tmp(0,0);
      double y1 = tmp(1,0);
      double r = m_boundValues[DiscTrapezoidalBounds::bv_rMax];
      // check if ellipse and circle overlap and return result
      return test.collide(x0, y0, w, h, x1, y1, r);
    }

  inline bool DiscTrapezoidalBounds::insideLoc1(const Vector2D &locpo, double tol1) const {
    double alpha = fabs(locpo[Acts::eLOC_PHI]-m_boundValues[DiscTrapezoidalBounds::bv_averagePhi]);
    if ( alpha>M_PI) alpha = 2*M_PI - alpha;

    return ( locpo[Acts::eLOC_R] > (m_boundValues[DiscTrapezoidalBounds::bv_rMin]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha) - tol1) &&
	     locpo[Acts::eLOC_R] < (m_boundValues[DiscTrapezoidalBounds::bv_rMax]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha) + tol1)); }

  inline bool DiscTrapezoidalBounds::insideLoc2(const Vector2D &locpo, double tol2) const {
    double alpha = fabs(locpo[Acts::eLOC_PHI]-m_boundValues[DiscTrapezoidalBounds::bv_averagePhi]);
    if ( alpha > M_PI ) alpha  = 2.*M_PI-alpha;
    return (alpha <= (m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector] + tol2));
  }

  inline double DiscTrapezoidalBounds::rMin() const { return m_boundValues[DiscTrapezoidalBounds::bv_rMin]; }

  inline double DiscTrapezoidalBounds::rMax() const { return m_boundValues[DiscTrapezoidalBounds::bv_rMax]; }

  inline double DiscTrapezoidalBounds::r() const { return m_boundValues[DiscTrapezoidalBounds::bv_rMax]; }

  inline double DiscTrapezoidalBounds::averagePhi() const { return m_boundValues[DiscTrapezoidalBounds::bv_averagePhi]; }

  inline double DiscTrapezoidalBounds::rCenter() const { return m_boundValues[DiscTrapezoidalBounds::bv_rCenter]; }

  inline double DiscTrapezoidalBounds::stereo() const { return m_boundValues[DiscTrapezoidalBounds::bv_stereo]; }

  inline double DiscTrapezoidalBounds::halfPhiSector() const { return m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector]; }

  inline double DiscTrapezoidalBounds::minHalflengthX() const { return m_boundValues[DiscTrapezoidalBounds::bv_minHalfX]; }

  inline double DiscTrapezoidalBounds::maxHalflengthX() const { return m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]; }

  inline double DiscTrapezoidalBounds::halflengthY() const {return m_boundValues[DiscTrapezoidalBounds::bv_halfY]; }

} // end of namespace

#endif // ACTS_SURFACES_DISCTRAPEZOIDALBOUNDS_H
