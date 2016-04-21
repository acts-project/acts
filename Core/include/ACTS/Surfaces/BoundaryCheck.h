///////////////////////////////////////////////////////////////////
// BoundaryCheck.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_BOUNDARYCHECK_H
#define ACTS_SURFACES_BOUNDARYCHECK_H 1

// STL include(s)
#include <cmath>

// Core module
#include "ACTS/Utilities/ParameterDefinitions.h"
#include "ACTS/Utilities/AlgebraDefinitions.h"

namespace Acts {

    /**
     @class BoundaryCheck

     The BoundaryCheck class allows to steer the way surface
     boundaries are used for inside/outside checks of
     parameters.

     These checks are performed in the local 2D frame of the
     surface and can either be:
     - inside/outside with and without tolerance
     - inside/outside according to a given chi2 value

	 It also provides all the necessary tools for the individual implementations in the different SurfaceBounds classes.

     @author Andreas.Salzburger@cern.ch, Roland.Jansky@cern.ch
     */

	// maint struct for comparing the extent of arbitrary convex polygons relative to prefixed axes
    struct KDOP {
       	float min;
       	// Minimum distance (from origin) along axes
        float max;
       	// Maximum distance (from origin) along axes
    };

	// struct needed for FastSinCos method (see below)
	struct sincosCache {
       	double sinC;
        double cosC;
    };

    class BoundaryCheck {
	  // saves us a lot of function calls in the EllipseToPoly method
	  static constexpr double s_cos22 = 0.923879532511286756128183189396788286822416625863642486115097;
	  static constexpr double s_cos45 = 0.707106781186547524400844362104849039284835937688474036588339;
	  static constexpr double s_cos67 = 0.382683432365089771728459984030398866761344562485627041433800;

      public:
        enum BoundaryCheckType {
            absolute  = 0,  //!< absolute check including tolerances
            chi2corr  = 1   //!< relative (chi2 based) with full correlations
        };

        bool                checkLoc1;          //!< check local 1 coordinate
        bool                checkLoc2;          //!< check local 2 coordinate

        double              toleranceLoc1;      //!< absolute tolerance in local 1 coordinate
        double              toleranceLoc2;      //!< absolute tolerance in local 2 coordinate

        int                 nSigmas;            //!< allowed sigmas for chi2 boundary check
        ActsSymMatrixD<2>     lCovariance;        //!< local covariance matrix

        BoundaryCheckType   bcType;

        /** Constructor for single boolean behavious */
        BoundaryCheck(bool sCheck) :
          checkLoc1(sCheck),
          checkLoc2(sCheck),
          toleranceLoc1(0.),
          toleranceLoc2(0.),
          nSigmas(-1),
          lCovariance(ActsSymMatrixD<2>::Identity()),
          bcType(absolute)
        {}

        /** Constructor for tolerance based check */
        BoundaryCheck(bool chkL1, bool chkL2, double tloc1=0., double tloc2=0.) :
          checkLoc1(chkL1),
          checkLoc2(chkL2),
          toleranceLoc1(tloc1),
          toleranceLoc2(tloc2),
          nSigmas(-1),
          lCovariance(ActsSymMatrixD<2>::Identity()),
          bcType(absolute)
        {}

        /** Constructor for chi2 based check */
        BoundaryCheck(const ActsSymMatrixD<2>& lCov,
					  int nsig=1,
                      bool chkL1=true,
                      bool chkL2=true) :
          checkLoc1(chkL1),
          checkLoc2(chkL2),
          toleranceLoc1(0.),
          toleranceLoc2(0.),
          nSigmas(nsig),
          lCovariance(lCov),
          bcType(chi2corr)
        {}

        /** Conversion operator to bool */
        operator bool() const { return (checkLoc1 || checkLoc2); }


      	/** Each Bounds has a method inside, which checks if a LocalPosition is inside the bounds.
        	  Inside can be called without/with boundary check */
      	void ComputeKDOP(std::vector<Vector2D> v, std::vector<Vector2D> KDOPAxes, std::vector<KDOP> &kdop) const;

      	std::vector<Vector2D> EllipseToPoly(int resolution = 3) const;

      	bool TestKDOPKDOP(std::vector<KDOP> &a, std::vector<KDOP> &b) const;

		double FastArcTan(double x) const;

		sincosCache FastSinCos(double x) const;
    };

	// should have maximum (average) error of 0.0015 (0.00089) radians or 0.0859 (0.0509) degrees, fine for us and much faster (>4 times)
  	inline double BoundaryCheck::FastArcTan(double x) const
	{
		double y;
		bool complement = false; // true if arg was >1
		bool sign = false; // true if arg was < 0
		if (x<0.){
			x = -x;
			sign = true; // arctan(-x)=-arctan(x)
		}
		if (x>1.){
			x = 1./x; // keep arg between 0 and 1
			complement = true;
		}
		y = M_PI_4*x - x*(fabs(x) - 1)*(0.2447 + 0.0663*fabs(x));
		if (complement) y = M_PI_2 - y; // correct for 1/x if we did that
		if (sign) y = -y; // correct for negative arg
		return y;
	}

	// should have maximum (average) error of 0.001 (0.0005) radians or 0.0573 (0.029) degrees, fine for us and much faster (>8 times)
	inline sincosCache BoundaryCheck::FastSinCos(double x) const
	{
		sincosCache tmp;
		//always wrap input angle to -PI..PI
		if (x < -M_PI)
			x += 2.*M_PI;
		else
		if (x >  M_PI)
			x -= 2.*M_PI;

		//compute sine
		if (x < 0.)
		{
			tmp.sinC = 1.27323954 * x + .405284735 * x * x;

			if (tmp.sinC < 0.)
				tmp.sinC = .225 * (tmp.sinC *-tmp.sinC - tmp.sinC) + tmp.sinC;
			else
				tmp.sinC = .225 * (tmp.sinC * tmp.sinC - tmp.sinC) + tmp.sinC;
		}
		else
		{
			tmp.sinC = 1.27323954 * x - 0.405284735 * x * x;

			if (tmp.sinC < 0.)
				tmp.sinC = .225 * (tmp.sinC *-tmp.sinC - tmp.sinC) + tmp.sinC;
			else
				tmp.sinC = .225 * (tmp.sinC * tmp.sinC - tmp.sinC) + tmp.sinC;
		}

		//compute cosine: sin(x + PI/2) = cos(x)
		x += M_PI_2;
		if (x >  M_PI)
			x -= 2.*M_PI;

		if (x < 0.)
		{
			tmp.cosC = 1.27323954 * x + 0.405284735 * x * x;

			if (tmp.cosC < 0.)
				tmp.cosC = .225 * (tmp.cosC *-tmp.cosC - tmp.cosC) + tmp.cosC;
			else
				tmp.cosC = .225 * (tmp.cosC * tmp.cosC - tmp.cosC) + tmp.cosC;
		}
		else
		{
			tmp.cosC = 1.27323954 * x - 0.405284735 * x * x;

			if (tmp.cosC < 0.)
				tmp.cosC = .225 * (tmp.cosC *-tmp.cosC - tmp.cosC) + tmp.cosC;
			else
				tmp.cosC = .225 * (tmp.cosC * tmp.cosC - tmp.cosC) + tmp.cosC;
		}
		return tmp;
	}

	// does the conversion of an ellipse of height h and width w to an polygon with 4 + 4*resolution points
  	inline std::vector<Vector2D> BoundaryCheck::EllipseToPoly(int resolution) const
	{
		const double h = nSigmas*sqrt(lCovariance(1,1));
		const double w = nSigmas*sqrt(lCovariance(0,0));

		// first add the four vertices
		std::vector<Vector2D> v((1+resolution)*4);
		Vector2D p;
		p << w,	0;	v[0] = p;
		p << -w,0;	v[1] = p;
		p << 0,	h;	v[2] = p;
		p << 0,-h;	v[3] = p;

		// now add a number, equal to the resolution, of equidistant points  in each quadrant
		// resolution == 3 seems to be a solid working point, but possibility open to change to any number in the future
		Vector2D t(1,1);
		ActsSymMatrixD<2> t1; t1 << 1,0,0,-1;
		ActsSymMatrixD<2> t2;	t2 << -1,0,0,-1;
		ActsSymMatrixD<2> t3;	t3 << -1,0,0,1;
		if(resolution != 3){
			sincosCache scResult;
			for(int i=1; i<=resolution; i++) {
			    scResult = FastSinCos(M_PI_2*i/(resolution+1));
				t << w*scResult.sinC,h*scResult.cosC;
				v[i*4+0] = t;
				v[i*4+1] = t1*t;
				v[i*4+2] = t2*t;
				v[i*4+3] = t3*t;
			}
		}
		else{
			t << w*s_cos22,h*s_cos67;
			v[4] = t;
			v[5] = t1*t;
			v[6] = t2*t;
			v[7] = t3*t;
			t << w*s_cos45,h*s_cos45;
			v[8] = t;
			v[9] = t1*t;
			v[10] = t2*t;
			v[11] = t3*t;
			t << w*s_cos67,h*s_cos22;
			v[12] = t;
			v[13] = t1*t;
			v[14] = t2*t;
			v[15] = t3*t;
		}
	return v;
	}

	// calculates KDOP object from given polygon and set of axes
    inline void BoundaryCheck::ComputeKDOP(std::vector<Vector2D> v, std::vector<Vector2D> KDOPAxes, std::vector<KDOP> &kdop) const
	{
		// initialize KDOP to first point
		unsigned int k = KDOPAxes.size();
		for(unsigned int i=0; i<k; i++) {
			kdop[i].max = KDOPAxes[i](0,0)*v[0](0,0)+KDOPAxes[i](1,0)*v[0](1,0);
			kdop[i].min = KDOPAxes[i](0,0)*v[0](0,0)+KDOPAxes[i](1,0)*v[0](1,0);
		}
		// now for each additional point, update KDOP bounds if necessary
		float value;
		for(unsigned int i=1; i<v.size(); i++) {
		    for(unsigned int j=0; j<k; j++) {
				value = KDOPAxes[j](0,0)*v[i](0,0)+KDOPAxes[j](1,0)*v[i](1,0);
				if (value < kdop[j].min) kdop[j].min = value;
				else if (value > kdop[j].max) kdop[j].max = value;
			}
		}
	}

	// this is the method to actually check if two KDOPs overlap
	inline bool BoundaryCheck::TestKDOPKDOP(std::vector<KDOP> &a, std::vector<KDOP> &b) const
	{
		int k = a.size();
		// check if any intervals are non-overlapping, return if so
		for(int i=0;i<k;i++) if(a[i].min > b[i].max || a[i].max < b[i].min) return false;
		// all intervals are overlapping, so KDOPs must intersect
		return true;
	}

}

#endif // ACTS_SURFACES_BOUNDARYCHECK_H
