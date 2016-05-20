// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TriangleBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACESTRIANGLEBOUNDS_H
#define ACTS_SURFACESTRIANGLEBOUNDS_H

#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include <utility>

namespace Acts {

   /**
    @class TriangleBounds

    Bounds for a triangular, planar surface.

    @image html TriangularBounds.gif


  class TriangleBounds : public PlanarBounds {


    public:
      /** @enum BoundValues for readability */
      enum BoundValues {
          bv_x1     = 0,
          bv_y1     = 1,
          bv_x2     = 2,
          bv_y2     = 3,
          bv_x3     = 4,
          bv_y3     = 5,
          bv_length = 6
      };

      /**Default Constructor - needed for persistency*/
      TriangleBounds();

      /**Constructor with coordinates of vertices - floats*/
      TriangleBounds( std::vector< std::pair<float,float> >  );

      /**Constructor with coordinates of vertices - double*/
      TriangleBounds( std::vector< std::pair<double,double> >  );

      /**Constructor from three 2 Vectors */
      TriangleBounds( const Vector2D& p1, const Vector2D& p2, const Vector2D& p3);

      /**Copy constructor*/
      TriangleBounds(const TriangleBounds& tribo);

      /**Destructor*/
      virtual ~TriangleBounds();

      /**Assignment Operator*/
      TriangleBounds& operator=(const TriangleBounds& recbo);

      /**Equality operator*/
      virtual bool operator==(const SurfaceBounds& sbo) const override;

      /**Virtual constructor*/
      virtual TriangleBounds* clone() const override;

      /** Return the type of the bounds for persistency */
      virtual BoundsType type() const override { return SurfaceBounds::Triangle; }

      /**This method checks if the provided local coordinates are inside the surface bounds*/
      virtual bool inside(const Vector2D &locpo, double tol1 = 0., double tol2 = 0.) const override;
      virtual bool inside(const Vector2D& locpo, const BoundaryCheck& bchk) const override;

      /** This method checks inside bounds in loc1
        - loc1/loc2 correspond to the natural coordinates of the surface */
      virtual bool insideLoc1(const Vector2D& locpo, double tol1=0.) const override;

      /** This method checks inside bounds in loc2
        - loc1/loc2 correspond to the natural coordinates of the surface */
      virtual bool insideLoc2(const Vector2D& locpo, double tol2=0.) const override;

      /** Minimal distance to boundary ( > 0 if outside and <=0 if inside) */
      virtual double minDistance(const Vector2D& pos) const override;

      /**This method returns the coordinates of vertices */
      const std::vector< Vector2D > vertices() const override;

      /**This method returns the maximal extension on the local plane, i.e. @f$s\sqrt{h_{\phi}^2 + h_{\eta}^2}\f$*/
      virtual double r() const override;

      /** Output Method for std::ostream */
      virtual std::ostream& dump(std::ostream& sl) const override;

    private:
      std::vector<TDD_real_t> m_boundValues;

  };

  inline TriangleBounds* TriangleBounds::clone() const
    { return new TriangleBounds(*this); }

  inline bool TriangleBounds::inside(const Vector2D &locpo, double tol1, double tol2) const {
    std::pair<double,double> locB(m_boundValues.at(TriangleBounds::bv_x2)-m_boundValues.at(TriangleBounds::bv_x1),
				  m_boundValues.at(TriangleBounds::bv_y2)-m_boundValues.at(TriangleBounds::bv_y1));
    std::pair<double,double> locT(m_boundValues.at(TriangleBounds::bv_x3) -locpo[0],
				  m_boundValues.at(TriangleBounds::bv_y3)-locpo[1]);
    std::pair<double,double> locV(m_boundValues.at(TriangleBounds::bv_x1) -locpo[0],
				  m_boundValues.at(TriangleBounds::bv_y1)-locpo[1]);

    // special case :: third vertex ?
    if (locT.first*locT.first+locT.second*locT.second<tol1*tol1) return true;

    // special case : lies on base ?
    double db = locB.first*locV.second-locB.second*locV.first;
    if ( fabs(db)<tol1 ) {
      double a = (locB.first!=0) ? -locV.first/locB.first : -locV.second/locB.second;
      if ( a>-tol2 && a-1.<tol2 ) return true;
      return false;
    }

    double dn = locB.first*locT.second-locB.second*locT.first;

    if ( fabs(dn) > fabs(tol1) ) {
      double t = (locB.first*locV.second-locB.second*locV.first)/dn;
      if ( t > 0.) return false;

      double a = (locB.first!=0.) ? (t*locT.first - locV.first)/locB.first :
                                    (t*locT.second - locV.second)/locB.second ;
      if ( a < -tol2 || a-1.>tol2 ) return false;
    } else {
      return false;
    }
    return true;
  }

  inline bool TriangleBounds::inside(const Vector2D& locpo, const BoundaryCheck& bchk) const
  {
	if (bchk.bcType==0)	return TriangleBounds::inside(locpo, bchk.toleranceLoc1, bchk.toleranceLoc2);

	// a fast FALSE
	double fabsR = sqrt(locpo[Acts::eLOC_X]*locpo[Acts::eLOC_X]+locpo[Acts::eLOC_Y]*locpo[Acts::eLOC_Y]);
	double max_ell = bchk.lCovariance(0,0) > bchk.lCovariance(1,1) ? bchk.lCovariance(0,0) :bchk.lCovariance(1,1);
	double limit = bchk.nSigmas*sqrt(max_ell);
	double r_max = TriangleBounds::r();
	if (fabsR > ( r_max + limit)) return false;

	// compute KDOP and axes for surface polygon
    std::vector<KDOP> elementKDOP(3);
    std::vector<Vector2D> elementP(3);
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
    p << m_boundValues.at(TriangleBounds::bv_x1),m_boundValues.at(TriangleBounds::bv_y1);
    elementP.at(0) =( rotMatrix * (p - locpo) );
    p << m_boundValues.at(TriangleBounds::bv_x2),m_boundValues.at(TriangleBounds::bv_y2);
    elementP.at(1) =( rotMatrix * (p - locpo) );
    p << m_boundValues.at(TriangleBounds::bv_x3),m_boundValues.at(TriangleBounds::bv_y3);
    elementP.at(2) =( rotMatrix * (p - locpo) );
    std::vector<Vector2D> axis = {normal*(elementP.at(1)-elementP.at(0)), normal*(elementP.at(2)-elementP.at(1)), normal*(elementP.at(2)-elementP.at(0))};
    bchk.ComputeKDOP(elementP, axis, elementKDOP);
	// compute KDOP for error ellipse
    std::vector<KDOP> errelipseKDOP(3);
	bchk.ComputeKDOP(bchk.EllipseToPoly(3), axis, errelipseKDOP);
	// check if KDOPs overlap and return result
	return bchk.TestKDOPKDOP(elementKDOP, errelipseKDOP);
  }

  inline bool TriangleBounds::insideLoc1(const Vector2D &locpo, double tol1) const
    { return inside(locpo,tol1,tol1); }

  inline bool TriangleBounds::insideLoc2(const Vector2D &locpo, double tol2) const
    { return inside(locpo,tol2,tol2); }
  
  inline const std::vector< Vector2D > TriangleBounds::vertices() const 
  { 
    std::vector< Vector2D > vertices;
    vertices.resize(3);
    for (size_t iv = 0; iv < 3 ; iv++) 
       vertices.push_back(Vector2D(m_boundValues.at(2*iv),m_boundValues.at(2*iv+1)));
   return vertices;  
  }

  inline double TriangleBounds::r() const {
    double rmax = 0.;
    for (size_t iv = 0; iv < 3 ; iv++)
      rmax = fmax(rmax, m_boundValues.at(2*iv)*m_boundValues.at(2*iv) + m_boundValues.at(2*iv+1)*m_boundValues.at(2*iv+1));
    return sqrt(rmax);
  }

} // end of namespace

#endif // ACTS_SURFACESRECTANGLEBOUNDS_H
