// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// EllipseBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_ELLIPSEBOUNDS_H
#define ACTS_SURFACES_ELLIPSEBOUNDS_H 1

#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include <math.h>
#include <stdlib.h>

namespace Acts {

  /**
   @class EllipseBounds

   Class to describe the bounds for a planar EllipseSurface,
   i.e. the surface between two ellipses.
   By providing an argument for hphisec, the bounds can
   be restricted to a phirange around the center position.

   */

  class EllipseBounds : public PlanarBounds {
    public:
      /** @enum for readibility */
      enum BoundValues {
        bv_rMinX          = 0,
        bv_rMinY          = 1,
        bv_rMaxX          = 2,
        bv_rMaxY          = 3,
        bv_averagePhi     = 4,
        bv_halfPhiSector  = 5,
        bv_length         = 6
      };

      /** Default Constructor*/
      EllipseBounds();

      /** Constructor for full of an ellipsoid disc around phi=0*/
      EllipseBounds(double minrad1, double minrad2, double maxrad1, double maxrad2, double hphisec=M_PI);

      /** Constructor for an ellipsoid disc around phi != 0*/
      EllipseBounds(double minrad1, double minrad2, double maxrad1, double maxrad2, double avephi, double hphisec);

      /** Copy constructor*/
      EllipseBounds(const EllipseBounds& discbo);

      /** Destructor*/
      virtual ~EllipseBounds();

      /** Assignment operator*/
      EllipseBounds& operator=(const EllipseBounds& discbo);

      /** Move assignment operator*/
      EllipseBounds& operator=(EllipseBounds&& discbo);

      /** Equality operator*/
      virtual bool operator==(const SurfaceBounds& sbo) const override;

      /** Virtual constructor*/
      virtual EllipseBounds* clone() const override;

      /** Return the type of the bounds for persistency */
      virtual BoundsType type() const override { return SurfaceBounds::Ellipse; }

      /** This method checks if the point given in the local coordinates is between two ellipsoids
         if only tol1 is given and additional in the phi sector is tol2 is given */
      virtual bool inside(const Vector2D &locpo, double tol1 = 0., double tol2 = 0.) const override;
      virtual bool inside(const Vector2D& locpo, const BoundaryCheck& bchk) const override;

      /** Check for inside first local coordinate */
      virtual bool insideLoc1(const Vector2D &locpo, double tol1 = 0.) const override;

      /** Check for inside second local coordinate */
      virtual bool insideLoc2(const Vector2D &locpo, double tol2 = 0.) const override;

      /** Minimal distance to boundary ( > 0 if outside and <=0 if inside) */
      virtual double minDistance(const Vector2D& pos) const override;

      /** This method returns first inner radius*/
      double rMinX() const;

      /** This method returns second inner radius*/
      double rMinY() const;

      /** This method returns first outer radius*/
      double rMaxX() const;

      /** This method returns second outer radius*/
      double rMaxY() const;

      /** This method returns the maximum expansion on the plane (=max(rMaxX,rMaxY))*/
      virtual double r() const override;

      /** This method returns the average phi*/
      double averagePhi() const;
      
      /** Return the vertices - or, the points of the extremas */
      virtual const std::vector< Vector2D > vertices() const override;

      /** This method returns the halfPhiSector which is covered by the disc*/
      double halfPhiSector() const;

      /** Output Method for std::ostream */
      virtual std::ostream& dump(std::ostream& sl) const override;

    private:
      /** The internal storage of the bounds can be float/double*/
      std::vector<TDD_real_t>           m_boundValues;

      /** helper function for squaring */
      inline double square(double x) const {return x*x;};
 };

 inline EllipseBounds* EllipseBounds::clone() const
 { return new EllipseBounds(*this); }

 inline bool EllipseBounds::inside(const Vector2D &locpo, double tol1, double tol2) const
 { double alpha = acos(cos(locpo[Acts::eLOC_PHI]-m_boundValues.at(EllipseBounds::bv_averagePhi)));
   bool insidePhi = ( m_boundValues.at(EllipseBounds::bv_halfPhiSector)+tol2 < M_PI ) ? (alpha <= (m_boundValues.at(EllipseBounds::bv_halfPhiSector)+tol2)) : 1;
   bool insideInner =( m_boundValues.at(EllipseBounds::bv_rMinX)<=tol1 ) || ( m_boundValues.at(EllipseBounds::bv_rMinY) <= tol1 )
                   || ( square(locpo[Acts::eLOC_X]/(m_boundValues.at(EllipseBounds::bv_rMinX)-tol1)) + square(locpo[Acts::eLOC_Y]/(m_boundValues.at(EllipseBounds::bv_rMinY)-tol1)) >1);
   bool insideOuter = ( square(locpo[Acts::eLOC_X]/(m_boundValues.at(EllipseBounds::bv_rMaxX)+tol1)) + square(locpo[Acts::eLOC_Y]/(m_boundValues.at(EllipseBounds::bv_rMaxY)+tol1)) <1);
   return ( insideInner && insideOuter && insidePhi );

 }

 inline bool EllipseBounds::inside(const Vector2D& locpo, const BoundaryCheck& bchk) const
  {
	return EllipseBounds::inside(locpo,bchk.toleranceLoc1,bchk.toleranceLoc2);
  }

 inline bool EllipseBounds::insideLoc1(const Vector2D &locpo, double tol1) const
 {
   bool insideInner =( m_boundValues[EllipseBounds::bv_rMinX]<=tol1 ) || ( m_boundValues.at(EllipseBounds::bv_rMinY) <= tol1 )
                   || ( square(locpo[Acts::eLOC_X]/(m_boundValues.at(EllipseBounds::bv_rMinX)-tol1)) + square(locpo[Acts::eLOC_Y]/(m_boundValues.at(EllipseBounds::bv_rMinY)-tol1)) >1);
   bool insideOuter = ( square(locpo[Acts::eLOC_X]/(m_boundValues.at(EllipseBounds::bv_rMaxX)+tol1)) + square(locpo[Acts::eLOC_Y]/(m_boundValues.at(EllipseBounds::bv_rMaxY)+tol1)) <1);
   return ( insideInner && insideOuter );
}

 inline bool EllipseBounds::insideLoc2(const Vector2D &locpo, double tol2) const
 {
   double alpha = acos(cos(locpo[Acts::eLOC_PHI]-m_boundValues.at(EllipseBounds::bv_averagePhi)));
   bool insidePhi = ( m_boundValues.at(EllipseBounds::bv_halfPhiSector)+tol2 < M_PI ) ? (alpha <= (m_boundValues.at(EllipseBounds::bv_halfPhiSector)+tol2)) : 1;
   return insidePhi;
 }

 inline double EllipseBounds::rMinX() const { return m_boundValues.at(EllipseBounds::bv_rMinX); }
 inline double EllipseBounds::rMinY() const { return m_boundValues.at(EllipseBounds::bv_rMinY); }
 inline double EllipseBounds::rMaxX() const { return m_boundValues.at(EllipseBounds::bv_rMaxX); }
 inline double EllipseBounds::rMaxY() const { return m_boundValues.at(EllipseBounds::bv_rMaxY); }
 inline double EllipseBounds::r() const { return std::max(m_boundValues.at(EllipseBounds::bv_rMaxX),m_boundValues.at(EllipseBounds::bv_rMaxY)); }
 inline double EllipseBounds::averagePhi() const { return m_boundValues.at(EllipseBounds::bv_averagePhi); }
 inline double EllipseBounds::halfPhiSector() const { return m_boundValues.at(EllipseBounds::bv_halfPhiSector); }
 
 inline const std::vector< Vector2D > EllipseBounds::vertices() const {
     // create the return vector
     std::vector< Vector2D > vertices;
     // fill the vertices
     vertices.reserve(4);
     vertices.push_back(Vector2D(m_boundValues.at(EllipseBounds::bv_rMaxX),0.));   // [0]
     vertices.push_back(Vector2D(0., m_boundValues.at(EllipseBounds::bv_rMaxY)));   // [1]
     vertices.push_back(Vector2D(-m_boundValues.at(EllipseBounds::bv_rMaxX),0.));  // [2]
     vertices.push_back(Vector2D(0., -m_boundValues.at(EllipseBounds::bv_rMaxY)));  // [3]
     return vertices;
     
 }

} // end of namespace

#endif // ACTS_SURFACES_ELLIPSEBOUNDS_H
