///////////////////////////////////////////////////////////////////
// BevelledCylinderVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BEVELLEDCYLINDERVOLUMESBOUNDS_H
#define ACTS_VOLUMES_BEVELLEDCYLINDERVOLUMESBOUNDS_H

// Geometry module
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Volumes/VolumeBounds.h"
#include "ACTS/Utilities/Definitions.h"
// Core module

namespace Acts {

   /** @struct BevelledBoundaryIntersector
          - identical code to ExtrapolationUtils/RealLinearEquation (but forbidden dependency!)
      */
   struct BevelledBoundaryIntersector {

   double yOfX;            //!< the result of x
   double segLength;      //!< length of the line segment

    BevelledBoundaryIntersector(double px, double py, double k, double xprime) {
      double deltax = xprime-px;
      yOfX = py + k*(deltax);
      double deltay = yOfX-py;
      segLength = sqrt(deltax*deltax+deltay*deltay);
    }

   };

   class Surface;
   class EllipseBounds;
   class TrapezoidBounds;
   class RectangleBounds;
   class CylinderBounds;
   class RadialBounds;

  /**
   @class BevelledCylinderVolumeBounds

   Bounds for a cylindrical Volume, the decomposeToSurfaces method creates a
   vector of up to 6 surfaces:

   case A) 3 Surfaces (full cylindrical tube):
    BoundarySurfaceFace [index]:
        - negativeFaceXY [0] : Acts::DiscSurface with \f$ r_{inner}=0 \f$,
                               parallel to \f$ xy \f$ plane at negative \f$ z \f$
        - positiveFaceXY [1] : Acts::DiscSurface with \f$ r_{inner}=0 \f$,
                               parallel to \f$ xy \f$ plane at positive \f$ z \f$
        - cylinderCover  [2] : Acts::CylinderSurface confining the Acts::Volume

   case B) 4 Surfaces (tube with inner and outer radius):
    BoundarySurfaceFace [index]:
        - negativeFaceXY [0] : Acts::DiscSurface with \f$ r_{inner}>0 \f$,
                               parallel to \f$ xy \f$ plane at negative \f$ z \f$
        - positiveFaceXY [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$,
                               parallel to \f$ xy \f$ plane at positive \f$ z \f$
        - tubeOuterCover [2] : Acts::CylinderSurface with \f$ r = r_{outer} \f$
        - tubeInnerCover [3] : Acts::CylinderSurface with \f$ r = r_{inner} \f$

   case C) 6 Surfaces (sectoral tube with inner and outer radius):
    BoundarySurfaceFace [index]:
        - negativeFaceXY        [0] : Acts::DiscSurface with \f$ r_{inner}>0  \f$ and \f$ \phi < \pi \f$,
                                      parallel to \f$ xy \f$ plane at negative \f$ z \f$
        - positiveFaceXY        [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$ and \f$ \phi < \pi \f$,
                                      parallel to \f$ xy \f$ plane at positive \f$ z \f$
        - tubeSectorOuterCover  [2] : Acts::CylinderSurface with \f$ r = r_{outer} \f$
        - tubeSectorInnerCover  [3] : Acts::CylinderSurface with \f$ r = r_{inner} \f$
        - tubeSectorNegativePhi [4] : Rectangular Acts::PlaneSurface attached to [0] and [1] at negative \f$ \phi \f$
        - tubeSectorNegativePhi [5] : Rectangular Acts::PlaneSurface attached to [0] and [1] at positive \f$ \phi \f$

   case D) 6 Surfaces (sectoral bevelled tube with inner and/or outer radius replaced by plane surface):
    BoundarySurfaceFace [index]:
        - negativeFaceXY        [0] : Acts::DiscSurface with \f$ r_{inner}>0  \f$ and \f$ \phi < \pi \f$,
                                      parallel to \f$ xy \f$ plane at negative \f$ z \f$
        - positiveFaceXY        [1] : Acts::DiscSurface with \f$ r_{inner}>0 \f$ and \f$ \phi < \pi \f$,
                                      parallel to \f$ xy \f$ plane at positive \f$ z \f$
        - tubeSectorOuterCover  [2] : Acts::CylinderSurface with \f$ r = r_{outer} \f$
                                   OR rectangular plane surface
        - tubeSectorInnerCover  [3] : Acts::CylinderSurface with \f$ r = r_{inner} \f$
                                   OR rectangular plane surface
        - tubeSectorNegativePhi [4] : Rectangular Acts::PlaneSurface attached to [0] and [1] at negative \f$ \phi \f$
        - tubeSectorNegativePhi [5] : Rectangular Acts::PlaneSurface attached to [0] and [1] at positive \f$ \phi \f$

    @image html CylinderVolumeBounds_decomp.gif

    @author Andreas.Salzburger@cern.ch
    */

 class BevelledCylinderVolumeBounds : public VolumeBounds {

  public:

    /** @enum BoundValues for readibility */
    enum BoundValues {
          bv_innerRadius    = 0,
          bv_outerRadius    = 1,
          bv_halfPhiSector  = 2,
          bv_halfZ          = 3,
          bv_thetaMinus     = 4,
          bv_thetaPlus      = 5,
          bv_length         = 6
    };

    /**Default Constructor*/
    BevelledCylinderVolumeBounds();

    /**Constructor - cylinder segment bevelled in R */
    BevelledCylinderVolumeBounds(double rinner, double router, double halfPhiSector, double halez, int type);

    /**Copy Constructor */
    BevelledCylinderVolumeBounds(const BevelledCylinderVolumeBounds& cylbo);

    /**Destructor */
    virtual ~BevelledCylinderVolumeBounds();

    /**Assignment operator*/
    BevelledCylinderVolumeBounds& operator=(const BevelledCylinderVolumeBounds& cylbo);

    /**Virtual constructor */
    BevelledCylinderVolumeBounds* clone() const override;

    /**This method checks if position in the 3D volume frame is inside the cylinder*/
    bool inside(const Vector3D& , double tol=0.) const override;

    /** Method to decompose the Bounds into boundarySurfaces */
    const std::vector<const Acts::Surface*>* decomposeToSurfaces(std::shared_ptr<Transform3D> transformPtr) const override;

    /**This method returns the inner radius*/
    double innerRadius() const;

    /**This method returns the outer radius*/
    double outerRadius() const;

    /**This method returns the medium radius*/
    double mediumRadius() const;

    /**This method returns the delta radius*/
    double deltaRadius() const;

    /**This method returns the halfPhiSector angle*/
    double halfPhiSector() const;

    /**This method returns the halflengthZ*/
    double halflengthZ() const;

     /**This method returns the thetaMinus*/
    double thetaMinus() const;

     /**This method returns the thetaPlus*/
    double thetaPlus() const;

     /**This method returns the type*/
    int type() const;

    /** Output Method for std::ostream */
    std::ostream& dump(std::ostream& sl) const override;

  private:
    template <class T> T& dumpT(T& tstream) const;

    /** This method returns the associated BevelledCylinderBounds of the inner BevelledCylinderSurfaces. */
    CylinderBounds* innerBevelledCylinderBounds() const;
    /** This method returns the associated BevelledCylinderBounds of the outer BevelledCylinderSurfaces. */
    CylinderBounds* outerBevelledCylinderBounds() const;
    /** This method returns the associated plane surface bounds of the inner bevelled surface*/
    RectangleBounds* innerBevelledPlaneBounds() const;
    /** This method returns the associated BevelledCylinderBounds of the outer BevelledCylinderSurfaces. */
    RectangleBounds* outerBevelledPlaneBounds() const;
    /** This method returns the associated EllipseBounds for the bottom/top EllipseSurface. */
    EllipseBounds* bottomEllipseBounds() const;
    /** This method returns the associated EllipseBounds for the bottom/top EllipseSurface. */
    EllipseBounds* topEllipseBounds() const;
    /** This method returns the associated CylinderBounds of the inner CylinderSurfaces. */
    CylinderBounds* innerCylinderBounds() const;
    /** This method returns the associated CylinderBounds of the outer CylinderSurfaces. */
    CylinderBounds* outerCylinderBounds() const;
    /** This method returns the associated RadialBounds for the bottom/top DiscSurface. */
    RadialBounds* discBounds() const;
    /** This method returns the bevelled area volume. */
    Volume* subtractedVolume() const;
    /** This method returns the associated PlaneBounds limiting a sectoral BevelledCylinderVolume. */
    TrapezoidBounds* sectorTrdBounds() const;
    RectangleBounds* sectorPlaneBounds() const;

    std::vector<TDD_real_t>               m_boundValues;
    int                                   m_type;

    /** numerical stability */
    static double s_numericalStable;

    mutable Acts::Volume*                  m_subtractedVolume;
 };

 inline BevelledCylinderVolumeBounds* BevelledCylinderVolumeBounds::clone() const
 { return new BevelledCylinderVolumeBounds(*this); }

 inline bool BevelledCylinderVolumeBounds::inside(const Vector3D &pos, double tol) const
 {
   double ros = pos.perp();
   bool insidePhi =  fabs(pos.phi()) <= m_boundValues[bv_halfPhiSector] + tol;
   double cphi = cos(pos.phi());
   bool insideR = insidePhi;
   if (m_type < 1)  insideR  = insidePhi ? ((ros >=  m_boundValues[bv_innerRadius] - tol ) && (ros <= m_boundValues[bv_outerRadius] + tol)) : false;
   else if (m_type == 1 ) insideR = insidePhi ? ((ros>= m_boundValues[bv_innerRadius]/cphi-tol)&&(ros<=m_boundValues[bv_outerRadius]+tol)):false;
   else if (m_type == 2 ) insideR = insidePhi ? ((ros>= m_boundValues[bv_innerRadius]-tol)&&(ros<=m_boundValues[bv_outerRadius]/cphi+tol)):false;
   else if (m_type == 3 ) insideR = insidePhi ? ((ros>= m_boundValues[bv_innerRadius]/cphi-tol)&&(ros<=m_boundValues[bv_outerRadius]/cphi+tol)):false;
//   bool insideZ = insideR ? (fabs(pos.z()) <= m_boundValues[bv_halfZ] + tol ) : false ;
   bool insideZ = insideR ? ((pos.z()<=m_boundValues[bv_halfZ]-(pos.x()+m_boundValues[bv_outerRadius])*tan(m_boundValues[bv_thetaPlus]) + tol)
			     && ( pos.z()>=-m_boundValues[bv_halfZ]+(pos.x()+m_boundValues[bv_outerRadius])*tan(m_boundValues[bv_thetaMinus]) - tol)) : false ;

   return insideZ;
 }

 inline double BevelledCylinderVolumeBounds::innerRadius() const { return m_boundValues[bv_innerRadius]; }
 inline double BevelledCylinderVolumeBounds::outerRadius() const { return m_boundValues[bv_outerRadius]; }
 inline double BevelledCylinderVolumeBounds::mediumRadius() const { return 0.5*(m_boundValues[bv_innerRadius]+m_boundValues[bv_outerRadius]); }
 inline double BevelledCylinderVolumeBounds::deltaRadius() const { return (m_boundValues[bv_outerRadius]-m_boundValues[bv_innerRadius]); }
 inline double BevelledCylinderVolumeBounds::halfPhiSector() const { return m_boundValues[bv_halfPhiSector]; }
 inline double BevelledCylinderVolumeBounds::halflengthZ() const { return m_boundValues[bv_halfZ]; }
 inline double BevelledCylinderVolumeBounds::thetaMinus() const { return m_boundValues[bv_thetaMinus]; }
 inline double BevelledCylinderVolumeBounds::thetaPlus() const { return m_boundValues[bv_thetaPlus]; }
 inline int BevelledCylinderVolumeBounds::type() const { return m_type; }

 template <class T> T& BevelledCylinderVolumeBounds::dumpT(T& tstream) const
 {
     tstream << std::setiosflags(std::ios::fixed);
     tstream << std::setprecision(7);
     tstream << "Acts::BevelledCylinderVolumeBounds: (innerR, outerR, halfPhiSector, halflengthInZ, thetaMinus, thetaPlus) = ";
     tstream << "(" << m_boundValues[bv_innerRadius] << ", " << m_boundValues[bv_outerRadius] << ", " << m_boundValues[bv_halfPhiSector] << ", " << m_boundValues[bv_halfZ] << m_boundValues[bv_thetaMinus]<< ", " << m_boundValues[bv_thetaPlus] << ")";
     return tstream;
 }

}

#endif // ACTS_VOLUMES_BEVELLEDCYLINDERVOLUMESBOUNDS_H



