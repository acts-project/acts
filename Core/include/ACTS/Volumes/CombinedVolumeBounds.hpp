///////////////////////////////////////////////////////////////////
// CombinedVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_COMBINEDVOLUMEBOUNDS_H
#define ACTS_VOLUMES_COMBINEDVOLUMEBOUNDS_H 1

// Geometry module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"
#include "ACTS/Volumes/Volume.hpp"
// Core module

namespace Acts {

   class SurfaceBounds;
   class Volume;
   class Surface;

  /**
   @class CombinedVolumeBounds

   Bounds for a generic combined volume, the decomposeToSurfaces method creates a
   vector of n surfaces (n1+n2-nshared):

   BoundarySurfaceFace [index]: [n1+n2-nshared] combined surfaces

   designed to allow transcript of GeoShapeUnion and GeoShapeIntersection

   @author Sarka.Todorova@cern.ch
  */

 class CombinedVolumeBounds : public VolumeBounds {

  public:
    /**Default Constructor*/
    CombinedVolumeBounds();

    /**Constructor - the box boundaries */
    CombinedVolumeBounds( Volume* first,Volume* second, bool intersection);

    /**Copy Constructor */
    CombinedVolumeBounds(const CombinedVolumeBounds& bobo);

    /**Destructor */
    virtual ~CombinedVolumeBounds();

    /**Assignment operator*/
    CombinedVolumeBounds& operator=(const CombinedVolumeBounds& bobo);

    /**Virtual constructor */
    CombinedVolumeBounds* clone() const override;

    /**This method checks if position in the 3D volume frame is inside the volume*/
    bool inside(const Vector3D&, double tol=0.) const override;

    /** Method to decompose the Bounds into boundarySurfaces */
    const std::vector<const Acts::Surface*>* decomposeToSurfaces(std::shared_ptr<Transform3D>transformPtr) const override;

    /**This method returns the first VolumeBounds*/
    Volume* first() const;

    /**This method returns the second VolumeBounds*/
    Volume* second() const;

    /**This method distinguishes between Union(0) and Intersection(1)*/
    bool intersection() const;

    /**This method returns bounds orientation*/
    const std::vector<bool> boundsOrientation() const;

    /** Output Method for std::ostream */
    std::ostream& dump(std::ostream& sl) const override;

  private:

    Acts::Volume* createSubtractedVolume(const Transform3D& transf, Acts::Volume* subtrVol) const;

    Volume*                     m_first;                   //!< first volume of the combination
    Volume*                     m_second;                  //!< second volume of the combination
    bool                        m_intersection;            //!< boolean if intersection is needed or not
    mutable std::vector<bool>   m_boundsOrientation;       //!< orientation of the bounds

 };

 inline CombinedVolumeBounds* CombinedVolumeBounds::clone() const
 { return new CombinedVolumeBounds(*this); }

 inline bool CombinedVolumeBounds::inside(const Vector3D &pos, double tol) const
 {
   if (m_intersection) return (m_first->inside(pos,tol) && m_second->inside(pos,tol) );
   return (m_first->inside(pos,tol) || m_second->inside(pos,tol) );
 }

 inline Volume* CombinedVolumeBounds::first() const { return m_first; }

 inline Volume* CombinedVolumeBounds::second() const { return m_second; }

 inline bool CombinedVolumeBounds::intersection() const { return m_intersection; }

 inline const std::vector<bool> CombinedVolumeBounds::boundsOrientation() const
  { return(m_boundsOrientation); }

}

#endif // ACTS_VOLUMES_COMBINEDVOLUMEBOUNDS_H
