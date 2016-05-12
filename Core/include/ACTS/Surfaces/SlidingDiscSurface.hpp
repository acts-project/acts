///////////////////////////////////////////////////////////////////
// SlidingDiscSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SLIDINGDISCSURFACE_H
#define ACTS_SURFACES_SLIDINGDISCSURFACE_H 1

// Geometry module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Surfaces/NoBounds.hpp"
#include "ACTS/Utilities/BinUtility.hpp"

namespace Acts {

  class DiscBounds;
  class DetectorElementBase;

  template<int DIM, class T, class S> class ParametersT;

  /**
   @class SlidingDiscSurface

   Class for a Calo DiscSurface with variable depth in the ATLAS detector.
   The variable depth is stored as a binned vector of shifts of the center along the normal.
   Local eta bin is defined by Acts::eLOC_R and z position for base transform ( corrected for misalignement ).
   Inherits from DiscSurface.

   @author sarka.todorova@cern.ch
   */

  class SlidingDiscSurface : public DiscSurface {

    public:
      /** Default Constructor*/
      SlidingDiscSurface();

      /** Constructor */
      SlidingDiscSurface(const DiscSurface& surf,
                         BinUtility* bu = nullptr,
                         const std::vector<float>* offset = nullptr,
                         Transform3D* align = nullptr);

      /** Copy Constructor*/
      SlidingDiscSurface(const SlidingDiscSurface& psf);

      /** Copy Constructor with shift*/
      SlidingDiscSurface(const SlidingDiscSurface& psf, const Transform3D& transf);

      /** Destructor*/
      virtual ~SlidingDiscSurface();

      /** Assignement operator*/
      SlidingDiscSurface& operator=(const SlidingDiscSurface&dsf);

      /** Equality operator*/
      bool operator==(const Surface& sf) const;

      /** Virtual constructor - shift can be optionally applied */
      virtual SlidingDiscSurface* clone(const Transform3D* shift = nullptr) const;

      /** This method returns true if the GlobalPosition is on the Surface for both, within
        or without check of whether the local position is inside boundaries or not */
      bool isOnSurface(const Vector3D& glopo, const BoundaryCheck& bchk=true) const;

      /** Specialized for DiscSurface: LocalToGlobal method without dynamic memory allocation */
      void localToGlobal(const Vector2D& locp, const Vector3D& mom, Vector3D& glob) const;

      /** Specialized for DiscSurface: GlobalToLocal method without dynamic memory allocation - boolean checks if on surface */
      bool globalToLocal(const Vector3D& glob, const Vector3D& mom, Vector2D& loc) const;

      /**This method allows access to the bin utility*/
      const BinUtility* binUtility() const { return m_etaBin; }

      /**This method allows access to the radial offset values*/
      const std::vector<float>* offset() const { return m_depth; }

      /** Return properly formatted class name for screen output */
      std::string name() const { return "Acts::SlidingDiscSurface"; }

    protected: //!< data members
      const std::vector<float>*                         m_depth;
      BinUtility*                                  m_etaBin;
      Transform3D*                                      m_align;
 };

  inline SlidingDiscSurface* SlidingDiscSurface::clone(const Transform3D* shift) const
  {
    if (shift) new SlidingDiscSurface(*this,*shift);
    return new SlidingDiscSurface(*this);
  }

} // end of namespace

#endif // ACTS_SURFACES_SLIDINGDISCSURFACE_H


