///////////////////////////////////////////////////////////////////
// GeometryID.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYID_H
#define ACTS_GEOMETRYUTILS_GEOMETRYID_H 1

// STD/STL
#include <iostream>

namespace Acts {

  typedef uint64_t geo_id_value;

  /** @class GeometryID

      Identifier for Geometry nodes - packing the
      - (Sensitive) Surfaces    - uses IdentiferHash
      - (Approach)  Surfaces    - uses simple counting through confinedLayers() or confinedArbitraryLayers()
      - (Layer)     Surfaces    - uses simple counting through approachSurfaces()
      - (Boundary)  Surfaces    - uses simple counting through boundarySurfaces()
      -  Volumes                - AtlasDetector region and AtlasSubDetector enum and static counter

      @author Andreas.Salzburger@cern.ch      */

  class GeometryID {
    public:
      /** constructor from a ready-made value */
      GeometryID(geo_id_value id_value = 0) :
        m_value(id_value)
      {}

      /** Copy constructor */
      GeometryID(const GeometryID& tddID) :
        m_value(tddID.m_value)
      {}

      /** Assignement operator */
      GeometryID& operator=(const GeometryID& tddID)
      {
         if (&tddID != this){
             m_value = tddID.m_value;
         }
         return (*this);
      }

      /** return the value */
      geo_id_value value() const;

    private:
      geo_id_value m_value;

  };

  inline geo_id_value GeometryID::value() const { return m_value; }

  /** Overload of operator< | <= | > | >=  for the usage in a map */
  bool operator< ( const GeometryID& one, const GeometryID& two );
  bool operator<=( const GeometryID& one, const GeometryID& two );
  bool operator> ( const GeometryID& one, const GeometryID& two );
  bool operator>=( const GeometryID& one, const GeometryID& two );

  /**Overload of << operator for std::ostream for debug output*/
  std::ostream& operator<<( std::ostream& sl, const GeometryID& tddID);

}

#endif // ACTS_GEOMETRYUTILS_GEOMETRYID_H
