///////////////////////////////////////////////////////////////////
// GeometryObject.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYOBJECT_H
#define ACTS_GEOMETRYUTILS_GEOMETRYOBJECT_H

// Geometry module
#include "ACTS/GeometryUtils/GeometryID.h"
#include "ACTS/GeometryUtils/BinningType.h"
// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"


namespace Acts {


  /** @class GeometryObject 

      Base class to provide GeometryID interface:
      - simple set and get
    
      It also provides the binningPosition method for 
      Geometry geometrical object to be binned in BinnedArrays

      @author Andreas.Salzburger@cern.ch 
    */

  class GeometryObject {
    public:
      /** constructor from a ready-made value */    
      GeometryObject() :
        m_geoID()
      {}

      /** return the value */     
      const GeometryID& geoID() const;

      /** set the value */
      void assignGeoID(const GeometryID& geoID) const; 
      
      /** force a binning position method */
      virtual Vector3D binningPosition(BinningValue bValue) const = 0;
    
      /** implement the binningValue */
      double binningPositionValue(BinningValue bValue) const;
          
    protected:      
      mutable GeometryID m_geoID;
  };

  inline const GeometryID& GeometryObject::geoID() const  { return m_geoID; }

  inline void GeometryObject::assignGeoID(const GeometryID& geoID) const  { m_geoID = geoID; }
  
  inline double GeometryObject::binningPositionValue(BinningValue bValue) const {
      // now switch     
      switch (bValue) {
          // case x
          case Acts::binX : { 
              return binningPosition(bValue).x();
          } break;
          // case y
          case Acts::binY : { 
              return binningPosition(bValue).y();
          } break;
          // case z
          case Acts::binZ : { 
              return binningPosition(bValue).z();
          } break;
          // case 
          case Acts::binR : {
              return binningPosition(bValue).perp();
          } break;
          // do nothing for the default
          default : return 0;
       }          
  }


}

#endif
