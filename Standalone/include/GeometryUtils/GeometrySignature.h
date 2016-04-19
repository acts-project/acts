///////////////////////////////////////////////////////////////////
// GeometrySignature.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYSIGNATURE_H 
#define ACTS_GEOMETRYUTILS_GEOMETRYSIGNATURE_H 1

namespace Acts {

  /** @class GeometrySignature
  
     An enumeration object that puts the signature
     of a GeometryBuilder to all subvolumes
    
    @TODO will be in the future be replace by GeometryID mechanism
      
    @author Andreas.Salzburger@cern.ch
    
    */
    
    enum GeometrySignature {
        Global                 =  0,       
        ID                     =  1,
        BeamPipe               =  2,
        Calo                   =  3,
        MS                     =  4,
        Cavern                 =  5,
        NumberOfSignatures     =  6,
        Unsigned               = 99
    };

    enum GeometryType {
        Static                  = 0,
        Dense                   = 1,
        DenseWithLayers         = 1,
        Detached                = 2,
        Master                  = 3,
        NumberOfGeometryTypes   = 3   
    };


} // end of namespace

#endif // ACTS_GEOMETRYUTILS_GEOMETRYSIGNATURE_H
