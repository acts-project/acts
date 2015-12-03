///////////////////////////////////////////////////////////////////
// GeoPrimitvesToStringConverter.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef GEOPRIMITIVESTOSTRINGCONVERTER_H_
#define GEOPRIMITIVESTOSTRINGCONVERTER_H_

#include "EventPrimitives/EventPrimitivesToStringConverter.h"
#include "GeoPrimitives/GeoPrimitives.h"

#ifndef ATS_GAUDI_BUILD
#   include "GeoPrimitives/CLHEPtoEigenConverter.h"
#   ifndef XAOD_STANDALONE
#       include "CLHEP/Geometry/Transform3D.h"
#       include "CLHEP/Geometry/Point3D.h"
#       include "CLHEP/Vector/TwoVector.h"
#   endif // not XAOD_STANDALONE
#endif // ATS_GAUDI_BUILD
namespace Amg {

  /** GeoPrimitvesToStringConverter

      static methods for conversion of GeoPrimitives and will call the EventPrimitives converter (Matrix)
      to std::string.

      This is to enhance formatted screen ouput and for ASCII based
      testing.

      The offset can be used to offset the lines (starting from line 2) wrt to the
      zero position for formatting reasons.

      @author Niels.Van.Eldik@cern.ch 

  */

  inline std::string toString( const Amg::Translation3D& translation, int precision = 4 ){
    Amg::Vector3D trans;
    trans[0] = translation.x();
    trans[1] = translation.y();
    trans[2] = translation.z();
    return toString( trans, precision );
  }


  inline std::string toString( const Amg::Transform3D& transform, int precision = 4, std::string offset="" ){
    std::ostringstream sout;
    sout << "Translation : " << toString( transform.translation(), precision ) << std::endl;
    std::string rotationOffset = offset + "              ";
    sout << offset << "Rotation    : " << toString( transform.rotation(), precision+2, rotationOffset );
    return sout.str();
  }

#ifndef ATS_GAUDI_BUILD
#   ifndef XAOD_STANDALONE

  inline std::string toString( const CLHEP::HepRotation& rot, int precision = 4, std::string offset="" ){
    std::ostringstream sout;

    sout << std::setiosflags(std::ios::fixed) << std::setprecision(precision);
    for( int i=0;i<3;++i ){
      for( int j=0;j<3;++j ){
	if( j == 0 ) sout << "(";
	double val = roundWithPrecision(rot(i,j),precision);
	sout << val;
	if( j == 2 ) sout << ")";
	else         sout << ", ";
      }
      if( i != 2 ) {
	sout << std::endl;
	sout << offset;    
      }       
    }    
    return sout.str();
  }


  inline std::string toString( const CLHEP::Hep3Vector& translation, int precision = 4){
    std::ostringstream sout;
    sout << std::setiosflags(std::ios::fixed) << std::setprecision(precision);
    for( int j=0;j<3;++j ){
      if( j == 0 ) sout << "(";
      double val = roundWithPrecision( translation[j],precision);
      sout << val;
      if( j == 2 ) sout << ")";
      else         sout << ", ";
    }
    return sout.str();
  }

  inline std::string toString( const CLHEP::Hep2Vector& translation, int precision = 4){
    std::ostringstream sout;
    sout << std::setiosflags(std::ios::fixed) << std::setprecision(precision);
    for( int j=0;j<2;++j ){
      if( j == 0 ) sout << "(";
      double val = roundWithPrecision( translation[j],precision);
      sout << val;
      if( j == 1 ) sout << ")";
      else         sout << ", ";
    }
    return sout.str();
  }

  inline std::string toString( const HepGeom::Transform3D& transf, int precision = 4, std::string offset=""){
    std::ostringstream sout;
    sout << "Translation : " << toString( transf.getTranslation(), precision ) << std::endl;
    std::string rotationOffset = offset + "              ";
    sout << offset << "Rotation    : " << toString( transf.getRotation(), precision+2, rotationOffset );
    return sout.str();
  }

#   endif // not XAOD_STANDALONE
#endif //ATS_GAUDI_BUILD

}

#endif /* GEOPRIMITIVESTOSTRINGCONVERTER_H_ */
