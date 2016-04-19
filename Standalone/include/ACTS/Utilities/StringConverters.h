///////////////////////////////////////////////////////////////////
// EventPrimitivesToStringConverter.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef EVENTPRIMITIVESTOSTRINGCONVERTER_H_
#define EVENTPRIMITIVESTOSTRINGCONVERTER_H_

#include "Core/AlgebraDefinitions.h"
#include <iostream>
#include <iomanip>
#include <string>

namespace Acts {

  /** EventPrimitvesToStringConverter

      inline methods for conversion of EventPrimitives (Matrix)
      to std::string.

      This is to enhance formatted screen ouput and for ASCII based
      testing.

      The offset can be used to offset the lines (starting from line 2) wrt to the
      zero position for formatting reasons.

      @author Niels.Van.Eldik@cern.ch, Andreas.Salzburger@cern.ch 

      */




     inline double roundWithPrecision( double val, int precision ) {
       if( val < 0 && fabs(val)*std::pow(10,precision) < 1. ) return -val;
       return val;
     }

     inline std::string toString( const ActsMatrixXd& matrix, int precision = 4, std::string offset="" ){
         std::ostringstream sout;

         sout << std::setiosflags(std::ios::fixed) << std::setprecision(precision);
         if( matrix.cols() == 1 ){
             sout << "(";
             for( int i=0;i<matrix.rows();++i ){
   	    double val = roundWithPrecision(matrix(i,0),precision);
   	    sout << val;
   	    if( i != matrix.rows() - 1 ) sout << ", ";
             }
             sout << ")";
         }else{
             for( int i=0;i<matrix.rows();++i ){
                 for( int j=0;j<matrix.cols();++j ){
                     if( j == 0 ) sout << "(";
   		          double val = roundWithPrecision(matrix(i,j),precision);
                     sout << val;
                     if( j == matrix.cols() - 1 ) sout << ")";
                     else                         sout << ", ";
                 }
                 if( i != matrix.rows() - 1 ) 
                 {   // make the end line and the offset in the next line
                     sout << std::endl;
                     sout << offset;
                 }                  
             }
         }
         return sout.str();
     }

     
    inline std::string toString( const Acts::Translation3D& translation, int precision = 4 ){
      Acts::Vector3D trans;
      trans[0] = translation.x();
      trans[1] = translation.y();
      trans[2] = translation.z();
      return toString( trans, precision );
    }
     
    inline std::string toString( const Acts::Transform3D& transform, int precision = 4, std::string offset="" ){
      std::ostringstream sout;
      sout << "Translation : " << toString( transform.translation(), precision ) << std::endl;
      std::string rotationOffset = offset + "              ";
      sout << offset << "Rotation    : " << toString( transform.rotation(), precision+2, rotationOffset );
      return sout.str();
    }
   

}

#endif /* EVENTPRIMITIVESTOSTRINGCONVERTER_H_ */
