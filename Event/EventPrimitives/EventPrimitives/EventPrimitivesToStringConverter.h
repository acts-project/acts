///////////////////////////////////////////////////////////////////
// EventPrimitivesToStringConverter.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef EVENTPRIMITIVESTOSTRINGCONVERTER_H_
#define EVENTPRIMITIVESTOSTRINGCONVERTER_H_

#include "EventPrimitives/EventPrimitives.h"
#ifndef ATS_GAUDI_BUILD
#   ifndef XAOD_STANDALONE
#       include "CLHEP/Matrix/GenMatrix.h"
#   endif // not XAOD_STANDALONE
#endif // ATS_GAUDI_BUILD
#include <iostream>
#include <iomanip>
#include <string>

namespace Amg {

  /** EventPrimitvesToStringConverter

      inline methods for conversion of EventPrimitives (Matrix)
      to std::string.

      This is to enhance formatted screen ouput and for ASCII based
      testing.

      The offset can be used to offset the lines (starting from line 2) wrt to the
      zero position for formatting reasons.

      @author Niels.Van.Eldik@cern.ch 

      */

  inline double roundWithPrecision( double val, int precision ) {
    if( val < 0 && fabs(val)*std::pow(10,precision) < 1. ) return -val;
    return val;
  }


  inline std::string toString( const MatrixX& matrix, int precision = 4, std::string offset="" ){
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
#ifndef ATS_GAUDI_BUILD
#   ifndef XAOD_STANDALONE
  inline std::string toString( const CLHEP::HepGenMatrix& matrix, int precision = 4, std::string offset="" ){
    std::ostringstream sout;

    sout << std::setiosflags(std::ios::fixed) << std::setprecision(precision);
    for( int i=0;i<matrix.num_row();++i ){
      for( int j=0;j<matrix.num_col();++j ){
	if( j == 0 ) sout << "(";
	double val = roundWithPrecision(matrix(i+1,j+1),precision);
	sout << val;
	if( j == matrix.num_col() - 1 ) sout << ")";
	else                         sout << ", ";
      }
      if( i != matrix.num_row() - 1 ) 
	{   // make the end line and the offset in the next line
	  sout << std::endl;
	  sout << offset;
	}                  
    }
    return sout.str();
  }
#   endif // not XAOD_STANDALONE
#endif //ATS_GAUDI_BUILD
}

#endif /* EVENTPRIMITIVESTOSTRINGCONVERTER_H_ */
