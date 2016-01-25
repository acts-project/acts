///////////////////////////////////////////////////////////////////
// EventPrimitives.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////


#ifndef EVENT_EVENTPRIMITIVES_H
#define EVENT_EVENTPRIMITIVES_H

#define EIGEN_MATRIXBASE_PLUGIN "EventPrimitives/AmgMatrixPlugin.h"
#define EIGEN_MATRIX_PLUGIN "EventPrimitives/SymmetricMatrixHelpers.h"
#define EIGEN_TRANSFORM_PLUGIN "EventPrimitives/AmgTransformPlugin.h"

#include <unistd.h>
#include <Eigen/Core>
#include <Eigen/Dense>

// These are the typedefs from Eigen to AMG ( Atlas Math and Geometry )
// some of those can be refined when switching to C++11
//  e.g. typedef Eigen::Matrix<scalar, int r, 1 0, int mr, 1> Vector<scalar, int r, int mr>;
// @author: Eigen
// @responsible: Andreas.Salzburger -at- cern.ch

namespace Amg {

    /** Dynamic Matrix - not recommended */
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>       MatrixX;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>       SymMatrixX;
                                                                        
    /** Dynamic Vector - not recommended */                             
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1>                    VectorX;
     
    /** Macros for fixed size - very recommended */ 
    #ifndef AmgMatrixDef
    #define AmgMatrixDef
    #define AmgMatrix(rows, cols) Eigen::Matrix<double,rows,cols,0,rows,cols>
    #define AmgSymMatrix(dim) Eigen::Matrix<double,dim,dim,0,dim,dim>
    #endif
                    
    #ifndef AmgVectorDef
    #define AmgVectorDef
    #define AmgVector(rows) Eigen::Matrix<double, rows, 1, 0, rows, 1>
    #define AmgRowVector(cols) Eigen::Matrix<double, 1, cols, 0, 1, cols>
    #endif          
    

}

#endif /* EVENT_EVENTPRIMITIVES_H */
