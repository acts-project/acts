///////////////////////////////////////////////////////////////////
// GeoPrimitives.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef GEOPRIMITIVES_GEOPRIMITIVES_H
#define GEOPRIMITIVES_GEOPRIMITIVES_H

#define EIGEN_MATRIXBASE_PLUGIN "EventPrimitives/AmgMatrixPlugin.h"
#define EIGEN_MATRIX_PLUGIN "EventPrimitives/SymmetricMatrixHelpers.h"
#define EIGEN_TRANSFORM_PLUGIN "EventPrimitives/AmgTransformPlugin.h"

#include <unistd.h>
#include <Eigen/Geometry>

/** Definition of ATLAS Math & Geometry primitives (Amg) 

    This is based on the Eigen geometry module:
    http://eigen.tuxfamily.org/dox/group__Geometry__Module.html

    @author  Robert Johannes Langenberg <robert.langenberg@cern.ch>
    @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
        
*/

namespace Amg {

    /** elment for code readability
        - please use these for access to the member variables if needed, e.g.
            double z  = position[Amg::z];
            double px = momentum[Amg::px];
    */
    enum AxisDefs {
        // position access
        x = 0,
        y = 1,
        z = 2,
        // momentum access
        px = 0,
        py = 1,
        pz = 2
    };

    typedef Eigen::Quaternion<double>                   Rotation3D;
    typedef Eigen::Translation<double, 3>               Translation3D;
    typedef Eigen::AngleAxisd                           AngleAxis3D;
    typedef Eigen::Affine3d                             Transform3D;
    typedef Eigen::Matrix<double, 3, 1>                 Vector3D;
    typedef Eigen::Matrix<double, 2, 1>                 Vector2D;
    typedef Eigen::Matrix<double, 3, 3>                 RotationMatrix3D;
    
   


}
#endif /* GEOPRIMITIVES_GEOPRIMITIVES_H */
