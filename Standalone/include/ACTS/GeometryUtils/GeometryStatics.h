///////////////////////////////////////////////////////////////////
// GeometryStatics.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H
#define ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H 1

#include "Core/AlgebraDefinitions.h"

/** Define statics for Geometry in Tracking
 */
namespace Acts {

// transformations

  static Transform3D s_idTransform = Transform3D::Identity();          //!< idendity transformation
  static Rotation3D s_idRotation = Acts::Rotation3D::Identity();            //!< idendity rotation
  
  // axis system
  
  static Vector3D s_xAxis(1,0,0);        //!< global x Axis;
  static Vector3D s_yAxis(0,1,0);        //!< global y Axis;
  static Vector3D s_zAxis(0,0,1);        //!< global z Axis;
  
  // origin
  
  static Vector3D s_origin(0,0,0);       //!< origin position
  
  static double helper[9] = {0.,1.,0.,1.,0.,0.,0.,0.,-1.};
  
  static Acts::RotationMatrix3D s_idRotationZinverse(helper);

}

#endif // ACTS_GEOMETYRUTILS_GEOMETRYSTATICS_H
