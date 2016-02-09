#ifndef ATS_ALGEBRADEFINITIONS_H
#define ATS_ALGEBRADEFINITIONS_H

// API addons 
#define EIGEN_MATRIXBASE_PLUGIN "Algebra/MatrixBasePlugin.h"
#define EIGEN_MATRIX_PLUGIN "Algebra/MatrixPlugin.h"
#define EIGEN_TRANSFORM_PLUGIN "Algebra/TransformPlugin.h"

// external include(s)
#include <Eigen/Dense>

namespace Ats
{
  template<typename T, unsigned int rows, unsigned int cols>
  using AtsMatrix = Eigen::Matrix<T,rows,cols>;

  template<unsigned int rows, unsigned int cols>
  using AtsMatrixD = AtsMatrix<double,rows,cols>;

  template<unsigned int rows, unsigned int cols>
  using AtsMatrixF = AtsMatrix<float,rows,cols>;

  template<typename T, unsigned int rows>
  using AtsSymMatrix = Eigen::Matrix<T,rows,rows>;

  template<unsigned int rows>
  using AtsSymMatrixD = AtsSymMatrix<double,rows>;

  template<unsigned int rows>
  using AtsSymMatrixF = AtsSymMatrix<float,rows>;

  template<typename T, unsigned int rows>
  using AtsVector = Eigen::Matrix<T,rows,1>;

  template<unsigned int rows>
  using AtsVectorD = AtsVector<double,rows>;

  template<unsigned int rows>
  using AtsVectorF = AtsVector<float,rows>;

  template<typename T, unsigned int cols>
  using AtsRowVector = Eigen::Matrix<T,1,cols>;

  template<unsigned int cols>
  using AtsRowVectorD = AtsRowVector<double,cols>;

  template<unsigned int cols>
  using AtsRowVectorF = AtsRowVector<float,cols>;

  template<typename T>
  using AtsMatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  using AtsMatrixXd = AtsMatrixX<double>;
  using AtsMatrixXf = AtsMatrixX<float>;

  template<typename T>
  using AtsVectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  using AtsVectorXd = AtsVectorX<double>;
  using AtsVectorXf = AtsVectorX<float>;

  template<typename T>
  using AtsRowVectorX = Eigen::Matrix<T, 1, Eigen::Dynamic>;

  using AtsRowVectorXd = AtsRowVectorX<double>;
  using AtsRowVectorXf = AtsRowVectorX<float>;
  
 
  /** elment for code readability
      - please use these for access to the member variables if needed, e.g.
          double z  = position[Ats::z];
          double px = momentum[Ats::px];
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
  

}  // end of namespace Ats

#endif
