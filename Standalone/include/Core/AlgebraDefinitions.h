#ifndef ACTS_ALGEBRADEFINITIONS_H
#define ACTS_ALGEBRADEFINITIONS_H

// API addons
#define EIGEN_MATRIXBASE_PLUGIN "Core/Algebra/MatrixBasePlugin.h"
#define EIGEN_MATRIX_PLUGIN "Core/Algebra/MatrixPlugin.h"
#define EIGEN_TRANSFORM_PLUGIN "Core/Algebra/TransformPlugin.h"

// external include(s)
#include <Eigen/Dense>

namespace Acts
{
  template<typename T, unsigned int rows, unsigned int cols>
  using ActsMatrix = Eigen::Matrix<T,rows,cols>;

  template<unsigned int rows, unsigned int cols>
  using ActsMatrixD = ActsMatrix<double,rows,cols>;

  template<unsigned int rows, unsigned int cols>
  using ActsMatrixF = ActsMatrix<float,rows,cols>;

  template<typename T, unsigned int rows>
  using ActsSymMatrix = Eigen::Matrix<T,rows,rows>;

  template<unsigned int rows>
  using ActsSymMatrixD = ActsSymMatrix<double,rows>;

  template<unsigned int rows>
  using ActsSymMatrixF = ActsSymMatrix<float,rows>;

  template<typename T, unsigned int rows>
  using ActsVector = Eigen::Matrix<T,rows,1>;

  template<unsigned int rows>
  using ActsVectorD = ActsVector<double,rows>;

  template<unsigned int rows>
  using ActsVectorF = ActsVector<float,rows>;

  template<typename T, unsigned int cols>
  using ActsRowVector = Eigen::Matrix<T,1,cols>;

  template<unsigned int cols>
  using ActsRowVectorD = ActsRowVector<double,cols>;

  template<unsigned int cols>
  using ActsRowVectorF = ActsRowVector<float,cols>;

  template<typename T>
  using ActsMatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  using ActsMatrixXd = ActsMatrixX<double>;
  using ActsMatrixXf = ActsMatrixX<float>;

  template<typename T>
  using ActsVectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  using ActsVectorXd = ActsVectorX<double>;
  using ActsVectorXf = ActsVectorX<float>;

  template<typename T>
  using ActsRowVectorX = Eigen::Matrix<T, 1, Eigen::Dynamic>;

  using ActsRowVectorXd = ActsRowVectorX<double>;
  using ActsRowVectorXf = ActsRowVectorX<float>;


  /** elment for code readability
      - please use these for access to the member variables if needed, e.g.
          double z  = position[Acts::eZ];
          double px = momentum[Acts::ePX];
  */
  enum AxisDefs {
      // position access
      eX = 0,
      eY = 1,
      eZ = 2,
      // momentum access
      ePX = 0,
      ePY = 1,
      ePZ = 2
  };

  typedef Eigen::Quaternion<double>                   Rotation3D;
  typedef Eigen::Translation<double, 3>               Translation3D;
  typedef Eigen::AngleAxisd                           AngleAxis3D;
  typedef Eigen::Affine3d                             Transform3D;
  typedef Eigen::Matrix<double, 3, 1>                 Vector3D;
  typedef Eigen::Matrix<double, 2, 1>                 Vector2D;
  typedef Eigen::Matrix<double, 3, 3>                 RotationMatrix3D;


}  // end of namespace Acts

#endif
