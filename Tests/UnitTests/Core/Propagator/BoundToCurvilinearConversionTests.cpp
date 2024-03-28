// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log_formatter.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <typeinfo>

using namespace Acts;
using namespace Acts::UnitLiterals;

bool isDebugOutputEnabled() {
  std::array<boost::unit_test::output_format, 1> formats{
      boost::unit_test::OF_CLF};
  for (auto a_format : formats) {
    auto formatter = ::boost::unit_test::unit_test_log.get_formatter(a_format);
    if (formatter != nullptr) {
      return formatter->get_log_level() < boost::unit_test::log_test_units;
    }
  }
  return false;
}

#define MSG_DEBUG(a)               \
  if (isDebugOutputEnabled()) {    \
    std::stringstream msg;         \
    msg << a;                      \
    BOOST_TEST_MESSAGE(msg.str()); \
  }                                \
  do {                             \
  } while (0)

namespace {
/** Helper function to compute dt/ds
 * Helper function to compute the derivative of the time as function of the path
 * length
 */
template <class ParticleHypothesis>
double computeDtDs(const ParticleHypothesis &hypothesis, double qop) {
  return std::hypot(1., hypothesis.mass() / hypothesis.extractMomentum(qop));
}

/** Compute the path length derivatives for the free/bound to curvilinear
 * parameter transform.
 * @param direction the direction of trajectory at the location in question
 * @param qop q/p of the particle at the location in question in Acts units
 * @param bfield the magnetic field at the location in question in Acts units
 * @param particle_hypothesis the particle hypothesis e.g. Acts::ParticleHypothesis::pion()
 * @return path length derivatives ( dr(...)/ds, dt/ds, dr(...)/ds2, d qop/ds [== 0] )
 */
template <class ParticleHypothesis>
FreeToPathMatrix computeFreeToPathDerivatives(
    const Vector3 &direction, double qop, const Vector3 &bField,
    const ParticleHypothesis &particleHypothis) {
  FreeToPathMatrix pathLengthDeriv;
#if defined(EIGEN_HAS_CONSTEXPR) && EIGEN_VERSION_AT_LEAST(3, 4, 0)
  static_assert(pathLengthDeriv.cols() ==
                8);  // ensure that all elements are initialized
#endif
  pathLengthDeriv.segment<3>(eFreePos0) = direction;
  pathLengthDeriv(0, eFreeTime) = computeDtDs(particleHypothis, qop);
  pathLengthDeriv.segment<3>(eFreeDir0) =
      (qop * direction.cross(bField)).transpose();
  pathLengthDeriv(0, Acts::eFreeQOverP) = 0.;
  return pathLengthDeriv;
}
template <typename T, std::size_t Rows, std::size_t Cols>
inline constexpr Eigen::Matrix<T, Rows, Cols> makeMatrix(
    std::initializer_list<double> elements) {
  // static_assert( elements.size() == Rows*Cols )
  if (!(elements.size() == Rows * Cols)) {
    // throw std::range_error("Initializer list size does not match matrix
    // dimensions.");
    std::abort();
  }
  Eigen::Matrix<T, Rows, Cols> matrix;
  auto iter = elements.begin();
  for (unsigned int row = 0; row < matrix.rows(); ++row) {
    for (unsigned int col = 0; col < matrix.cols(); ++col) {
      matrix(row, col) = *iter;
      ++iter;
    }
  }
  return matrix;
}
template <typename T, std::size_t Rows>
inline constexpr Eigen::Matrix<T, Rows, 1> makeVector(
    std::initializer_list<double> elements) {
  return makeMatrix<T, Rows, 1>(elements);
}

template <typename Matrix>
Matrix matrixRatio(const Matrix &a, const Matrix &b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    std::abort();
  }
  Matrix ret;
  for (unsigned int row = 0; row < a.rows(); ++row) {
    for (unsigned int col = 0; col < a.cols(); ++col) {
      if (b(row, col) == 0.) {
        ret(row, col) = a(row, col) - b(row, col);
      } else {
        ret(row, col) = a(row, col) / b(row, col);
      }
    }
  }
  return ret;
}

}  // namespace

struct TestData {
  enum ESurfaceType { kPlane, kPolarDisk, kCylinder };
  TestData(Vector3 &&a_surfaceCenter, ActsMatrix<3, 3> &&a_surfaceRot,
           ESurfaceType a_surfaceType, BoundVector &&a_paramVec,
           BoundSquareMatrix &&a_paramCov, Vector3 &&a_bField)
      : surfaceCenter(std::move(a_surfaceCenter)),
        surfaceRot(std::move(a_surfaceRot)),
        surfaceType(a_surfaceType),
        paramVec(std::move(a_paramVec)),
        paramCov(std::move(a_paramCov)),
        bfield(std::move(a_bField)) {}

  Vector3 surfaceCenter;
  ActsMatrix<3, 3> surfaceRot;
  ESurfaceType surfaceType;
  BoundVector paramVec;
  BoundSquareMatrix paramCov;
  Vector3 bfield;
};

template <typename StepperCreator>
void test_bound_to_curvilinear(const std::vector<TestData> &testDataList,
                               const StepperCreator &stepperCreator) {
  GeometryContext geoCtx;
  MagneticFieldContext magFieldContext;

  for (const auto &testData : testDataList) {
    // create a constant magnetic field provider for the test_data
    std::shared_ptr<Acts::MagneticFieldProvider> bField =
        std::dynamic_pointer_cast<Acts::MagneticFieldProvider>(
            std::make_shared<ConstantBField>(testData.bfield));

    // create bound parameters from test data
    const Vector3 &surface_center = testData.surfaceCenter;
    const ActsMatrix<3, 3> &surfaceRot = testData.surfaceRot;
    const BoundVector &paramVec = testData.paramVec;
    const BoundSquareMatrix &cov = testData.paramCov;

    AngleAxis3 surfaceTransform0;
    surfaceTransform0 = surfaceRot;

    std::shared_ptr<Surface> surface;
    switch (testData.surfaceType) {
      case TestData::kPlane: {
        surface = std::dynamic_pointer_cast<Surface>(
            Surface::makeShared<PlaneSurface>(Translation3(surface_center) *
                                              surfaceTransform0));
        break;
      }
      case TestData::kPolarDisk: {
        surface =
            std::dynamic_pointer_cast<Surface>(Surface::makeShared<DiscSurface>(
                Translation3(surface_center) * surfaceTransform0));
        break;
      }
      default: {
        throw std::runtime_error("Unhandled surface type.");
        std::abort();
      }
    }

    FreeVector freeParamVec = Acts::detail::transformBoundToFreeParameters(
        *surface, geoCtx, paramVec);

    Vector3 direction = freeParamVec.segment<3>(eFreeDir0);
    Vector3 position = surface->localToGlobal(
        geoCtx, Vector2{paramVec[0], paramVec[1]}, direction);
    BoundTrackParameters params(surface, paramVec,
                                std::optional<BoundSquareMatrix>(cov),
                                ParticleHypothesis::pion());

    // compute curvilinear parameters by using the propagator with
    // a small step size : ==0, >0 but below path limit,  > path limit, *10 and
    // > path limit
    for (unsigned int i = 0; i < 4; ++i) {
      MagneticFieldProvider::Cache cache = bField->makeCache(magFieldContext);

      Result<Acts::Vector3> localBField = bField->getField(position, cache);
      assert(localBField.ok());

      auto pathLengthDerivatives = computeFreeToPathDerivatives(
          direction, params.parameters()[eBoundQOverP], localBField.value(),
          ParticleHypothesis::pion());
      MSG_DEBUG("derivatives : " << pathLengthDerivatives);

      // compute Jacobian for bound to curvilinear covariance transformation
      Acts::BoundMatrix b2c;
      Acts::detail::boundToCurvilinearTransportJacobian(
          direction, surface->boundToFreeJacobian(geoCtx, freeParamVec),
          Acts::FreeMatrix::Identity(),
          computeFreeToPathDerivatives(
              direction, params.parameters()[eBoundQOverP], localBField.value(),
              ParticleHypothesis::pion()),
          b2c);

      auto curviCovAlt = b2c * cov * b2c.transpose();

      MSG_DEBUG("curvilinear covariance alt.:" << std::endl << curviCovAlt);

      // configure propagator for tiny step size
      Acts::PropagatorOptions<> nullPropagationOptions(geoCtx, magFieldContext);

      nullPropagationOptions.pathLimit =
          i == 0 ? 0 : 1e-12 * 1_m * std::pow(10, i - 1);
      if (nullPropagationOptions.pathLimit > 0 && i > 1) {
        nullPropagationOptions.stepTolerance =
            nullPropagationOptions.pathLimit * .99;
        nullPropagationOptions.surfaceTolerance =
            nullPropagationOptions.pathLimit * .99;
      }
      auto stepper = stepperCreator(bField);
      MSG_DEBUG("Stepper type " << typeid(stepper).name());

      auto logLevel = (isDebugOutputEnabled() ? Acts::Logging::VERBOSE
                                              : Acts::Logging::INFO);

      // Use propagator with small step size to compute parameters in
      // curvilinear parameterisation
      Propagator<decltype(stepper)> propagator(
          std::move(stepper), Acts::VoidNavigator(),
          Acts::getDefaultLogger("Propagator", logLevel));
      auto result = propagator.propagate(params, nullPropagationOptions, true);
      {
        const auto &curvilinearParameters = result.value().endParameters;
        if (curvilinearParameters.has_value() &&
            curvilinearParameters.value().covariance().has_value()) {
          MSG_DEBUG(i << " | "
                      << "limit: " << nullPropagationOptions.pathLimit
                      << " tolerance: " << nullPropagationOptions.stepTolerance
                      << std::endl);

          Acts::BoundSquareMatrix curviCov =
              curvilinearParameters.value().covariance().value();
          MSG_DEBUG("curvilinear covariance:" << std::endl
                                              << curviCov << std::endl);
          if (isDebugOutputEnabled()) {
            Acts::BoundSquareMatrix b(curviCovAlt);
            auto ratio = matrixRatio(curviCov, b);
            MSG_DEBUG("ratio:" << std::endl << ratio << std::endl);
          }
          // test that result from propagation and explicit computation are
          // compatible.
          BOOST_CHECK(curviCovAlt.isApprox(curviCov));
        }
      }
    }
  }
}

std::vector<TestData> make_test_data(double magFieldScale = 1.) {
  std::vector<TestData> testDataList{
      TestData(
          makeVector<ActsScalar, 3>(
              {-442.883, -624.094, 857.272}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {0.677197, 0.0176111, -0.735591, -0.735342, -0.0191232, -0.677426,
               -0.0259971, 0.999662, -2.22045e-16}),  // surface rot
          TestData::kPlane,
          makeVector<ActsScalar, eBoundSize>({46.5758, 4.5564, -2.38067,
                                              0.72974, 0.73159,
                                              1163.57}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.00025406,   0.00334274,  2.9713e-06,   -6.40317e-06,
               2.52229e-05,  0.00291208,  0.00334274,   9.77017,
               0.00109081,   -0.0106064,  -0.00340842,  7.25206,
               2.9713e-06,   0.00109081,  4.86984e-07,  -1.68459e-06,
               2.0707e-06,   0.000848693, -6.40317e-06, -0.0106064,
               -1.68459e-06, 1.58516e-05, 4.40043e-06,  -0.00786289,
               2.52229e-05,  -0.00340842, 2.0707e-06,   4.40043e-06,
               2.786e-05,    -0.00210611, 0.00291208,   7.25206,
               0.000848693,  -0.00786289, -0.00210611,  89880.9}),  // param cov
          makeVector<ActsScalar, 3>({-2.64634e-05 * 1000_T,
                                     -4.38183e-05 * 1000_T,
                                     0.00197353 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-215.895, 979.521, 808.928}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {-0.999319, -0.0259882, -0.0261772, -0.0261683, -0.00068053,
               0.999657, -0.0259971, 0.999662, -2.22045e-16}),  // surface rot
          TestData::kPlane,
          makeVector<ActsScalar, eBoundSize>({-47.2414, -20.7881, 1.46297,
                                              0.926114, 0.723167,
                                              1318.63}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.000299382,  -0.00331811,  -9.93116e-06,
               4.79934e-06,  9.50183e-06,  -0.00184948,
               -0.00331811,  0.212531,     -3.5517e-05,
               -0.00030374,  3.77471e-05,  0.129021,
               -9.93116e-06, -3.5517e-05,  1.26087e-06,
               2.63359e-08,  -1.11054e-06, -4.07474e-05,
               4.79934e-06,  -0.00030374,  2.63359e-08,
               2.42802e-06,  1.77196e-07,  -0.000180912,
               9.50183e-06,  3.77471e-05,  -1.11054e-06,
               1.77196e-07,  2.13352e-05,  0.000394159,
               -0.00184948,  0.129021,     -4.07474e-05,
               -0.000180912, 0.000394159,  89875.6}),  // param cov
          makeVector<ActsScalar, 3>({-7.28154e-06 * 1000_T,
                                     4.91679e-05 * 1000_T,
                                     0.00200021 * 1000_T})),  // magnetic field

      TestData(makeVector<ActsScalar, 3>(
                   {-100.1, 9.9476e-14, 2623}),  // surface center
               makeMatrix<ActsScalar, 3, 3>({-9.4369e-16, -1, -2.22045e-16, -1,
                                             9.4369e-16, 2.09541e-31, 0,
                                             2.22045e-16, -1}),  // surface rot
               TestData::kPlane,
               makeVector<ActsScalar, eBoundSize>({5.1223, 16.6267, -3.08166,
                                                   0.0439704, -0.0358564,
                                                   2625.58}),  // param_vec
               makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
                   {0.00017846,   6.32301e-09,  1.31319e-05,
                    3.28335e-08,  -1.37002e-07, 6.412e-07,
                    6.32301e-09,  0.000178573,  -7.46352e-07,
                    5.77821e-07,  1.22545e-08,  7.82577e-06,
                    1.31319e-05,  -7.46352e-07, 4.27034e-05,
                    -2.44549e-13, -2.98712e-09, 5.95667e-09,
                    3.28335e-08,  5.77821e-07,  -2.44549e-13,
                    8.26179e-08,  8.02087e-11,  2.53174e-08,
                    -1.37002e-07, 1.22545e-08,  -2.98712e-09,
                    8.02087e-11,  1.36315e-06,  -1.7853e-06,
                    6.412e-07,    7.82577e-06,  5.95667e-09,
                    2.53174e-08,  -1.7853e-06,  89875.5}),  // param cov
               makeVector<ActsScalar, 3>(
                   {-5.04066e-05 * 1000_T, -1.84572e-06 * 1000_T,
                    0.00107321 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-2.79072, 18.1615, 1962.71}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {-0.986831, -0.161755, -2.38917e-17, -0.161755, 0.986831,
               -2.35312e-18, 2.39577e-17, 1.54245e-18, -1}),  // surface rot
          TestData::kPolarDisk,
          makeVector<ActsScalar, eBoundSize>({874.522, -0.0199525, -2.87012,
                                              0.412785, -0.218474,
                                              2144.77}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.268052,     6.23529e-06,  -2.89316e-05,
               0.000334472,  6.69436e-06,  0.107219,
               6.23529e-06,  4.51554e-10,  -5.00139e-09,
               7.5944e-09,   1.32896e-09,  2.47693e-06,
               -2.89316e-05, -5.00139e-09, 1.10493e-06,
               -8.28542e-08, -1.11202e-07, -1.05201e-05,
               0.000334472,  7.5944e-09,   -8.28542e-08,
               7.43596e-07,  3.40043e-08,  0.00013338,
               6.69436e-06,  1.32896e-09,  -1.11202e-07,
               3.40043e-08,  2.28008e-06,  -1.3933e-05,
               0.107219,     2.47693e-06,  -1.05201e-05,
               0.00013338,   -1.3933e-05,  89875.6}),  // param cov
          makeVector<ActsScalar, 3>({-0.000238594 * 1000_T,
                                     -3.95897e-05 * 1000_T,
                                     0.00170904 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-1.04461, -18.345, -2232.72}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {-0.997764, 0.0668313, 1.22465e-16, 0.0668313, 0.997764,
               1.22465e-16, -1.14006e-16, 1.30375e-16, -1}),  // surface rot
          TestData::kPolarDisk,
          makeVector<ActsScalar, eBoundSize>({919.923, -0.0334805, 3.06771,
                                              2.75516, 0.106936,
                                              2410.21}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.249495,     -5.22193e-06, -3.13362e-06,
               -0.000237539, -8.44559e-08, 0.0938443,
               -5.22193e-06, 4.17316e-10,  -3.62456e-09,
               4.90101e-09,  2.09791e-10,  -1.95858e-06,
               -3.13362e-06, -3.62456e-09, 9.69964e-07,
               -1.38253e-08, -2.38732e-08, -1.43246e-06,
               -0.000237539, 4.90101e-09,  -1.38253e-08,
               4.54156e-07,  3.99992e-09,  -8.93313e-05,
               -8.44559e-08, 2.09791e-10,  -2.38732e-08,
               3.99992e-09,  5.6852e-07,   4.62431e-06,
               0.0938443,    -1.95858e-06, -1.43246e-06,
               -8.93313e-05, 4.62431e-06,  89875.6}),  // param cov
          makeVector<ActsScalar, 3>({0.00036598 * 1000_T, -5.73626e-06 * 1000_T,
                                     0.0015599 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {2.25897, 16.0914, 1962.71}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {-0.99163, 0.129111, 2.38917e-17, 0.129111, 0.99163, -2.35312e-18,
               -2.39955e-17, 7.51254e-19, -1}),  // surface rot
          TestData::kPolarDisk,
          makeVector<ActsScalar, eBoundSize>({855.523, -0.00206235, 2.78002,
                                              0.430255, 0.430811,
                                              2164.93}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.46359,      9.75743e-06,  0.000448092,
               0.000688107,  -0.000233571, 0.188268,
               9.75743e-06,  6.10067e-10,  -4.97007e-10,
               1.52396e-08,  3.93299e-09,  4.13675e-06,
               0.000448092,  -4.97007e-10, 4.44146e-06,
               8.69542e-07,  -2.28873e-06, 0.00014712,
               0.000688107,  1.52396e-08,  8.69542e-07,
               2.0094e-06,   -8.25987e-07, 0.000271189,
               -0.000233571, 3.93299e-09,  -2.28873e-06,
               -8.25987e-07, 4.58208e-05,  0.000645024,
               0.188268,     4.13675e-06,  0.00014712,
               0.000271189,  0.000645024,  89875.6}),  // param cov
          makeVector<ActsScalar, 3>({-0.000235096 * 1000_T,
                                     3.37809e-05 * 1000_T,
                                     0.00170337 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-956.977, -300.445, -563.272}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {0.113165, 0.00294296, -0.993572, -0.993236, -0.02583, -0.113203,
               -0.0259971, 0.999662, -9.95799e-17}),  // surface rot
          TestData::kPlane,
          makeVector<ActsScalar, eBoundSize>({43.1027, -20.0508, -2.58378,
                                              2.04794, -0.61896,
                                              1281.4}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.000415869,  -0.00624916,  4.30188e-06,
               1.11988e-05,  -8.34103e-07, 0.00297658,
               -0.00624916,  0.287414,     0.000262654,
               -0.000530088, -1.45312e-05, -0.131126,
               4.30188e-06,  0.000262654,  2.57536e-06,
               -4.80124e-07, -2.45959e-07, -0.000110478,
               1.11988e-05,  -0.000530088, -4.80124e-07,
               3.87906e-06,  6.73799e-08,  0.0002415,
               -8.34103e-07, -1.45312e-05, -2.45959e-07,
               6.73799e-08,  4.23156e-06,  -4.9642e-05,
               0.00297658,   -0.131126,    -0.000110478,
               0.0002415,    -4.9642e-05,  89875.6}),  // param cov
          makeVector<ActsScalar, 3>({3.50239e-05 * 1000_T, 1.16125e-05 * 1000_T,
                                     0.00201393 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-14.6729, -11.0605, 2860.72}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {0.609896, -0.792481, -1.01826e-16, -0.792481, -0.609896,
               -1.90502e-16, 8.88665e-17, 1.96882e-16, -1}),  // surface rot
          TestData::kPolarDisk,
          makeVector<ActsScalar, eBoundSize>({878.179, 0.0377489, -0.917117,
                                              0.298851, -0.155266,
                                              2993.62}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.247492,     5.75149e-06,  -7.25055e-06,
               0.000186384,  2.26831e-05,  0.0724092,
               5.75149e-06,  4.34292e-10,  -3.3005e-09,
               4.28514e-09,  1.10915e-10,  1.68264e-06,
               -7.25055e-06, -3.3005e-09,  5.13085e-07,
               -2.15014e-08, 1.38502e-07,  -3.13833e-06,
               0.000186384,  4.28514e-09,  -2.15014e-08,
               2.64484e-07,  4.58027e-08,  5.43642e-05,
               2.26831e-05,  1.10915e-10,  1.38502e-07,
               4.58027e-08,  5.0017e-06,   -3.05689e-05,
               0.0724092,    1.68264e-06,  -3.13833e-06,
               5.43642e-05,  -3.05689e-05, 89875.5}),  // param cov
          makeVector<ActsScalar, 3>({0.000246188 * 1000_T,
                                     -0.000361053 * 1000_T,
                                     0.000756348 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-748.018, 127.04, 1249.27}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {-0.368695, 0.00958824, -0.929501, -0.929187, 0.0241643, 0.36882,
               0.0259971, 0.999662, -2.22045e-16}),  // surface rot
          TestData::kPlane,
          makeVector<ActsScalar, eBoundSize>({-12.7967, -0.0111021, 2.81269,
                                              0.548845, 0.365828,
                                              1466.51}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.000407068, 0.00585003,  -5.06759e-06, -1.4351e-06, 3.49367e-07,
               0.00501982,  0.00585003,  0.220933,     -5.8514e-05, 7.21617e-06,
               3.69039e-05, 0.189194,    -5.06759e-06, -5.8514e-05, 6.93556e-07,
               3.87361e-07, 2.09743e-07, -4.8068e-05,  -1.4351e-06, 7.21617e-06,
               3.87361e-07, 1.61864e-06, -3.89113e-08, 5.41601e-06, 3.49367e-07,
               3.69039e-05, 2.09743e-07, -3.89113e-08, 3.38809e-06, 6.8792e-05,
               0.00501982,  0.189194,    -4.8068e-05,  5.41601e-06, 6.8792e-05,
               89875.7}),  // param cov
          makeVector<ActsScalar, 3>({-8.32767e-05 * 1000_T,
                                     1.42907e-05 * 1000_T,
                                     0.00191073 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {-522.319, -215.962, 1139.19}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {0.182174, 0.00473759, -0.983255, -0.982923, -0.0255617,
               -0.182236, -0.0259971, 0.999662, -2.22045e-16}),  // surface rot
          TestData::kPlane,
          makeVector<ActsScalar, eBoundSize>({-15.9238, -9.71785, -2.58798,
                                              0.458544, -0.490701,
                                              1263.01}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.000388298,  -0.00663415, -4.45582e-06, -2.98583e-06,
               -1.79271e-07, -0.00594178, -0.00663415,  0.28281,
               0.000192569,  3.50465e-05, 4.73666e-05,  0.253947,
               -4.45582e-06, 0.000192569, 4.58913e-07,  -3.91615e-07,
               3.05147e-07,  0.000170095, -2.98583e-06, 3.50465e-05,
               -3.91615e-07, 9.06356e-07, 1.00807e-08,  3.10251e-05,
               -1.79271e-07, 4.73666e-05, 3.05147e-07,  1.00807e-08,
               3.51684e-06,  4.83158e-06, -0.00594178,  0.253947,
               0.000170095,  3.10251e-05, 4.83158e-06,  89875.7}),  // param cov
          makeVector<ActsScalar, 3>({-5.25166e-05 * 1000_T,
                                     -2.12894e-05 * 1000_T,
                                     0.00191577 * 1000_T})),  // magnetic field

      TestData(
          makeVector<ActsScalar, 3>(
              {889.91, 463.263, 269.272}),  // surface center
          makeMatrix<ActsScalar, 3, 3>(
              {-0.28392, -0.00738357, 0.95882, 0.958496, 0.0249265, 0.284016,
               -0.0259971, 0.999662, -2.22045e-16}),  // surface rot
          TestData::kPlane,
          makeVector<ActsScalar, eBoundSize>({-31.3391, 22.7156, 0.668348,
                                              1.23223, -0.672912,
                                              887.012}),  // param_vec
          makeMatrix<ActsScalar, eBoundSize, eBoundSize>(
              {0.000401596,  -0.00580898,  3.84845e-06,
               1.07378e-05,  -7.00551e-06, -0.00176611,
               -0.00580898,  0.260857,     0.000203293,
               -0.000491391, 1.75445e-05,  0.0869615,
               3.84845e-06,  0.000203293,  2.13361e-06,
               -3.99989e-07, -1.30627e-06, 8.83175e-05,
               1.07378e-05,  -0.000491391, -3.99989e-07,
               3.41664e-06,  1.46081e-07,  -0.000166376,
               -7.00551e-06, 1.75445e-05,  -1.30627e-06,
               1.46081e-07,  2.06909e-05,  -0.000261632,
               -0.00176611,  0.0869615,    8.83175e-05,
               -0.000166376, -0.000261632, 89875.6}),  // param cov
          makeVector<ActsScalar, 3>({1.43009e-05 * 1000_T, 7.04031e-06 * 1000_T,
                                     0.00202663 * 1000_T}))  // magnetic field

  };
  if (magFieldScale != 1.) {
    for (TestData &testData : testDataList) {
      testData.bfield *= magFieldScale;
    }
  }
  return testDataList;
}

BOOST_AUTO_TEST_CASE(BoundToCurvilinearEigenStepper) {
  // Compare covariance in curvilinear parameterisation:
  // explicit computation vs. dummy propagation using EigenStepper
  std::vector<TestData> testDataList(make_test_data());
  test_bound_to_curvilinear(
      testDataList,
      [](const std::shared_ptr<Acts::MagneticFieldProvider> &bField) {
        return EigenStepper<>(bField);
      });
}

BOOST_AUTO_TEST_CASE(BoundToCurvilinearStraightLineStepper) {
  // Compare covariance in curvilinear parameterisation for vanishing magnetic
  // field: explicit computation vs. dummy propagation using StraightLineStepper
  std::vector<TestData> testDataList(
      make_test_data(0.));  // scale magnetic field in test data to zero.
  test_bound_to_curvilinear(
      testDataList,
      []([[maybe_unused]] const std::shared_ptr<Acts::MagneticFieldProvider>
             &bField) { return StraightLineStepper(); });
}
