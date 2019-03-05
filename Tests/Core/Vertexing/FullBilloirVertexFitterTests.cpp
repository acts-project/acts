// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE FullBilloirVertexFitter Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/IVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

namespace bdata = boost::unit_test::data;

namespace Acts {
namespace Test {

  template <typename InputTrack_t, typename Propagator_t>
  Vertex<InputTrack_t>
  myFitWrapper(IVertexFitter<InputTrack_t, Propagator_t>* fitter,
               std::vector<InputTrack_t>& tracks,
               const Propagator_t&        propagator,
               Vertex<InputTrack_t>*      constraint = nullptr)
  {
    if (constraint != nullptr) {
      return fitter->fit(tracks, propagator, *constraint);
    } else {
      return fitter->fit(tracks, propagator);
    }
  }

  ///
  /// @brief Unit test for FullBilloirVertexFitter
  ///
  BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_empty_input_test)
  {

    // Set up constant B-Field
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    Propagator<EigenStepper<ConstantBField>> propagator(stepper);

    // Set up Billoir Vertex Fitter
    FullBilloirVertexFitter<ConstantBField,
                            BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>::Config
        vertexFitterCfg(bField);
    FullBilloirVertexFitter<ConstantBField,
                            BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>
        billoirFitter(vertexFitterCfg);

    // Constraint for vertex fit
    Vertex<BoundParameters> myConstraint;
    // Some abitrary values
    ActsSymMatrixD<3> myCovMat = ActsSymMatrixD<3>::Zero();
    myCovMat(0, 0) = 30.;
    myCovMat(1, 1) = 30.;
    myCovMat(2, 2) = 30.;
    myConstraint.setCovariance(std::move(myCovMat));
    myConstraint.setPosition(Vector3D(0, 0, 0));

    std::vector<BoundParameters> emptyVector;

    Vertex<BoundParameters> fittedVertex
        = billoirFitter.fit(emptyVector, propagator, myConstraint);
    Vector3D origin(0., 0., 0.);
    BOOST_CHECK_EQUAL(fittedVertex.position(), origin);

    ActsSymMatrixD<3> zeroMat = ActsSymMatrixD<3>::Zero();
    BOOST_CHECK_EQUAL(fittedVertex.covariance(), zeroMat);

    fittedVertex = billoirFitter.fit(emptyVector, propagator);
    BOOST_CHECK_EQUAL(fittedVertex.position(), origin);
    BOOST_CHECK_EQUAL(fittedVertex.covariance(), zeroMat);
  }

  // Vertex x/y position distribution
  std::uniform_real_distribution<> vXYDist(-0.1 * units::_mm, 0.1 * units::_mm);
  // Vertex z position distribution
  std::uniform_real_distribution<> vZDist(-20 * units::_mm, 20 * units::_mm);
  // Track d0 distribution
  std::uniform_real_distribution<> d0Dist(-0.01 * units::_mm,
                                          0.01 * units::_mm);
  // Track z0 distribution
  std::uniform_real_distribution<> z0Dist(-0.2 * units::_mm, 0.2 * units::_mm);
  // Track pT distribution
  std::uniform_real_distribution<> pTDist(0.4 * units::_GeV, 10. * units::_GeV);
  // Track phi distribution
  std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
  // Track theta distribution
  std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
  // Track charge helper distribution
  std::uniform_real_distribution<> qDist(-1, 1);
  // Track IP resolution distribution
  std::uniform_real_distribution<> resIPDist(0., 100. * units::_um);
  // Track angular distribution
  std::uniform_real_distribution<> resAngDist(0., 0.1);
  // Track q/p resolution distribution
  std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
  // Number of tracks distritbution
  std::uniform_int_distribution<> nTracksDist(3, 10);

  ///
  /// @brief Unit test for FullBilloirVertexFitter
  /// with default input track type (= BoundParameters)
  ///
  BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_defaulttrack_test)
  {

    bool debugMode = false;

    // Set up RNG
    std::random_device rd;
    std::mt19937       gen(rd());

    // Set up constant B-Field
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    Propagator<EigenStepper<ConstantBField>> propagator(stepper);

    // Number of events
    const int nEvents = 10;

    for (int eventIdx = 0; eventIdx < nEvents; ++eventIdx) {

      // Number of tracks
      unsigned int nTracks = nTracksDist(gen);

      // Set up Billoir Vertex Fitter
      FullBilloirVertexFitter<ConstantBField,
                              BoundParameters,
                              Propagator<EigenStepper<ConstantBField>>>::Config
          vertexFitterCfg(bField);
      FullBilloirVertexFitter<ConstantBField,
                              BoundParameters,
                              Propagator<EigenStepper<ConstantBField>>>
          billoirFitter(vertexFitterCfg);

      // Constraint for vertex fit
      Vertex<BoundParameters> myConstraint;
      // Some abitrary values
      ActsSymMatrixD<3> myCovMat = ActsSymMatrixD<3>::Zero();
      myCovMat(0, 0) = 30.;
      myCovMat(1, 1) = 30.;
      myCovMat(2, 2) = 30.;
      myConstraint.setCovariance(std::move(myCovMat));
      myConstraint.setPosition(Vector3D(0, 0, 0));

      // Create position of vertex and perigee surface
      double x = vXYDist(gen);
      double y = vXYDist(gen);
      double z = vZDist(gen);

      Vector3D                        vertexPosition(x, y, z);
      std::shared_ptr<PerigeeSurface> perigeeSurface
          = Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

      // Calculate d0 and z0 corresponding to vertex position
      double d0V = sqrt(x * x + y * y);
      double z0V = z;

      // Start constructing nTracks tracks in the following
      std::vector<BoundParameters> tracks;

      // Construct random track emerging from vicinity of vertex position
      // Vector to store track objects used for vertex fit
      for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
        // Construct positive or negative charge randomly
        double q = qDist(gen) < 0 ? -1. : 1.;

        // Construct random track parameters
        TrackParametersBase::ParVector_t paramVec;
        paramVec << d0V + d0Dist(gen), z0V + z0Dist(gen), phiDist(gen),
            thetaDist(gen), q / pTDist(gen);

        // Fill vector of track objects with simple covariance matrix
        std::unique_ptr<ActsSymMatrixD<5>> covMat
            = std::make_unique<ActsSymMatrixD<5>>();

        // Resolutions
        double resD0 = resIPDist(gen);
        double resZ0 = resIPDist(gen);
        double resPh = resAngDist(gen);
        double resTh = resAngDist(gen);
        double resQp = resQoPDist(gen);

        (*covMat) << resD0 * resD0, 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
            0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., resTh * resTh, 0.,
            0., 0., 0., 0., resQp * resQp;
        tracks.push_back(
            BoundParameters(std::move(covMat), paramVec, perigeeSurface));
      }

      // Do the actual fit with 4 tracks without constraint
      Vertex<BoundParameters> fittedVertex
          = billoirFitter.fit(tracks, propagator);
      if (fittedVertex.tracks().size() > 0) {
        CHECK_CLOSE_ABS(
            fittedVertex.position(), vertexPosition, 1 * units::_mm);
      }
      // Do the fit with a constraint
      Vertex<BoundParameters> fittedVertexConstraint
          = billoirFitter.fit(tracks, propagator, myConstraint);
      if (fittedVertexConstraint.tracks().size() > 0) {
        CHECK_CLOSE_ABS(
            fittedVertexConstraint.position(), vertexPosition, 1 * units::_mm);
      }
      // Test the IVertexFitter interface
      Vertex<BoundParameters> testVertex
          = myFitWrapper(&billoirFitter, tracks, propagator);
      if (testVertex.tracks().size() > 0) {
        CHECK_CLOSE_ABS(testVertex.position(), vertexPosition, 1 * units::_mm);
      }

      if (debugMode) {
        std::cout << "Fitting nTracks: " << nTracks << std::endl;
        std::cout << "True Vertex: " << x << ", " << y << ", " << z
                  << std::endl;
        std::cout << "Fitted Vertex: " << fittedVertex.position() << std::endl;
        std::cout << "Fitted constraint Vertex: "
                  << fittedVertexConstraint.position() << std::endl;
      }
    }
  }

  // Dummy user-defined InputTrack type
  struct InputTrack
  {
    InputTrack(const BoundParameters& params) : m_parameters(params) {}

    const BoundParameters&
    parameters() const
    {
      return m_parameters;
    }

    // store e.g. link to original objects here

  private:
    BoundParameters m_parameters;
  };

  ///
  /// @brief Unit test for FullBilloirVertexFitter with user-defined InputTrack
  /// type
  ///
  BOOST_AUTO_TEST_CASE(billoir_vertex_fitter_usertrack_test)
  {

    bool debugMode = false;

    // Set up RNG
    std::random_device rd;
    std::mt19937       gen(rd());

    // Set up constant B-Field
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    Propagator<EigenStepper<ConstantBField>> propagator(stepper);

    const int nEvents = 10;

    for (int eventIdx = 0; eventIdx < nEvents; ++eventIdx) {

      unsigned int nTracks = nTracksDist(gen);

      // Create a custom std::function to extract BoundParameters from
      // user-defined InputTrack
      std::function<BoundParameters(InputTrack)> extractParameters
          = [](InputTrack params) { return params.parameters(); };

      // Set up Billoir Vertex Fitter
      FullBilloirVertexFitter<ConstantBField,
                              InputTrack,
                              Propagator<EigenStepper<ConstantBField>>>::Config
          vertexFitterCfg(bField);
      FullBilloirVertexFitter<ConstantBField,
                              InputTrack,
                              Propagator<EigenStepper<ConstantBField>>>
          billoirFitter(vertexFitterCfg, extractParameters);

      // Constraint for vertex fit
      Vertex<InputTrack> myConstraint;
      // Some abitrary values
      ActsSymMatrixD<3> myCovMat = ActsSymMatrixD<3>::Zero();
      myCovMat(0, 0) = 30.;
      myCovMat(1, 1) = 30.;
      myCovMat(2, 2) = 30.;
      myConstraint.setCovariance(std::move(myCovMat));
      myConstraint.setPosition(Vector3D(0, 0, 0));

      // Create position of vertex and perigee surface
      double x = vXYDist(gen);
      double y = vXYDist(gen);
      double z = vZDist(gen);

      Vector3D                        vertexPosition(x, y, z);
      std::shared_ptr<PerigeeSurface> perigeeSurface
          = Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

      // Calculate d0 and z0 corresponding to vertex position
      double d0V = sqrt(x * x + y * y);
      double z0V = z;

      // Start constructing nTracks tracks in the following
      std::vector<InputTrack> tracks;

      // Construct random track emerging from vicinity of vertex position
      // Vector to store track objects used for vertex fit
      for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
        // Construct positive or negative charge randomly
        double q = qDist(gen) < 0 ? -1. : 1.;

        // Construct random track parameters
        TrackParametersBase::ParVector_t paramVec;
        paramVec << d0V + d0Dist(gen), z0V + z0Dist(gen), phiDist(gen),
            thetaDist(gen), q / pTDist(gen);

        // Fill vector of track objects with simple covariance matrix
        std::unique_ptr<ActsSymMatrixD<5>> covMat
            = std::make_unique<ActsSymMatrixD<5>>();

        // Resolutions
        double resD0 = resIPDist(gen);
        double resZ0 = resIPDist(gen);
        double resPh = resAngDist(gen);
        double resTh = resAngDist(gen);
        double resQp = resQoPDist(gen);

        (*covMat) << resD0 * resD0, 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
            0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., resTh * resTh, 0.,
            0., 0., 0., 0., resQp * resQp;
        tracks.push_back(InputTrack(
            BoundParameters(std::move(covMat), paramVec, perigeeSurface)));
      }

      // Do the actual fit with 4 tracks without constraint
      Vertex<InputTrack> fittedVertex = billoirFitter.fit(tracks, propagator);
      if (fittedVertex.tracks().size() > 0) {
        CHECK_CLOSE_ABS(
            fittedVertex.position(), vertexPosition, 1 * units::_mm);
      }
      // Do the fit with a constraint
      Vertex<InputTrack> fittedVertexConstraint
          = billoirFitter.fit(tracks, propagator, myConstraint);
      if (fittedVertexConstraint.tracks().size() > 0) {
        CHECK_CLOSE_ABS(
            fittedVertexConstraint.position(), vertexPosition, 1 * units::_mm);
      }
      // Test the IVertexFitter interface
      Vertex<InputTrack> testVertex
          = myFitWrapper(&billoirFitter, tracks, propagator);
      if (testVertex.tracks().size() > 0) {
        CHECK_CLOSE_ABS(testVertex.position(), vertexPosition, 1 * units::_mm);
      }
      if (debugMode) {
        std::cout << "Fitting nTracks: " << nTracks << std::endl;
        std::cout << "True Vertex: " << x << ", " << y << ", " << z
                  << std::endl;
        std::cout << "Fitted Vertex: " << fittedVertex.position() << std::endl;
        std::cout << "Fitted constraint Vertex: "
                  << fittedVertexConstraint.position() << std::endl;
      }
    }
  }

}  // namespace Test
}  // namespace Acts
