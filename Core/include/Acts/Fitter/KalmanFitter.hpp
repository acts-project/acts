// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/variant.hpp>
#include <memory>
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief Kalman fitter implementation of Acts
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updator_t Type of the updated class
/// @tparam calibrator_t Type of the calibrator class
/// @tparam input_converter_t Type of the input converter class
/// @tparam output_converter_t Type of the output converter class
///
/// The Kalman forward filter is a subclass of the KalmanFilter and
/// implemented as an Actor to the propagator engine. This takes the
/// measurements and provides it to the KalmanSequencer that has to
/// be part of the navigator.
///
/// Measurements are not required to be ordered for the KalmanFilter,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The Updator is the implemented kalman updator formalism, it
/// runs via a visitor pattern through the measurements.
///
/// The Calibrator is a dedicated calibration algorithm that allows
/// to calibrate measurements using track information, this could be
/// e.g. sagging for wires, module deformations, etc.
///
/// The Input converter is a converter that transforms the input
/// measurement/track/segments into a set of FittableMeasurements
///
/// The Output converter is a converter that transforms the
/// set of track states into a given track/track particle class
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t,
          typename updator_t          = VoidKalmanComponents,
          typename calibrator_t       = VoidKalmanComponents,
          typename input_converter_t  = VoidKalmanComponents,
          typename output_converter_t = VoidKalmanComponents>
class KalmanFitter
{

public:
  /// Constructor from arguments
  ///
  KalmanFitter(propagator_t       pPropagator,
               updator_t          pUpdator    = updator_t(),
               calibrator_t       pCalibrator = calibrator_t(),
               input_converter_t  pInputCnv   = input_converter_t(),
               output_converter_t pOutputCnv  = output_converter_t())
    : m_propagator(std::move(pPropagator))
    , m_updator(std::move(pUpdator))
    , m_calibrator(std::move(pCalibrator))
    , m_inputConverter(std::move(pInputCnv))
    , m_outputConverter(std::move(pOutputCnv))
  {
  }

  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam input_measurements_t Type of the fittable measurements
  /// @tparam parameters_t Type of the initial parameters
  /// @tparam surface_t Type of the reference surface
  ///
  /// @param measurements are the fittable measurements
  /// @param initalParameters is the initial track parameters
  ///
  /// @return the output as an output track
  template <typename input_measurements_t,
            typename parameters_t,
            typename surface_t>
  auto
  fit(const input_measurements_t& measurements,
      const parameters_t& /*initalParameters*/,
      const surface_t* /*pReferenceSurface = nullptr*/) const
  {
    // Bring the measurements into Acts style
    auto trackStates = m_inputConverter(measurements);
    // Run the forward filter as plugin to the Propagator

    // Apply the smoothing

    // Return the converted Track
    return m_outputConverter(trackStates);
  }

private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// The updator for measurement updates
  updator_t m_updator;

  /// The measurement calibrator
  calibrator_t m_calibrator;

  /// The input converter to Fittable measurements
  input_converter_t m_inputConverter;

  /// The output converter into a given format
  output_converter_t m_outputConverter;
};

}  // namespace Acts
