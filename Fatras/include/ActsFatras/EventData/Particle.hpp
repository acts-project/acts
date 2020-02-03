// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

namespace ActsFatras {

/// Typedef the pdg code
typedef int pdg_type;

/// Typedef the process code
typedef unsigned int process_code;

/// Typedef barcode
typedef unsigned int barcode_type;

/// Simulation particle information and kinematic state.
class Particle {
 public:
  /// Default
  Particle() = default;

  /// @brief Construct a particle consistently
  ///
  /// @param pposition The particle position at construction
  /// @param pmomentum The particle momentum at construction
  /// @param pm The particle mass
  /// @param pq The partilce charge
  /// @param pbarcode The particle barcode
  Particle(const Acts::Vector3D &position, const Acts::Vector3D &momentum,
           double m, double q, pdg_type pdg = 0, barcode_type barcode = 0,
           double startTime = 0.)
      : m_position(position),
        m_momentum(momentum),
        m_m(m),
        m_q(q),
        m_p(momentum.norm()),
        m_pT(Acts::VectorHelpers::perp(momentum)),
        m_pdg(pdg),
        m_barcode(barcode),
        m_timeStamp(startTime) {
    m_E = std::sqrt(m_p * m_p + m_m * m_m);
    m_beta = (m_p / m_E);
    m_gamma = (m_E / m_m);
  }

  /// @brief Set the limits
  ///
  /// @param x0Limit the limit in X0 to be passed
  /// @param l0Limit the limit in L0 to be passed
  /// @param timeLimit the readout time limit to be passed
  void setLimits(double x0Limit, double l0Limit,
                 double timeLimit = std::numeric_limits<double>::max()) {
    m_limitInX0 = x0Limit;
    m_limitInL0 = l0Limit;
    m_timeLimit = timeLimit;
  }

  /// @brief Update the particle with applying energy loss
  ///
  /// @param deltaE is the energy loss to be applied
  void scatter(Acts::Vector3D nmomentum) {
    m_momentum = std::move(nmomentum);
    m_pT = Acts::VectorHelpers::perp(m_momentum);
  }

  /// @brief Update the particle with applying energy loss
  ///
  /// @param deltaE is the energy loss to be applied
  void energyLoss(double deltaE) {
    // particle falls to rest
    if (m_E - deltaE < m_m) {
      m_E = m_m;
      m_p = 0.;
      m_pT = 0.;
      m_beta = 0.;
      m_gamma = 1.;
      m_momentum = Acts::Vector3D(0., 0., 0.);
      m_alive = false;
    }
    // updatet the parameters
    m_E -= deltaE;
    m_p = std::sqrt(m_E * m_E - m_m * m_m);
    m_momentum = m_p * m_momentum.normalized();
    m_pT = Acts::VectorHelpers::perp(m_momentum);
    m_beta = (m_p / m_E);
    m_gamma = (m_E / m_m);
  }

  /// @brief Update the particle with a new position and momentum,
  /// this corresponds to a step update
  ///
  /// @param position New position after update
  /// @param momentum New momentum after update
  /// @param deltaPathX0 passed since last step
  /// @param deltaPathL0 passed since last step
  /// @param deltaTime The time elapsed
  ///
  /// @return break condition
  bool update(const Acts::Vector3D &position, const Acts::Vector3D &momentum,
              double deltaPahtX0 = 0., double deltaPahtL0 = 0.,
              double deltaTime = 0.) {
    m_position = position;
    m_momentum = momentum;
    m_p = momentum.norm();
    if (m_p) {
      m_pT = Acts::VectorHelpers::perp(momentum);
      m_E = std::sqrt(m_p * m_p + m_m * m_m);
      m_timeStamp += deltaTime;
      m_beta = (m_p / m_E);
      m_gamma = (m_E / m_m);

      // set parameters and check limits
      m_pathInX0 += deltaPahtX0;
      m_pathInL0 += deltaPahtL0;
      m_timeStamp += deltaTime;
      if (m_pathInX0 >= m_limitInX0 || m_pathInL0 >= m_limitInL0 ||
          m_timeStamp > m_timeLimit) {
        m_alive = false;
      }
    }
    return !m_alive;
  }

  /// @brief Access methods: position
  const Acts::Vector3D &position() const { return m_position; }

  /// @brief Access methods: momentum
  const Acts::Vector3D &momentum() const { return m_momentum; }

  /// @brief Access methods: p
  double p() const { return m_p; }

  /// @brief Access methods: pT
  double pT() const { return m_pT; }

  /// @brief Access methods: E
  double E() const { return m_E; }

  /// @brief Access methods: m
  double m() const { return m_m; }

  /// @brief Access methods: beta
  double beta() const { return m_beta; }

  /// @brief Access methods: gamma
  double gamma() const { return m_gamma; }

  /// @brief Access methods: charge
  double q() const { return m_q; }

  /// @brief Access methods: pdg code
  pdg_type pdg() const { return m_pdg; }

  /// @brief Access methods: barcode
  barcode_type barcode() const { return m_barcode; }

  /// @brief Access methods: path/X0
  double pathInX0() const { return m_pathInX0; }

  /// @brief Access methods: limit/X0
  double limitInX0() const { return m_limitInX0; }

  /// @brief Access methods: pdg code
  double pathInL0() const { return m_limitInX0; }

  /// @brief Access methods: barcode
  double limitInL0() const { return m_limitInL0; }

  /// @brief boolean operator indicating the particle to be alive
  operator bool() { return m_alive; }

 private:
  Acts::Vector3D m_position = Acts::Vector3D(0., 0., 0.);  //!< kinematic info
  Acts::Vector3D m_momentum = Acts::Vector3D(0., 0., 0.);  //!< kinematic info

  double m_m = 0.;             //!< particle mass
  double m_E = 0.;             //!< total energy
  double m_q = 0.;             //!< the charge
  double m_beta = 0.;          //!< relativistic beta factor
  double m_gamma = 1.;         //!< relativistic gamma factor
  double m_p = 0.;             //!< momentum magnitude
  double m_pT = 0.;            //!< transverse momentum magnitude
  pdg_type m_pdg = 0;          //!< pdg code of the particle
  barcode_type m_barcode = 0;  //!< barcode of the particle

  double m_pathInX0 = 0.;  //!< passed path in X0
  double m_limitInX0 =
      std::numeric_limits<double>::max();  //!< path limit in X0

  double m_pathInL0 = 0.;  //!< passed path in L0
  double m_limitInL0 =
      std::numeric_limits<double>::max();  //!< path limit in X0

  double m_timeStamp = 0.;  //!< passed time elapsed
  double m_timeLimit = std::numeric_limits<double>::max();  // time limit

  bool m_alive = true;  //!< the particle is alive
};

}  // namespace ActsFatras
