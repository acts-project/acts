// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialProperties.hpp"

namespace Acts {

/// @brief AccumulatedMaterialProperties
///
/// This class uses the  MaterialProperties class and adds
/// the necessary information for the SurfaceMaterial
/// Mapper.
///
/// There are two stores for material averaging:
///
/// - the event store collects material steps accumulated
///   during an event. Material can be assigned to the same bin
///   multiple time by one particle, e.g. if the simulation had created
///   more than one step in the material, or if several components are
///   compressed into one description
///
/// - the total store collects accumulated material properties
///   of the run, at the end of a run, an average over the material
///   information from all mapped events per bin is taken.
///
/// The averaging is always done to unit thickness
class AccumulatedMaterialProperties {
 public:
  /// Default constructor sets everything to zero
  AccumulatedMaterialProperties() = default;

  /// Default destructor sets
  ~AccumulatedMaterialProperties() = default;

  /// Default copy constructor
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties(const AccumulatedMaterialProperties& amp) =
      default;

  /// Default copy move constructor
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties(AccumulatedMaterialProperties&& amp) = default;

  /// Default assignment operator
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties& operator=(
      const AccumulatedMaterialProperties& amp) = default;

  /// Default move assignment operator
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties& operator=(
      AccumulatedMaterialProperties&& amp) = default;

  /// Accumulation operator
  /// @brief this adds the material properties to event store
  ///
  /// MaterialProperties accumulation is done as follows:
  ///   * t/X0 and l/X0 are additive
  ///   * rho is averaged by the total thickness per event
  ///   * A,Z is averaged by the sum(rho*thickness) per event
  ///
  /// Vacuum steps with step length different from zero are
  /// also used in order to count for holes in material structures
  ///
  /// @param amp the source Accumulated properties
  /// @param pathCorreciton is the correction to nominal incident
  void accumulate(const MaterialProperties& amp, double pathCorrection = 1.);

  /// Average the information accumulated for one track
  /// - in case of an empty hit this will just increase the event counter
  /// @param emtpyHit indicate an empty hit
  void trackAverage(bool emptyHit = false);

  /// Average the information accumulated during the entire
  /// mapping process
  ///
  /// The total average takesthe full run into account
  ///
  /// @returns the  total avarage material properties and the
  /// total events used for mapping these properties
  /// The latter can be used when parallelising the mapping
  /// to several jobs
  ///
  std::pair<MaterialProperties, unsigned int> totalAverage();

 private:
  double m_eventPathInX0{0.};  //!< event: accumulate the thickness in X0
  double m_eventPathInL0{0.};  //!< event: accumulate the thickness in L0
  double m_eventAr{0.};        //!< event: accumulate the contribution to A
  double m_eventZ{0.};         //!< event: accumulate the contribution to Z
  double m_eventRho{0.};       //!< event: accumulate the contribution to rho
  double m_eventPath{0.};      //!< event: the event path for normalisation
  double m_eventPathCorrection{0.};  //!< event: remember the path correction

  double m_totalPathInX0{0.};  //!< total: accumulate the thickness in X0
  double m_totalPathInL0{0.};  //!< total: accumulate the thickness in L0
  double m_totalAr{0.};        //!< total: accumulate the contribution to A
  double m_totalZ{0.};         //!< total: accumulate the contribution to Z
  double m_totalRho{0.};       //!< total: accumulate the contribution to rho

  unsigned int m_totalEvents{0};  //!< the number of events
};

inline void AccumulatedMaterialProperties::accumulate(
    const MaterialProperties& amp, double pathCorrection) {
  m_eventPathInX0 += amp.thicknessInX0();
  m_eventPathInL0 += amp.thicknessInL0();
  double t = amp.thickness();
  double r = amp.material().massDensity();
  m_eventPath += t;
  m_eventRho += r * t;

  m_eventAr += amp.material().Ar() * r * t;
  m_eventZ += amp.material().Z() * r * t;

  m_eventPathCorrection += pathCorrection * t;
}

inline void AccumulatedMaterialProperties::trackAverage(bool emptyHit) {
  // Increase the total event count either for empty hit or by touched hit
  if (emptyHit || m_eventPath > 0.) {
    ++m_totalEvents;
  }

  // Average the event quantities
  if (m_eventPath > 0. && m_eventRho > 0.) {
    m_eventPathCorrection /= m_eventPath;
    m_totalPathInX0 += m_eventPathInX0 / m_eventPathCorrection;
    m_totalPathInL0 += m_eventPathInL0 / m_eventPathCorrection;
    m_totalAr += (m_eventAr / m_eventRho);
    m_totalZ += (m_eventZ / m_eventRho);
    m_totalRho += (m_eventRho);
  }
  m_eventPathInX0 = 0.;
  m_eventPathInL0 = 0.;
  m_eventAr = 0.;
  m_eventZ = 0.;
  m_eventRho = 0.;
  m_eventPath = 0.;
  m_eventPathCorrection = 0.;
}

inline std::pair<MaterialProperties, unsigned int>
AccumulatedMaterialProperties::totalAverage() {
  if (m_totalEvents > 0 && m_totalPathInX0 > 0.) {
    double eventScalor = 1. / m_totalEvents;
    m_totalPathInX0 *= eventScalor;
    m_totalPathInL0 *= eventScalor;
    m_totalRho *= eventScalor;
    m_totalAr *= eventScalor;
    m_totalZ *= eventScalor;
    // Create the material
    double X0 = 1. / m_totalPathInX0;
    double L0 = 1. / m_totalPathInL0;
    // Create the material properties - fixed to unit path length
    MaterialProperties averageMat(X0, L0, m_totalAr, m_totalZ, m_totalRho, 1.);
    return std::pair<MaterialProperties, unsigned int>(std::move(averageMat),
                                                       m_totalEvents);
  }
  return std::pair<MaterialProperties, unsigned int>(
      MaterialProperties{Material(), 0.}, m_totalEvents);
}

}  // end of namespace Acts
