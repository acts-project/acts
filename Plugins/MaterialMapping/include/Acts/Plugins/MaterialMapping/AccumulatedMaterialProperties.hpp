// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedMaterialProperties.hpp, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

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
/// - the event store collects material steps accumulated
///   during an event
/// - the total store colles accumulated material properties
///   of the run
class AccumulatedMaterialProperties
{

public:
  /// Default constructor sets everything to zero
  AccumulatedMaterialProperties() = default;

  /// Default destructor sets
  ~AccumulatedMaterialProperties() = default;

  /// Default copy constructor
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties(const AccumulatedMaterialProperties& amp)
      = default;

  /// Default copy move constructor
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties(AccumulatedMaterialProperties&& amp) = default;

  /// Default assignment operator
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties&
  operator=(const AccumulatedMaterialProperties& amp)
      = default;

  /// Default move assignment operator
  /// @param amp the source Accumulated properties
  AccumulatedMaterialProperties&
  operator=(AccumulatedMaterialProperties&& amp)
      = default;

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
  void
  operator+=(const MaterialProperties& amp);

  /// Average the information accumulated during one event
  /// using the event weights
  void
  eventAverage();

  /// Average the information accumulated during the entire
  /// mapping process
  ///
  /// The total average takes
  ///
  /// @returns the  total avarage material properties and the
  /// total events used for mapping these properties
  /// The latter can be used when parallelising the mapping
  /// to several jobs
  ///
  /// @param unitThickness sets the material to unit thickness
  std::pair<MaterialProperties, unsigned int>
  totalAverage(bool unitThickness = true);

private:
  double m_eventA{0.};     //!< event: accumulate the contribution to A
  double m_eventZ{0.};     //!< event: accumulate the contribution to Z
  double m_eventRho{0.};   //!< event: accumulate the contribution to rho
  double m_eventPath{0.};  //!< event: the event path for normalisation

  double m_totalPathInX0{0.};  //!< total: accumulate the thickness in X0
  double m_totalPathInL0{0.};  //!< total: accumulate the thickness in L0
  double m_totalA{0.};         //!< total: accumulate the contribution to A
  double m_totalZ{0.};         //!< total: accumulate the contribution to Z
  double m_totalRho{0.};       //!< total: accumulate the contribution to rho
  double m_totalPath{0.};      //!< total: the event path for normalisation

  unsigned int m_totalEvents{0};  //!< the number of events
};

inline void
AccumulatedMaterialProperties::operator+=(const MaterialProperties& amp)
{
  m_totalPathInX0 += amp.thicknessInX0();
  m_totalPathInL0 += amp.thicknessInL0();

  double t = amp.thickness();
  double r = amp.averageRho();
  m_eventPath += t;
  m_eventRho += r * t;

  m_eventA += amp.averageA() * r * t;
  m_eventZ += amp.averageZ() * r * t;
}

inline void
AccumulatedMaterialProperties::eventAverage()
{
  // Always count a hit if a path length is registered
  if (m_eventPath > 0.) {
    m_totalPath += m_eventPath;
    ++m_totalEvents;
  }
  // Average the event quantities
  if (m_eventPath > 0. && m_eventRho > 0.) {
    m_totalA += (m_eventA / m_eventRho);
    m_totalZ += (m_eventZ / m_eventRho);
    m_totalRho += (m_eventRho / m_eventPath);
  }
  m_eventA    = 0.;
  m_eventZ    = 0.;
  m_eventRho  = 0.;
  m_eventPath = 0.;
}

inline std::pair<MaterialProperties, unsigned int>
AccumulatedMaterialProperties::totalAverage(bool unitThickness)
{
  if (m_totalEvents && m_totalPathInX0 > 0.) {
    double eventScalor = 1. / m_totalEvents;
    m_totalPathInX0 *= eventScalor;
    m_totalPathInL0 *= eventScalor;
    m_totalRho *= eventScalor;
    m_totalA *= eventScalor;
    m_totalZ *= eventScalor;
    m_totalPath *= eventScalor;
    // create the material
    double X0 = m_totalPath / m_totalPathInX0;
    double L0 = m_totalPath / m_totalPathInL0;

    // create the material properties
    MaterialProperties averageMat(
        X0, L0, m_totalA, m_totalZ, m_totalRho, m_totalPath);
    if (unitThickness) averageMat.scaleToUnitThickness();

    return std::pair<MaterialProperties, unsigned int>(std::move(averageMat),
                                                       m_totalEvents);
  }
  return std::pair<MaterialProperties, unsigned int>(
      MaterialProperties{Material(), 0.}, 0);
}

}  // end of namespace Acts
