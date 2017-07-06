// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ParticleDefinitions.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EVENTDATA_PARTICLEDEFINITIONS_H
#define ACTS_EVENTDATA_PARTICLEDEFINITIONS_H 1

#include <vector>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

// barcodes
typedef unsigned long barcode_type;
typedef int           process_type;
typedef int           pdg_type;

// define the particle hypotheses
#define PARTICLETYPES 11

namespace Acts {

/// @enum ParticleType
///
/// Enumeration for Particle type respecting
/// the interaction with the material
///
enum ParticleType {

  nonInteracting = 0,  //!< for non-interacting extrapolation
  geantino       = 0,  //!< for non-interacting extrapolation
  electron = 1,  //!< reconstruction + fatras : type as electron hypothesis
  muon     = 2,  //!< reconstruction + fatras : type as muon hypothesis
  pion     = 3,  //!< reconstruction + fatras : type as pion hypothesis
  kaon     = 4,  //!< reconstruction + fatras : type as kaon hypothesis
  proton   = 5,  //!< reconstruction + fatras : type as proton hypothesis
  photon   = 6,  //!< for fatras usage
  neutron  = 7,  //!< for fatras usage
  pi0      = 8,  //!< for fatras usage
  k0       = 9,  //!< for fatras usage
  nonInteractingMuon
  = 10,  //!< for material collection @todo check if this is still needed
  noHypothesis = 99,
  undefined    = 99
};

/// @struct ParticleMasses
/// Mass declaration of particles covered
/// in the ParticleType.
/// Masses are given in MeV and taken from:
///
/// Review of Particle Physics (2010)
/// K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010)
///
struct ParticleMasses
{
  /** the vector of masses - in MeV */
  std::vector<double> mass;

  /**Default constructor*/
  ParticleMasses()
  {
    mass.reserve(PARTICLETYPES);

    mass.push_back(0.);                        // non interacting mass
    mass.push_back(0.51099891 * units::_MeV);  // electron mass
    mass.push_back(105.658367 * units::_MeV);  // muon mass
    mass.push_back(139.57018 * units::_MeV);   // charged pion mass
    mass.push_back(493.677 * units::_MeV);     // kaon mass
    mass.push_back(938.272013 * units::_MeV);  // proton mass
    mass.push_back(0. * units::_MeV);          // photon rest mass
    mass.push_back(939.565346 * units::_MeV);  // neutron rest mass
    mass.push_back(134.9766 * units::_MeV);    // pi0 rest mass
    mass.push_back(497.614 * units::_MeV);     // K0 rest mass
    mass.push_back(105.658367 * units::_MeV);  // muon mass
  }
};

/// Static method : convert to ParticleType from pdg
static ParticleType
convertToParticleType(pdg_type pdg, bool& stable, bool& exiting, double charge)
{
  int pdgCode = abs(pdg);

  stable  = false;
  exiting = false;

  ParticleType particleType;

  // try to follow number of appearance
  switch (pdgCode) {
  // leptons
  case 11:  // e+/e-
    particleType = electron;
    stable       = true;
    exiting      = false;
    break;
  case 13:  // mu+/mu-
    particleType = muon;
    stable       = false;
    exiting      = false;
    break;
  case 12:  // e neutrino
  case 14:  // mu neutrino
  case 16:  // tau neutrino
    particleType = nonInteracting;
    stable       = true;
    exiting      = true;
    break;
  case 22:  // gamma
    particleType = photon;
    stable       = true;
    exiting      = false;
    break;
  case 211:  // pi+/pi-
    particleType = pion;
    stable       = false;
    exiting      = false;
    break;
  case 111:  // pi0
    particleType = pi0;
    stable       = false;
    exiting      = false;
    break;
  case 2212:  // proton
    particleType = proton;
    stable       = true;
    exiting      = false;
    break;
  case 2112:  // neutron
    particleType = neutron;
    stable       = true;
    exiting      = true;
    break;
  case 321:  // K
    particleType = kaon;
    stable       = false;
    exiting      = false;
    break;
  case 130:  // K_long
    particleType = k0;
    stable       = false;
    exiting      = false;
    break;
  case 310:  // K_short
    particleType = k0;
    stable       = false;
    exiting      = false;
    break;
  default:  // treat mesons as pions
    particleType = charge != 0. ? pion : pi0;
    stable       = false;
    exiting      = false;
  }

  // and all baryons as proton hypo
  if (pdgCode > 999 && pdgCode != 2112) {
    particleType = charge != 0. ? proton : neutron;
    stable       = false;
    exiting      = false;
  }

  // ignore SUSY particles for now
  if (pdgCode > 1000000) {
    particleType = nonInteracting;
    stable       = true;
    exiting      = true;
  }

  return particleType;
}

/**
 Static method : convert to pdg from ParticleType
static int
convertToPdg(ParticleType particleHypo, double charge, bool dist)
{
  int pdg = 0;

  switch (particleHypo) {
  // the electron case
  case electron: {
    pdg = 11;
    pdg *= charge > 0. ? -1 : 1;
  }
    return pdg;
  // the muon case
  case muon: {
    pdg = 13;
    pdg *= charge > 0. ? -1 : 1;
  }
    return pdg;
  // the kaon case
  case kaon: {
    pdg = 321;
    pdg *= charge > 0. ? 1 : -1;
  }
    return pdg;
  // the proton case
  case proton: {
    pdg = 2212;
    pdg *= charge > 0. ? 1 : -1;
    if (charge * charge < 0.0001) pdg = dist ? 2112 : -2112;
  }
    return pdg;
  // the photon case
  case photon: {
    pdg = 22;
  }
    return pdg;
  // the neutron case
  case neutron: {
    pdg = 2112;
  }
    return pdg;
  // the neutral pion case
  case pi0: {
    pdg = 111;
  }
    return pdg;
  // the neutral kaon case
  case k0: {
    pdg = dist ? 130 : 310;
  }
    return pdg;
  // the pion case - is the default
  default: {
    pdg = 211;
    pdg *= charge > 0. ? 1 : -1;
    if (charge * charge < 0.0001) pdg = 111;
  };
  }
  return pdg;
}
*/
/** @class ParticleProperties

    very simplistic class for particle properties,
    in order to be used in fast track simulation

*/
class ParticleProperties
{
public:
  /** constructor */
  ParticleProperties(const Vector3D& vertex,
                     const Vector3D& momentum,
                     double          mass    = 0.,
                     double          charge  = 0.,
                     pdg_type        pID     = 0.,
                     barcode_type    barcode = 0)
    : m_vertex(vertex)
    , m_momentum(momentum)
    , m_mass(mass)
    , m_charge(charge)
    , m_pdgID(pID)
    , m_particleType(pion)
    , m_barcode(barcode)
  {
    bool exiting, stable;
    m_particleType = convertToParticleType(pID, exiting, stable, charge);
  }

  /// constructor
  ParticleProperties(const ParticleProperties& pProperties)
    : m_vertex(pProperties.m_vertex)
    , m_momentum(pProperties.m_momentum)
    , m_mass(pProperties.m_mass)
    , m_charge(pProperties.m_charge)
    , m_pdgID(pProperties.m_pdgID)
    , m_particleType(pProperties.m_particleType)
    , m_barcode(pProperties.m_barcode)
  {
  }

  /// destructor
  ~ParticleProperties() {}

  /** assignment operator */
  ParticleProperties&
  operator=(const ParticleProperties& pProperties)
  {
    if (this != &pProperties) {
      m_vertex       = pProperties.m_vertex;
      m_momentum     = pProperties.m_momentum;
      m_mass         = pProperties.m_mass;
      m_charge       = pProperties.m_charge;
      m_pdgID        = pProperties.m_pdgID;
      m_particleType = pProperties.m_particleType;
      m_barcode      = pProperties.m_barcode;
    }
    return (*this);
  }

  /// shift the particle
  void
  shift(const Vector3D& shift)
  {
    m_vertex += shift;
  }

  /// assign the barcode
  void
  assign(barcode_type barcode)
  {
    m_barcode = barcode;
  }

  /// return the momentum
  const Vector3D&
  vertex() const
  {
    return m_vertex;
  }

  /// return the momentum
  const Vector3D&
  momentum() const
  {
    return m_momentum;
  }

  /// return the mass
  double
  mass() const
  {
    return m_mass;
  }

  /// return the charge
  double
  charge() const
  {
    return m_charge;
  }

  /// return the particle ID
  pdg_type
  pdgID() const
  {
    return m_pdgID;
  }

  /// return the particle type
  ParticleType
  particleType() const
  {
    return m_particleType;
  }

  /// return the particle barcode
  barcode_type
  barcode() const
  {
    return m_barcode;
  }

private:
  Vector3D     m_vertex;
  Vector3D     m_momentum;
  double       m_mass;
  double       m_charge;
  pdg_type     m_pdgID;
  ParticleType m_particleType;
  barcode_type m_barcode;
};

/// @class ProcessVertex
///
/// process vertex class for the fast track simulation
///

class ProcessVertex
{
public:
  /// constructor
  ProcessVertex(const Vector3D&                        pVertex,
                double                                 pTime,
                process_type                           pType,
                const std::vector<ParticleProperties>& pOutgoing)
    : m_vertex(pVertex), m_time(pTime), m_type(pType), m_outgoing(pOutgoing)
  {
  }

  /// destructor
  ~ProcessVertex() {}

  /// Add a particle
  void
  addOutgoing(const ParticleProperties& pProperties);

  /// Return the vertex position
  const Vector3D&
  vertex() const
  {
    return m_vertex;
  }

  /// Return the time of production
  double
  interactionTime() const
  {
    return m_time;
  }

  /// Return the type of production
  process_type
  processType() const
  {
    return m_type;
  }

  /// Return the outgoing properties
  const std::vector<ParticleProperties>&
  outgoingParticles() const
  {
    return m_outgoing;
  }

private:
  Vector3D                        m_vertex;
  double                          m_time;
  process_type                    m_type;
  std::vector<ParticleProperties> m_outgoing;
};

inline void
ProcessVertex::addOutgoing(const ParticleProperties& pProperties)
{
  m_outgoing.push_back(pProperties);
}
}

#endif  // ACTS_EVENTDATA_PARTICLEDEFINITIONS_H
