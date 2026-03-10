// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::Experimental {

/// GBTS region-of-interest descriptor with eta/phi/zed bounds.
class GbtsRoiDescriptor {
 public:
  /// @param eta eta of RoI
  /// @param etaMin eta at rear  of RoI
  /// @param etaMax eta at front of RoI
  /// @param phi phi of RoI
  /// @param phiMin minimum phi of RoI
  /// @param phiMax maximum phi of RoI
  /// @param z z of RoI
  /// @param zMin z at rear  of RoI
  /// @param zMax z at front of RoI
  GbtsRoiDescriptor(double eta, double etaMin, double etaMax, double phi,
                    double phiMin, double phiMax, double z = 0, double zMin = 0,
                    double zMax = 0);

  // Methods to retrieve data members

  /// Get phi coordinate of RoI center
  /// @return Phi coordinate
  double phi() const { return m_phi; }
  /// Get eta coordinate of RoI center
  /// @return Eta coordinate
  double eta() const { return m_eta; }
  /// Get z coordinate of RoI center
  /// @return Z coordinate
  double z() const { return m_z; }

  /// these quantities probably don't need to be used any more
  /// - they are implemented here only because we had them in
  ///   the original legacy interface

  /// Get z at the most backward end of the RoI
  /// @return Z at rear
  double zMin() const { return m_zMin; }
  /// Get z at the most forward end of the RoI
  /// @return Z at front
  double zMax() const { return m_zMax; }

  /// Gets eta at zMin
  /// @return Eta at rear
  double etaMin() const { return m_etaMin; }
  /// Gets eta at zMax
  /// @return Eta at front
  double etaMax() const { return m_etaMax; }

  /// Gets phiMinus
  /// @return Minimum phi
  double phiMin() const { return m_phiMin; }
  /// Gets phiPlus
  /// @return Maximum phi
  double phiMax() const { return m_phiMax; }

  // return the gradients
  /// Get dz/dr at the rear of the RoI
  /// @return Gradient dzdr at rear
  double dzdrMin() const { return m_dzdrMin; }
  /// Get dz/dr at the front of the RoI
  /// @return Gradient dzdr at front
  double dzdrMax() const { return m_dzdrMax; }

  /// Get dr/dz at the rear of the RoI
  /// @return Gradient drdz at rear
  double drdzMin() const { return m_drdzMin; }
  /// Get dr/dz at the front of the RoI
  /// @return Gradient drdz at front
  double drdzMax() const { return m_drdzMax; }

  /// Get z at the most backward end of the RoI at outer radius
  /// @return Z at rear outer radius
  double zOuterMin() const { return m_zOuterMin; }
  /// Get z at the most forward end of the RoI at outer radius
  /// @return Z at front outer radius
  double zOuterMax() const { return m_zOuterMax; }

 private:
  float m_phi{};  //!< phi of RoI center
  float m_eta{};  //!< eta of RoI center
  float m_z{};    //!< z of RoI center

  float m_phiMin{};  //!< most negative RoI in azimuthal
  float m_phiMax{};  //!< most positive RoI in azimuthal

  float m_etaMin{};  //!< eta of RoI at zMax
  float m_etaMax{};  //!< eta of RoI at zMin

  float m_zMin{};  //!< z position at most negative position along the beamline
  float m_zMax{};  //!< z position at most positive position along the beamline

  float m_dzdrMin{};  //!<  dz/dr at the rear of the RoI
  float m_dzdrMax{};  //!<  dz/dr at the front of the RoI

  float m_drdzMin{};  //!<  dr/dz at the rear of the RoI
  float m_drdzMax{};  //!<  dr/dz at the front of the RoI

  float m_zOuterMin{};  //!< z at rear of RoI at the outer radius ( = 1100 mm)
  float m_zOuterMax{};  //!< z at front of RoI at the outer radius ( = 1100 mm)
};

}  // namespace Acts::Experimental
