// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::Experimental {

/// GBTS region of interest, as the eta and z bounds the seeder cuts on.
///
/// The eta bounds are only kept as the dz/dr gradients they imply. There is no
/// phi range: the seeder searches the full azimuth.
class GbtsRoiDescriptor {
 public:
  /// @param etaMin eta at rear of RoI
  /// @param etaMax eta at front of RoI
  /// @param zMin z at rear of RoI
  /// @param zMax z at front of RoI
  GbtsRoiDescriptor(double etaMin, double etaMax, double zMin = 0,
                    double zMax = 0);

  /// Get z at the most backward end of the RoI
  /// @return Z at rear
  double zMin() const { return m_zMin; }
  /// Get z at the most forward end of the RoI
  /// @return Z at front
  double zMax() const { return m_zMax; }

  /// Get dz/dr at the rear of the RoI
  /// @return Gradient dzdr at rear
  double dzdrMin() const { return m_dzdrMin; }
  /// Get dz/dr at the front of the RoI
  /// @return Gradient dzdr at front
  double dzdrMax() const { return m_dzdrMax; }

 private:
  float m_zMin{};  //!< z position at most negative position along the beamline
  float m_zMax{};  //!< z position at most positive position along the beamline

  float m_dzdrMin{};  //!<  dz/dr at the rear of the RoI
  float m_dzdrMax{};  //!<  dz/dr at the front of the RoI
};

}  // namespace Acts::Experimental
