// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <vector>

namespace Acts {

/// Toroid magnetic field implementation
///
/// This class implements a toroid magnetic field configuration similar to
/// those used in ATLAS and other detector systems. It uses Biot-Savart
/// calculations with discrete current-carrying wire segments to compute
/// the field at any position.
class ToroidField final : public MagneticFieldProvider {
 public:
  /// Configuration for barrel toroid coils
  struct BarrelConfig {
    double R_in = 4.9 * UnitConstants::m;    ///< Inner radius
    double R_out = 10.0 * UnitConstants::m;  ///< Outer radius
    double c = 25.3 * UnitConstants::m;      ///< Coil length along z
    double b = 0.16 * UnitConstants::m;      ///< Coil width (radial)
    double I = 20500.0;                      ///< Current [A]
    int Nturns = 120;                        ///< Number of turns
  };

  /// Configuration for end-cap toroid (ECT) coils
  struct EctConfig {
    double R_in = 1.65 * 0.5 * UnitConstants::m;   ///< Inner radius ≈ 0.825 m
    double R_out = 10.7 * 0.5 * UnitConstants::m;  ///< Outer radius ≈ 5.35 m
    double c = 5.0 * UnitConstants::m;             ///< Coil length along z
    double b = 0.12 * UnitConstants::m;            ///< Coil width (radial)
    double I = 20500.0;                            ///< Current [A]
    int Nturns = 116;                              ///< Number of turns
    double gap = 0.5 * UnitConstants::m;  ///< Gap between barrel and endcap
  };

  /// Configuration for coil layout and discretization
  struct LayoutConfig {
    double theta0 = 22.5 * UnitConstants::degree;  ///< Initial azimuthal angle
    double thetaStep =
        45.0 * UnitConstants::degree;  ///< Angular spacing between coils
    int nCoils = 8;                    ///< Number of coils

    int nArc = 200;         ///< Number of segments in arc portions
    int nStraight = 160;    ///< Number of segments in straight portions
    bool closeLoop = true;  ///< Whether to close the coil loop

    double eps = 1e-18;  ///< Small epsilon for numerical stability
  };

  /// Full configuration for the toroid field
  struct Config {
    BarrelConfig barrel;  ///< Barrel coil configuration
    EctConfig ect;        ///< End-cap toroid configuration
    LayoutConfig layout;  ///< Layout and discretization parameters

    /// Per-coil current senses (applied to dl)
    /// Size must be nCoils; default filled in ctor if empty
    std::vector<int> barrelSigns;
    /// Per-coil current senses for endcaps
    /// Size must be 2*nCoils (first +z endcap [0..nCoils-1], then -z
    /// [nCoils..2*nCoils-1]); default filled in ctor if empty
    std::vector<int> ectSigns;
  };

  /// Cache for magnetic field provider
  struct Cache {};

  /// Construct with default configuration
  ToroidField();

  /// Construct with custom configuration
  /// @param cfg Configuration parameters
  explicit ToroidField(Config cfg);

  /// @copydoc MagneticFieldProvider::makeCache
  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override;

  /// @copydoc MagneticFieldProvider::getField
  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override;

  /// Get the configuration
  /// @return Configuration struct reference
  const Config& config() const { return m_cfg; }

 private:
  // Helpers (declared here, defined in .cpp)
  static std::vector<std::array<float, 2>> ectRacetrackRadial(
      float Lrho, float Lz, int nArc, int nStraight, bool close);

  static std::vector<std::array<double, 2>> racetrackRZ(double a, double b,
                                                        double Lz, int nArc,
                                                        int nStraight,
                                                        bool close);

  static void buildSegsMidsRZ(const std::vector<std::array<double, 2>>& rz,
                              std::vector<std::array<double, 2>>& d_rz,
                              std::vector<std::array<double, 2>>& m_rz);

  static void mapRingToXYZ(double l,
                           const std::vector<std::array<double, 2>>& m_rz,
                           const std::vector<std::array<double, 2>>& d_rz,
                           double phi, int sign, double zShift,
                           std::vector<std::array<double, 3>>& mids_out,
                           std::vector<std::array<double, 3>>& segs_out);

  void buildGeometry();

  Config m_cfg;

  // Precomputed geometry (float storage; double accumulation)
  std::vector<std::array<double, 3>> m_segs_barrel;
  std::vector<std::array<double, 3>> m_mids_barrel;

  std::vector<std::array<double, 3>> m_segs_ect;  // both endcaps combined
  std::vector<std::array<double, 3>> m_mids_ect;

  void accumulateBarrelField(double X, double Y, double Z, double eps,
                             double pref, double& bx, double& by,
                             double& bz) const;

  void accumulateEndcapField(double X, double Y, double Z, double eps,
                             double pref, double& bx, double& by,
                             double& bz) const;
};

}  // namespace Acts
