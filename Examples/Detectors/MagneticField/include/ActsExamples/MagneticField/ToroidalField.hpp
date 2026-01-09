// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <vector>

namespace Acts {

class ToroidalField final : public MagneticFieldProvider {
 public:
  struct BarrelConfig {
    float R_in = 4.9f;    // [m]
    float R_out = 10.0f;  // [m]
    float c = 25.3f;      // [m]
    float b = 0.16f;      // [m]
    float I = 20500.0f;   // [A]
    int Nturns = 120;     // [turns]
  };

  struct ECTConfig {
    float R_in = 1.65f * 0.5f;   // [m] ≈ 0.825
    float R_out = 10.7f * 0.5f;  // [m] ≈ 5.35
    float c = 5.0f;              // [m]
    float b = 0.12f;             // [m]
    float I = 20500.0f;          // [A]
    int Nturns = 116;            // [turns]
    float gap = 0.5f;            // [m]
  };

  struct LayoutConfig {
    float theta0_deg = 22.5f;
    float thetaStep_deg = 45.0f;
    int nCoils = 8;

    int nArc = 200;
    int nStraight = 160;
    bool closeLoop = true;

    float eps = 1e-18f;
  };

  struct Config {
    BarrelConfig barrel;
    ECTConfig ect;
    LayoutConfig layout;

    // Per-coil current senses (applied to dl):
    // - barrelSigns size must be nCoils
    // - ectSigns    size must be 2*nCoils (first +z endcap [0..nCoils-1], then
    // -z [nCoils..2*nCoils-1])
    std::vector<int> barrelSigns;  // default filled in ctor if empty
    std::vector<int> ectSigns;     // default filled in ctor if empty
  };

  struct Cache {
    explicit Cache(const MagneticFieldContext& /*ctx*/) {}
  };

  ToroidalField();                     // default config
  explicit ToroidalField(Config cfg);  // with config

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override;

  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override;

  bool isInside(const Vector3& /*position*/) const { return true; }

  const Config& config() const { return m_cfg; }

 private:
  static float deg2rad(float deg);

  // Helpers (declared here, defined in .cpp)
  static std::vector<std::array<float, 2>> ectRacetrackRadial(
      float Lrho, float Lz, int nArc, int nStraight, bool close);

  static std::vector<std::array<float, 2>> racetrackRZ(float a, float b,
                                                       float Lz, int nArc,
                                                       int nStraight,
                                                       bool close);

  static void buildSegsMidsRZ(const std::vector<std::array<float, 2>>& rz,
                              std::vector<std::array<float, 2>>& d_rz,
                              std::vector<std::array<float, 2>>& m_rz);

  static void mapRingToXYZ(float l,
                           const std::vector<std::array<float, 2>>& m_rz,
                           const std::vector<std::array<float, 2>>& d_rz,
                           float phi, int sign, float zShift,
                           std::vector<std::array<float, 3>>& mids_out,
                           std::vector<std::array<float, 3>>& segs_out);

  void buildGeometry();

 private:
  Config m_cfg;

  // Precomputed geometry (float storage; double accumulation)
  std::vector<std::array<float, 3>> m_segs_barrel;
  std::vector<std::array<float, 3>> m_mids_barrel;

  std::vector<std::array<float, 3>> m_segs_ect;  // both endcaps combined
  std::vector<std::array<float, 3>> m_mids_ect;

  void accumulateBarrelField(double X, double Y, double Z, double eps,
                             double pref, double& bx, double& by,
                             double& bz) const;

  void accumulateEndcapField(double X, double Y, double Z, double eps,
                             double pref, double& bx, double& by,
                             double& bz) const;
};

}  // namespace Acts
