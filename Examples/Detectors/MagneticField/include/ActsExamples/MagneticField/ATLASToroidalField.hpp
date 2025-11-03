// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// ATLASToroidalField.hpp  (header-only)
// Full ATLAS toroidal field (barrel + two endcaps) with per-coil current
// senses.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace Acts {

class ATLASToroidalField final : public MagneticFieldProvider {
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
    float gap = 0.5f;            // [m] small mechanical gap
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

  ATLASToroidalField() : ATLASToroidalField(Config{}) {}

  explicit ATLASToroidalField(Config cfg) : m_cfg(std::move(cfg)) {
    // Fill default signs if not provided
    const int nC = m_cfg.layout.nCoils;
    if (m_cfg.barrelSigns.empty()) {
      m_cfg.barrelSigns.resize(nC);
      for (int k = 0; k < nC; ++k)
        m_cfg.barrelSigns[k] = (k % 2 == 0) ? +1 : -1;
    }
    if (m_cfg.ectSigns.empty()) {
      // Updated default to match Python:
      // coil_signs_ect_front = coil_signs_ect_back = [-1, +1, -1, +1, -1, +1,
      // -1, +1]
      m_cfg.ectSigns.resize(2 * nC);
      for (int k = 0; k < nC; ++k) {
        const int alt = (k % 2 == 0) ? -1 : +1;
        m_cfg.ectSigns[k] = alt;       // +z endcap
        m_cfg.ectSigns[nC + k] = alt;  // -z endcap
      }
    }
    // Sanity on sizes
    if (static_cast<int>(m_cfg.barrelSigns.size()) != nC)
      throw std::invalid_argument(
          "ATLASToroidalField: barrelSigns size must equal nCoils");
    if (static_cast<int>(m_cfg.ectSigns.size()) != 2 * nC)
      throw std::invalid_argument(
          "ATLASToroidalField: ectSigns size must equal 2*nCoils");

    buildGeometry();
  }

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override {
    (void)cache;

    const float X = static_cast<float>(position.x());
    const float Y = static_cast<float>(position.y());
    const float Z = static_cast<float>(position.z());

    double bx = 0.0, by = 0.0, bz = 0.0;

    constexpr double mu0 = 4e-7 * std::numbers::pi;  // [T·m/A]
    const double prefBarrel = (mu0 * static_cast<double>(m_cfg.barrel.Nturns) *
                               static_cast<double>(m_cfg.barrel.I)) /
                              (4.0 * std::numbers::pi);
    const double prefECT = (mu0 * static_cast<double>(m_cfg.ect.Nturns) *
                            static_cast<double>(m_cfg.ect.I)) /
                           (4.0 * std::numbers::pi);
    const float eps = m_cfg.layout.eps;

    // Barrel
    for (std::size_t i = 0; i < m_mids_barrel.size(); ++i) {
      const auto& mid = m_mids_barrel[i];
      const auto& dl = m_segs_barrel[i];
      const float rx = X - mid[0], ry = Y - mid[1], rz = Z - mid[2];
      const double r2 = static_cast<double>(rx) * rx +
                        static_cast<double>(ry) * ry +
                        static_cast<double>(rz) * rz + static_cast<double>(eps);
      const double invr = 1.0 / std::sqrt(r2), invr3 = invr / r2;
      const double cx =
          static_cast<double>(dl[1]) * rz - static_cast<double>(dl[2]) * ry;
      const double cy =
          static_cast<double>(dl[2]) * rx - static_cast<double>(dl[0]) * rz;
      const double cz =
          static_cast<double>(dl[0]) * ry - static_cast<double>(dl[1]) * rx;
      bx += prefBarrel * cx * invr3;
      by += prefBarrel * cy * invr3;
      bz += prefBarrel * cz * invr3;
    }

    // Endcaps (both)
    for (std::size_t i = 0; i < m_mids_ect.size(); ++i) {
      const auto& mid = m_mids_ect[i];
      const auto& dl = m_segs_ect[i];
      const float rx = X - mid[0], ry = Y - mid[1], rz = Z - mid[2];
      const double r2 = static_cast<double>(rx) * rx +
                        static_cast<double>(ry) * ry +
                        static_cast<double>(rz) * rz + static_cast<double>(eps);
      const double invr = 1.0 / std::sqrt(r2), invr3 = invr / r2;
      const double cx =
          static_cast<double>(dl[1]) * rz - static_cast<double>(dl[2]) * ry;
      const double cy =
          static_cast<double>(dl[2]) * rx - static_cast<double>(dl[0]) * rz;
      const double cz =
          static_cast<double>(dl[0]) * ry - static_cast<double>(dl[1]) * rx;
      bx += prefECT * cx * invr3;
      by += prefECT * cy * invr3;
      bz += prefECT * cz * invr3;
    }

    return Result<Vector3>::success(Vector3(bx, by, bz));
  }

  bool isInside(const Vector3& /*position*/) const { return true; }
  const Config& config() const { return m_cfg; }

 private:
  static float deg2rad(float deg) {
    return deg * static_cast<float>(std::numbers::pi / 180.0);
  }

  // End-cap racetrack with LONG straights along z:
  //  - straights at ρ = ±(Lrho/2),    z ∈ [-Lz/2, +Lz/2]
  //  - arcs      at z = ±(Lz/2) (constant z), ρ sweeps via rr*sin(θ)
  // Mirrors the Python ect_racetrack_radial() exactly.
  static std::vector<std::array<float, 2>> ectRacetrackRadial(
      float Lrho, float Lz, int nArc, int nStraight, bool close) {
    const float rr = 0.5f * Lrho;  // radial half-span
    const float rz = 0.5f * Lz;    // axial half-length
    std::vector<std::array<float, 2>> pts;
    pts.reserve(
        static_cast<std::size_t>(2 * nArc + 2 * nStraight + (close ? 1 : 0)));

    // Straight at ρ = +rr : z from +rz -> -rz (exclude endpoint)
    for (int i = 0; i < nStraight; ++i) {
      const float t = static_cast<float>(i) / static_cast<float>(nStraight);
      const float z = (+rz) * (1.0f - t) + (-rz) * t;
      pts.push_back({+rr, z});
    }
    // Inner arc at z = -rz : θ: +π/2 → -π/2 (exclude endpoint) ; ρ = rr*sin(θ)
    for (int i = 0; i < nArc; ++i) {
      const float th = static_cast<float>(std::numbers::pi / 2) +
                       (static_cast<float>(i) / nArc) *
                           static_cast<float>(-std::numbers::pi);
      const float rho = rr * std::sin(th);
      const float z = -rz;
      pts.push_back({rho, z});
    }
    // Straight at ρ = -rr : z from -rz -> +rz (exclude endpoint)
    for (int i = 0; i < nStraight; ++i) {
      const float t = static_cast<float>(i) / static_cast<float>(nStraight);
      const float z = (-rz) * (1.0f - t) + (+rz) * t;
      pts.push_back({-rr, z});
    }
    // Outer arc at z = +rz : θ: -π/2 → +π/2 (close if requested)
    const int nLast = nArc;
    for (int i = 0; i < nLast; ++i) {
      const float th = static_cast<float>(-std::numbers::pi / 2) +
                       (static_cast<float>(i) / nLast) *
                           static_cast<float>(std::numbers::pi);
      const float rho = rr * std::sin(th);
      const float z = +rz;
      pts.push_back({rho, z});
    }
    if (close) {
      if (pts.front()[0] != pts.back()[0] || pts.front()[1] != pts.back()[1]) {
        pts.push_back(pts.front());
      }
    }
    return pts;
  }

  // Build a single racetrack in (ρ,z) as polyline
  static std::vector<std::array<float, 2>> racetrackRZ(float a, float b,
                                                       float Lz, int nArc,
                                                       int nStraight,
                                                       bool close) {
    const float r = 0.5f * b;
    std::vector<std::array<float, 2>> pts;
    pts.reserve(
        static_cast<std::size_t>(2 * nArc + 2 * nStraight + (close ? 1 : 0)));

    // Arc at x=+Lz/2 (x≡z, y≡ρ)
    for (int i = 0; i < nArc; ++i) {
      const float th =
          static_cast<float>(std::numbers::pi) / 2 +
          (static_cast<float>(i) / nArc) * static_cast<float>(std::numbers::pi);
      const float x = (+Lz * 0.5f) + r * std::cos(th);
      const float y = r * std::sin(th);
      pts.push_back({y, x});
    }
    // Straight at y=-r from +Lz/2 -> -Lz/2
    for (int i = 0; i < nStraight; ++i) {
      const float t = static_cast<float>(i) / static_cast<float>(nStraight);
      const float x = (+Lz * 0.5f) * (1.0f - t) + (-Lz * 0.5f) * t;
      pts.push_back({-r, x});
    }
    // Arc at x=-Lz/2
    for (int i = 0; i < nArc; ++i) {
      const float th = static_cast<float>(3.0 * std::numbers::pi / 2) +
                       (static_cast<float>(i) / nArc) *
                           static_cast<float>(-std::numbers::pi);
      const float x = (-Lz * 0.5f) + r * std::cos(th);
      const float y = r * std::sin(th);
      pts.push_back({y, x});
    }
    // Straight at y=+r from -Lz/2 -> +Lz/2
    const int nStraightLast = close ? nStraight : (nStraight + 1);
    for (int i = 0; i < nStraightLast; ++i) {
      const float t =
          static_cast<float>(i) / static_cast<float>(nStraightLast - 1);
      const float x = (-Lz * 0.5f) * (1.0f - t) + (+Lz * 0.5f) * t;
      pts.push_back({+r, x});
    }

    // Shift to place straight legs at ρ=±a/2 (keeping arc radius r)
    const float delta = 0.5f * a - r;
    for (auto& p : pts) {
      const float sgn = (p[0] >= 0.0f) ? 1.0f : -1.0f;
      p[0] = p[0] + sgn * delta;
    }

    if (close) {
      if (pts.front()[0] != pts.back()[0] || pts.front()[1] != pts.back()[1])
        pts.push_back(pts.front());
    }
    return pts;
  }

  static void buildSegsMidsRZ(const std::vector<std::array<float, 2>>& rz,
                              std::vector<std::array<float, 2>>& d_rz,
                              std::vector<std::array<float, 2>>& m_rz) {
    d_rz.clear();
    m_rz.clear();
    d_rz.reserve(rz.size() - 1);
    m_rz.reserve(rz.size() - 1);
    for (std::size_t i = 0; i + 1 < rz.size(); ++i) {
      const float dr = rz[i + 1][0] - rz[i][0];
      const float dz = rz[i + 1][1] - rz[i][1];
      d_rz.push_back({dr, dz});
      m_rz.push_back(
          {0.5f * (rz[i + 1][0] + rz[i][0]), 0.5f * (rz[i + 1][1] + rz[i][1])});
    }
  }

  // Map one ring with sign applied to dl. zShift for endcaps.
  static void mapRingToXYZ(float l,
                           const std::vector<std::array<float, 2>>& m_rz,
                           const std::vector<std::array<float, 2>>& d_rz,
                           float phi, int sign, float zShift,
                           std::vector<std::array<float, 3>>& mids_out,
                           std::vector<std::array<float, 3>>& segs_out) {
    const float ct = std::cos(phi), st = std::sin(phi);
    const float s = (sign >= 0) ? 1.0f : -1.0f;

    for (const auto& rm : m_rz) {
      const float rho = rm[0];
      const float zz = rm[1] + zShift;
      const float rxy = l + rho;
      mids_out.push_back({rxy * ct, rxy * st, zz});
    }
    for (const auto& dlrz : d_rz) {
      const float dr = dlrz[0];
      const float dz = dlrz[1];
      segs_out.push_back({s * (dr * ct), s * (dr * st), s * dz});
    }
  }

  void buildGeometry() {
    // ---- Barrel base curve ----
    const float rEndB = 0.5f * m_cfg.barrel.b;
    const float lB = 0.5f * (m_cfg.barrel.R_in + m_cfg.barrel.R_out);
    const float aB = (m_cfg.barrel.R_out - m_cfg.barrel.R_in) - 2.0f * rEndB;

    const auto rz_barrel =
        racetrackRZ(aB, m_cfg.barrel.b, m_cfg.barrel.c, m_cfg.layout.nArc,
                    m_cfg.layout.nStraight, m_cfg.layout.closeLoop);
    std::vector<std::array<float, 2>> d_rzB, m_rzB;
    buildSegsMidsRZ(rz_barrel, d_rzB, m_rzB);

    // ---- ECT base curve (use exact ect_racetrack_radial shape) ----
    // Lrho_ECT = R_out_ECT - R_in_ECT ; Lz = c_ECT
    const float lE = 0.5f * (m_cfg.ect.R_in + m_cfg.ect.R_out);
    const float LrhoE = (m_cfg.ect.R_out - m_cfg.ect.R_in);
    // straights at ρ=±Lrho/2, arcs at z=±c_ECT/2 (constant-z arcs)
    const auto rz_ect =
        ectRacetrackRadial(/*Lrho=*/LrhoE,
                           /*Lz=*/m_cfg.ect.c, m_cfg.layout.nArc,
                           m_cfg.layout.nStraight, m_cfg.layout.closeLoop);
    std::vector<std::array<float, 2>> d_rzE, m_rzE;
    buildSegsMidsRZ(rz_ect, d_rzE, m_rzE);

    // ---- Angles ----
    const int nC = m_cfg.layout.nCoils;
    const float th0 = deg2rad(m_cfg.layout.theta0_deg);
    const float dth = deg2rad(m_cfg.layout.thetaStep_deg);

    // ---- Reserve & fill ----
    m_mids_barrel.clear();
    m_segs_barrel.clear();
    m_mids_ect.clear();
    m_segs_ect.clear();

    m_mids_barrel.reserve(static_cast<std::size_t>(nC) * m_rzB.size());
    m_segs_barrel.reserve(static_cast<std::size_t>(nC) * d_rzB.size());
    m_mids_ect.reserve(static_cast<std::size_t>(2 * nC) * m_rzE.size());
    m_segs_ect.reserve(static_cast<std::size_t>(2 * nC) * d_rzE.size());

    // Barrel rings with alternating signs
    for (int k = 0; k < nC; ++k) {
      const float phi = th0 + k * dth;
      const int sign = m_cfg.barrelSigns[k];
      mapRingToXYZ(lB, m_rzB, d_rzB, phi, sign, /*zShift=*/0.0f, m_mids_barrel,
                   m_segs_barrel);
    }

    // Endcap centers: z = ±(c/2 + gap + c_ECT/2)
    const float zECT =
        0.5f * m_cfg.barrel.c + m_cfg.ect.gap + 0.5f * m_cfg.ect.c;

    // +z endcap (indices 0..nC-1 in ectSigns)
    for (int k = 0; k < nC; ++k) {
      const float phi = th0 + k * dth;
      const int sign = m_cfg.ectSigns[k];
      mapRingToXYZ(lE, m_rzE, d_rzE, phi, sign, /*zShift=*/+zECT, m_mids_ect,
                   m_segs_ect);
    }
    // -z endcap (indices nC..2*nC-1 in ectSigns)
    for (int k = 0; k < nC; ++k) {
      const float phi = th0 + k * dth;
      const int sign = m_cfg.ectSigns[nC + k];
      mapRingToXYZ(lE, m_rzE, d_rzE, phi, sign, /*zShift=*/-zECT, m_mids_ect,
                   m_segs_ect);
    }
  }

 private:
  Config m_cfg;

  // Precomputed geometry (float storage; double accumulation)
  std::vector<std::array<float, 3>> m_segs_barrel;
  std::vector<std::array<float, 3>> m_mids_barrel;

  std::vector<std::array<float, 3>> m_segs_ect;  // both endcaps combined
  std::vector<std::array<float, 3>> m_mids_ect;
};

}  // namespace Acts
