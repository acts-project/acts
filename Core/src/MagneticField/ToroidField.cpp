// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN
//
// MPL 2.0: https://mozilla.org/MPL/2.0/

#include "Acts/MagneticField/ToroidField.hpp"

#include "Acts/Definitions/Units.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace Acts {

//-------------------------
// Ctors
//-------------------------

ToroidField::ToroidField() : ToroidField(Config{}) {}

ToroidField::ToroidField(Config cfg) : m_cfg(std::move(cfg)) {
  // Fill default signs if not provided
  const int nC = m_cfg.layout.nCoils;
  if (m_cfg.barrelSigns.empty()) {
    m_cfg.barrelSigns.assign(nC, +1);
  }
  if (m_cfg.ectSigns.empty()) {
    m_cfg.ectSigns.assign(2 * nC, +1);
  }
  // Sanity on sizes
  if (static_cast<int>(m_cfg.barrelSigns.size()) != nC) {
    throw std::invalid_argument(
        "ToroidField: barrelSigns size must equal nCoils");
  }
  if (static_cast<int>(m_cfg.ectSigns.size()) != 2 * nC) {
    throw std::invalid_argument(
        "ToroidField: ectSigns size must equal 2*nCoils");
  }

  buildGeometry();
}

//-------------------------
// MagneticFieldProvider
//-------------------------

MagneticFieldProvider::Cache ToroidField::makeCache(
    const MagneticFieldContext& mctx) const {
  return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
}

Result<Vector3> ToroidField::getField(
    const Vector3& position, MagneticFieldProvider::Cache& cache) const {
  (void)cache;

  const double X = position.x();
  const double Y = position.y();
  const double Z = position.z();

  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;

  constexpr double mu0 = 4e-7 * std::numbers::pi;  // [T·m/A]
  const double prefBarrel =
      (mu0 * static_cast<double>(m_cfg.barrel.Nturns) * m_cfg.barrel.I) /
      (4.0 * std::numbers::pi);
  const double prefECT =
      (mu0 * static_cast<double>(m_cfg.ect.Nturns) * m_cfg.ect.I) /
      (4.0 * std::numbers::pi);
  const double eps = m_cfg.layout.eps;

  // Split computations
  accumulateBarrelField(X, Y, Z, eps, prefBarrel, bx, by, bz);
  accumulateEndcapField(X, Y, Z, eps, prefECT, bx, by, bz);

  return Result<Vector3>::success(Vector3(bx, by, bz));
}

//-------------------------
// Private helpers
//-------------------------

void ToroidField::accumulateBarrelField(double X, double Y, double Z,
                                        double eps, double pref, double& bx,
                                        double& by, double& bz) const {
  for (std::size_t i = 0; i < m_mids_barrel.size(); ++i) {
    const auto& mid = m_mids_barrel[i];
    const auto& dl = m_segs_barrel[i];

    const double rx = X - mid[0];
    const double ry = Y - mid[1];
    const double rz = Z - mid[2];

    const double r2 = rx * rx + ry * ry + rz * rz + eps;
    const double invr = 1.0 / std::sqrt(r2);
    const double invr3 = invr / r2;

    const double cx = dl[1] * rz - dl[2] * ry;
    const double cy = dl[2] * rx - dl[0] * rz;
    const double cz = dl[0] * ry - dl[1] * rx;

    bx += pref * cx * invr3;
    by += pref * cy * invr3;
    bz += pref * cz * invr3;
  }
}

void ToroidField::accumulateEndcapField(double X, double Y, double Z,
                                        double eps, double pref, double& bx,
                                        double& by, double& bz) const {
  for (std::size_t i = 0; i < m_mids_ect.size(); ++i) {
    const auto& mid = m_mids_ect[i];
    const auto& dl = m_segs_ect[i];

    const double rx = X - mid[0];
    const double ry = Y - mid[1];
    const double rz = Z - mid[2];

    const double r2 = rx * rx + ry * ry + rz * rz + eps;
    const double invr = 1.0 / std::sqrt(r2);
    const double invr3 = invr / r2;

    const double cx = dl[1] * rz - dl[2] * ry;
    const double cy = dl[2] * rx - dl[0] * rz;
    const double cz = dl[0] * ry - dl[1] * rx;

    bx += pref * cx * invr3;
    by += pref * cy * invr3;
    bz += pref * cz * invr3;
  }
}

// End-cap racetrack with LONG straights along z (kept for parity with
// header-only version)
std::vector<std::array<float, 2>> ToroidField::ectRacetrackRadial(
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
                     (static_cast<float>(i) / static_cast<float>(nArc)) *
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
  for (int i = 0; i < nArc; ++i) {
    const float th = -static_cast<float>(std::numbers::pi) / 2 +
                     (static_cast<float>(i) / static_cast<float>(nArc)) *
                         static_cast<float>(std::numbers::pi);
    const float rho = rr * std::sin(th);
    const float z = +rz;
    pts.push_back({rho, z});
  }
  if (close &&
      (pts.front()[0] != pts.back()[0] || pts.front()[1] != pts.back()[1])) {
    pts.push_back(pts.front());
  }
  return pts;
}

// Rounded-rectangle racetrack in local (ρ, z)
std::vector<std::array<double, 2>> ToroidField::racetrackRZ(double a, double b,
                                                            double Lz, int nArc,
                                                            int nStraight,
                                                            bool close) {
  const double r = 0.5 * b;     // corner radius
  const double rz = 0.5 * Lz;   // axial half-length
  const double zTop = rz - r;   // top straight z
  const double zBot = -rz + r;  // bottom straight z

  std::vector<std::array<double, 2>> pts;
  pts.reserve(
      static_cast<std::size_t>(2 * nArc + 2 * nStraight + (close ? 1 : 0)));

  // +ρ straight (top → bottom), ρ = +a/2
  for (int i = 0; i < nStraight; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(nStraight);
    const double z = zTop * (1.0 - t) + zBot * t;
    pts.push_back({+0.5 * a, z});
  }
  // bottom-right quarter
  for (int i = 0; i < nArc; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(nArc);
    const double th = -std::numbers::pi / 2 - t * (std::numbers::pi / 2);
    const double cx = +0.5 * a - r;  // ρ-center
    const double cz = -rz + r;       // z-center
    const double rho = cx + r * std::cos(th);
    const double z = cz + r * std::sin(th);
    pts.push_back({rho, z});
  }
  // −ρ straight (bottom → top), ρ = −a/2
  for (int i = 0; i < nStraight; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(nStraight);
    const double z = zBot * (1.0 - t) + zTop * t;
    pts.push_back({-0.5 * a, z});
  }
  // top-left quarter
  for (int i = 0; i < nArc; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(nArc);
    const double th = +std::numbers::pi / 2 - t * (std::numbers::pi / 2);
    const double cx = -0.5 * a + r;  // ρ-center
    const double cz = +rz - r;       // z-center
    const double rho = cx + r * std::cos(th);
    const double z = cz + r * std::sin(th);
    pts.push_back({rho, z});
  }
  if (close &&
      (pts.front()[0] != pts.back()[0] || pts.front()[1] != pts.back()[1])) {
    pts.push_back(pts.front());
  }
  return pts;
}

void ToroidField::buildSegsMidsRZ(const std::vector<std::array<double, 2>>& rz,
                                  std::vector<std::array<double, 2>>& d_rz,
                                  std::vector<std::array<double, 2>>& m_rz) {
  d_rz.clear();
  m_rz.clear();
  d_rz.reserve(rz.size() - 1);
  m_rz.reserve(rz.size() - 1);
  for (std::size_t i = 0; i + 1 < rz.size(); ++i) {
    const double dr = rz[i + 1][0] - rz[i][0];
    const double dz = rz[i + 1][1] - rz[i][1];
    d_rz.push_back({dr, dz});
    m_rz.push_back(
        {0.5f * (rz[i + 1][0] + rz[i][0]), 0.5 * (rz[i + 1][1] + rz[i][1])});
  }
}

void ToroidField::mapRingToXYZ(double l,
                               const std::vector<std::array<double, 2>>& m_rz,
                               const std::vector<std::array<double, 2>>& d_rz,
                               double phi, int sign, double zShift,
                               std::vector<std::array<double, 3>>& mids_out,
                               std::vector<std::array<double, 3>>& segs_out) {
  const double ct = std::cos(phi);
  const double st = std::sin(phi);
  const double s = (sign >= 0) ? 1.0 : -1.0;

  for (const auto& rm : m_rz) {
    const double rho = rm[0];
    const double zz = rm[1] + zShift;
    const double rxy = l + rho;
    mids_out.push_back({rxy * ct, rxy * st, zz});
  }
  for (const auto& dlrz : d_rz) {
    const double dr = dlrz[0];
    const double dz = dlrz[1];
    segs_out.push_back({s * (dr * ct), s * (dr * st), s * dz});
  }
}

void ToroidField::buildGeometry() {
  // ---- Barrel base curve ----
  const double rEndB = 0.5 * m_cfg.barrel.b;
  const double lB = 0.5 * (m_cfg.barrel.R_in + m_cfg.barrel.R_out);
  const double aB = (m_cfg.barrel.R_out - m_cfg.barrel.R_in) - 2.0 * rEndB;

  const auto rz_barrel =
      racetrackRZ(aB, m_cfg.barrel.b, m_cfg.barrel.c, m_cfg.layout.nArc,
                  m_cfg.layout.nStraight, m_cfg.layout.closeLoop);
  std::vector<std::array<double, 2>> d_rzB;
  std::vector<std::array<double, 2>> m_rzB;
  buildSegsMidsRZ(rz_barrel, d_rzB, m_rzB);

  // ---- ECT base curve ----
  const double lE = 0.5 * (m_cfg.ect.R_in + m_cfg.ect.R_out) / UnitConstants::m;
  const double rEndECT = 0.5 * m_cfg.ect.b / UnitConstants::m;
  const double aE =
      (m_cfg.ect.R_out - m_cfg.ect.R_in) / UnitConstants::m - 2.0 * rEndECT;
  const auto rz_ect = racetrackRZ(
      /*a=*/aE, /*b=*/m_cfg.ect.b / UnitConstants::m,
      /*Lz=*/m_cfg.ect.c / UnitConstants::m, m_cfg.layout.nArc,
      m_cfg.layout.nStraight, m_cfg.layout.closeLoop);
  std::vector<std::array<double, 2>> d_rzE;
  std::vector<std::array<double, 2>> m_rzE;
  buildSegsMidsRZ(rz_ect, d_rzE, m_rzE);

  // ---- Angles ----
  const int nC = m_cfg.layout.nCoils;
  const double th0 = static_cast<float>(m_cfg.layout.theta0);
  const double dth = static_cast<float>(m_cfg.layout.thetaStep);

  // ---- Reserve & fill ----
  m_mids_barrel.clear();
  m_segs_barrel.clear();
  m_mids_ect.clear();
  m_segs_ect.clear();

  m_mids_barrel.reserve(static_cast<std::size_t>(nC) * m_rzB.size());
  m_segs_barrel.reserve(static_cast<std::size_t>(nC) * d_rzB.size());
  m_mids_ect.reserve(static_cast<std::size_t>(2 * nC) * m_rzE.size());
  m_segs_ect.reserve(static_cast<std::size_t>(2 * nC) * d_rzE.size());

  // Barrel rings with signs
  for (int k = 0; k < nC; ++k) {
    const double phi = th0 + static_cast<double>(k) * dth;
    const int sign = m_cfg.barrelSigns[k];
    mapRingToXYZ(lB, m_rzB, d_rzB, phi, sign, /*zShift=*/0.0, m_mids_barrel,
                 m_segs_barrel);
  }

  // Endcap centers (overlap placement)
  const float zECT =
      0.5f * static_cast<float>(m_cfg.barrel.c / UnitConstants::m) -
      0.5f * static_cast<float>(m_cfg.ect.c / UnitConstants::m) +
      2.0f * static_cast<float>(m_cfg.ect.gap / UnitConstants::m);

  // +z endcap (indices 0..nC-1 in ectSigns)
  for (int k = 0; k < nC; ++k) {
    const double phi = th0 + static_cast<double>(k) * dth;
    const int sign = m_cfg.ectSigns[k];
    mapRingToXYZ(lE, m_rzE, d_rzE, phi, sign, /*zShift=*/+zECT, m_mids_ect,
                 m_segs_ect);
  }
  // -z endcap (indices nC..2*nC-1 in ectSigns)
  for (int k = 0; k < nC; ++k) {
    const double phi = th0 + static_cast<double>(k) * dth;
    const int sign = m_cfg.ectSigns[nC + k];
    mapRingToXYZ(lE, m_rzE, d_rzE, phi, sign, /*zShift=*/-zECT, m_mids_ect,
                 m_segs_ect);
  }
}

}  // namespace Acts
