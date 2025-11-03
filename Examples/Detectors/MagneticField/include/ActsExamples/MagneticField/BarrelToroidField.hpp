// This file is part of the ACTS project.
//
// BarrelToroidField.hpp
//
// A discretized Biot–Savart magnetic field for the ATLAS barrel toroid,
// using racetrack coils (straight legs along z joined by semicircles in the ρ–z plane)
// replicated at 8 azimuthal positions.
//
// The implementation mirrors the provided Python reference:
// - Geometry parameters (R_in, R_out, coil length c, end-arc diameter b, current I, turns N)
// - 8 racetrack coils at angles 22.5° + k*45°
// - Polyline discretization of the racetrack in local (ρ,z), mapped to (x,y,z)
// - Field evaluation: B = μ0 I / (4π) ∑ (dl × r) / |r|^3 using segment midpoints
//
// Notes:
// * Internal precomputation uses float for storage; accumulation is in double.
// * Returned Acts::Vector3 is in the global coordinate system (Tesla).
// * This is header-only for convenience; you can move the implementation to a .cpp if preferred.

#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

class BarrelToroidField final : public MagneticFieldProvider {
 public:
  /// Configuration for geometry, current, and discretization
  struct Config {
    // --- Geometry & current (ATLAS-like defaults from the Python snippet) ---
    float R_in   = 4.9f;      // [m]
    float R_out  = 10.0f;     // [m]
    float c      = 25.3f;     // [m] coil length along z (long straight legs)
    float b      = 0.16f;     // [m] end-arc diameter in ρ–z plane (r_end = b/2)
    float I      = 20500.0f;  // [A] current per turn
    int   Nturns = 120;       // total turns per coil

    // --- Discretization of one racetrack in local (ρ,z) ---
    int nArc      = 200;      // points on each semicircle
    int nStraight = 160;      // points per straight leg
    bool closeLoop = true;    // close polyline

    // --- Ring placement (as Python: start at 22.5°, step 45°, 8 coils) ---
    float theta0_deg   = 22.5f;
    float thetaStep_deg= 45.0f;
    int   nCoils       = 8;

    // Small epsilon to regularize 1/r^3
    float eps = 1e-18f;
  };

  struct Cache {
    explicit Cache(const MagneticFieldContext& /*ctx*/) {}
  };

  /// Default constructor with default config
  BarrelToroidField() : BarrelToroidField(Config{}) {}

  /// Constructor with explicit config
  explicit BarrelToroidField(const Config& cfg)
      : m_cfg(cfg) {
    buildGeometry();
  }

  /// Make the cache (opaque wrapper provided by base)
  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  /// Field query at a position (global coordinates, meters); returns Tesla.
  Result<Vector3> getField(const Vector3& position,
                           MagneticFieldProvider::Cache& cache) const override {
    (void)cache;

    // Position components (double -> float for computation)
    const double Xd = position.x();
    const double Yd = position.y();
    const double Zd = position.z();

    const float X = static_cast<float>(Xd);
    const float Y = static_cast<float>(Yd);
    const float Z = static_cast<float>(Zd);

    // Biot–Savart accumulation in double
    double bx = 0.0, by = 0.0, bz = 0.0;

    // Prefactor μ0 * (N * I) / (4π)
    constexpr double mu0 = 4e-7 * M_PI; // [T·m/A]
    const double Ieff = static_cast<double>(m_cfg.Nturns) * static_cast<double>(m_cfg.I);
    const double pref = mu0 * Ieff / (4.0 * M_PI);
    const float  eps  = m_cfg.eps;

    // Loop over all segments
    const std::size_t M = m_segs.size();
    for (std::size_t i = 0; i < M; ++i) {
      const auto& mid = m_mids[i];
      const auto& dl  = m_segs[i];

      const float rx = X - mid[0];
      const float ry = Y - mid[1];
      const float rz = Z - mid[2];

      const double r2 = static_cast<double>(rx)*rx
                      + static_cast<double>(ry)*ry
                      + static_cast<double>(rz)*rz
                      + static_cast<double>(eps);

      const double invr  = 1.0 / std::sqrt(r2);
      const double invr3 = invr / r2;

      // cross(dl, r) = dl x r
      const double cx = static_cast<double>(dl[1]) * rz - static_cast<double>(dl[2]) * ry;
      const double cy = static_cast<double>(dl[2]) * rx - static_cast<double>(dl[0]) * rz;
      const double cz = static_cast<double>(dl[0]) * ry - static_cast<double>(dl[1]) * rx;

      bx += cx * invr3;
      by += cy * invr3;
      bz += cz * invr3;
    }

    const Vector3 B(pref * bx, pref * by, pref * bz);
    return Result<Vector3>::success(B);
  }

  /// Coarse geometric inclusion check (optional); here always true for simplicity.
  bool isInside(const Vector3& /*position*/) const { return true; }

  /// Access to configuration actually used (after construction)
  const Config& config() const { return m_cfg; }

 private:
  // ===== Geometry construction (precompute segments & midpoints) =====
  void buildGeometry() {
    // Derived parameters (as in Python)
    const float R_in  = m_cfg.R_in;
    const float R_out = m_cfg.R_out;
    const float c     = m_cfg.c;
    const float b     = m_cfg.b;
    const float r_end = 0.5f * b;

    // Mean radius and radial spacing between straight legs (center-to-center)
    const float l = 0.5f * (R_in + R_out);
    const float a = (R_out - R_in) - 2.0f * r_end; // distance between straight legs (radial)
    // Build racetrack polyline in local (ρ,z)
    std::vector<std::array<float, 2>> rzLoop = racetrackRZ(a, b, c, m_cfg.nArc, m_cfg.nStraight, m_cfg.closeLoop);

    // Build segment vectors and midpoints in (ρ,z)
    std::vector<std::array<float, 2>> d_rz;
    std::vector<std::array<float, 2>> m_rz;
    d_rz.reserve(rzLoop.size() - 1);
    m_rz.reserve(rzLoop.size() - 1);
    for (std::size_t i = 0; i + 1 < rzLoop.size(); ++i) {
      const float dr = rzLoop[i + 1][0] - rzLoop[i][0];
      const float dz = rzLoop[i + 1][1] - rzLoop[i][1];
      d_rz.push_back({dr, dz});
      m_rz.push_back({0.5f * (rzLoop[i + 1][0] + rzLoop[i][0]),
                      0.5f * (rzLoop[i + 1][1] + rzLoop[i][1])});
    }

    // Place at nCoils azimuthal angles: θ = θ0 + k*Δθ
    const int nCoils = m_cfg.nCoils;
    const float th0  = deg2rad(m_cfg.theta0_deg);
    const float dth  = deg2rad(m_cfg.thetaStep_deg);

    m_segs.clear();
    m_mids.clear();
    m_segs.reserve(static_cast<std::size_t>(nCoils) * d_rz.size());
    m_mids.reserve(static_cast<std::size_t>(nCoils) * m_rz.size());

    for (int k = 0; k < nCoils; ++k) {
      const float t  = th0 + k * dth;
      const float ct = std::cos(t);
      const float st = std::sin(t);

      // Map local midpoints (ρ_mid, z_mid) -> global (x,y,z):
      // x = (l + ρ_mid) * cosθ; y = (l + ρ_mid) * sinθ; z = z_mid
      for (const auto& rm : m_rz) {
        const float rho = rm[0];
        const float zz  = rm[1];
        const float rxy = l + rho;
        m_mids.push_back({rxy * ct, rxy * st, zz});
      }

      // Map local segment vectors dℓ = dρ e_r(θ) + dz e_z
      // e_r(θ) = (cosθ, sinθ, 0)
      for (const auto& dlrz : d_rz) {
        const float dr = dlrz[0];
        const float dz = dlrz[1];
        m_segs.push_back({dr * ct, dr * st, dz});
      }
    }
  }

  static float deg2rad(float deg) {
    return deg * static_cast<float>(M_PI / 180.0);
  }

  // Build a single racetrack polyline in local (ρ, z)
  // Semicircles of radius r_end = b/2 at z = ±Lz/2; straight legs at ρ = ±a/2.
  static std::vector<std::array<float, 2>> racetrackRZ(float a, float b, float Lz,
                                                       int nArc, int nStraight, bool close) {
    const float r = 0.5f * b;
    const float z0 = 0.0f;

    std::vector<std::array<float, 2>> pts;
    pts.reserve(static_cast<std::size_t>(2 * nArc + 2 * nStraight + (close ? 1 : 0)));

    // We reuse the Python parameterization by swapping axes from a (x,y) racetrack:
    // Make a racetrack in (x≡z, y≡ρ): two semicircles centered at x=±Lz/2, y=0, radius r,
    // and straight segments connecting at y=±r between x from +Lz/2 to -Lz/2 and back.

    // Arc at x = +Lz/2: θ ∈ [π/2, 3π/2)  -> (x1=z, y1=ρ)
    for (int i = 0; i < nArc; ++i) {
      float th = static_cast<float>(M_PI_2) + (static_cast<float>(i) / nArc) * static_cast<float>(M_PI);
      float x = (+Lz * 0.5f) + r * std::cos(th);
      float y = 0.0f + r * std::sin(th);
      pts.push_back({y, x}); // (ρ, z)
    }

    // Straight leg at y = -r from x = +Lz/2 -> -Lz/2 (exclude endpoint to avoid duplication)
    for (int i = 0; i < nStraight; ++i) {
      float t = static_cast<float>(i) / static_cast<float>(nStraight);
      float x = (+Lz * 0.5f) * (1.0f - t) + (-Lz * 0.5f) * t;
      float y = -r;
      pts.push_back({y, x});
    }

    // Arc at x = -Lz/2: θ ∈ [3π/2, π/2)  -> (x3=z, y3=ρ)
    for (int i = 0; i < nArc; ++i) {
      float th = static_cast<float>(3.0 * M_PI_2) + (static_cast<float>(i) / nArc) * static_cast<float>(-M_PI);
      float x = (-Lz * 0.5f) + r * std::cos(th);
      float y = 0.0f + r * std::sin(th);
      pts.push_back({y, x});
    }

    // Straight leg at y = +r from x = -Lz/2 -> +Lz/2
    const int nStraightLast = close ? nStraight : (nStraight + 1);
    for (int i = 0; i < nStraightLast; ++i) {
      float t = static_cast<float>(i) / static_cast<float>(nStraightLast - 1);
      float x = (-Lz * 0.5f) * (1.0f - t) + (+Lz * 0.5f) * t;
      float y = +r;
      pts.push_back({y, x});
    }

    // Shift the straight legs to be at ρ = ±a/2 instead of ±r by adding (±a/2 ∓ r)
    // The Python racetrack in (x,y) used x ≡ z, y ≡ ρ with arcs radius r and straight legs at y = ±r.
    // We want the straight legs at ρ = ±a/2 while keeping arcs radius r about those leg centers.
    // That is achieved already by the construction where the *end-circle centers* are at ρ = 0 with radius r,
    // and the straight portions lie at ±r. To move them to ±a/2, we offset all ρ by (a/2 - r) in the + region
    // and by -(a/2 - r) in the - region. An easier equivalent: add a uniform offset Δ so that ±r → ±a/2:
    // Δ = a/2 - r. Then ρ := ρ + sign(ρ) * Δ.
    const float delta = 0.5f * a - r;
    for (auto& p : pts) {
      const float rho = p[0];
      const float sgn = (rho >= 0.0f) ? 1.0f : -1.0f;
      p[0] = rho + sgn * delta;
    }

    // Close loop explicitly if requested
    if (close) {
      if (pts.front()[0] != pts.back()[0] || pts.front()[1] != pts.back()[1]) {
        pts.push_back(pts.front());
      }
    }

    return pts;
  }

 private:
  Config m_cfg;

  // Precomputed for all coils: segment vectors dl and midpoints
  // Stored as float arrays {x,y,z}
  std::vector<std::array<float, 3>> m_segs;
  std::vector<std::array<float, 3>> m_mids;
};

}  // namespace Acts 