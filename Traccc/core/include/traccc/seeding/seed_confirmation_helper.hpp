/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"

namespace traccc {

/// Helper functions for the seed weights and the seed confirmation of the
/// Acts triplet seed filter, shared by the host and device seed selection
struct seed_confirmation_helper {
  /// Radius of a top spacepoint used to look for compatible seeds
  ///
  /// @param cfg the seed filter configuration
  /// @param spM_r the radius of the middle spacepoint
  /// @param spM_z the z of the middle spacepoint
  /// @param spT_r the radius of the top spacepoint
  /// @param spT_z the z of the top spacepoint
  ///
  static TRACCC_HOST_DEVICE scalar top_radius(const seedfilter_config& cfg,
                                              scalar spM_r, scalar spM_z,
                                              scalar spT_r, scalar spT_z) {
    if (cfg.useDeltaRinsteadOfTopRadius) {
      const scalar dr = spT_r - spM_r;
      const scalar dz = spT_z - spM_z;
      return math::sqrt(dr * dr + dz * dz);
    }
    return spT_r;
  }

  /// Weight terms that depend on the number of compatible top spacepoints
  /// and on the seed direction
  ///
  /// @param cfg the seed filter configuration
  /// @param n_compat the number of compatible top spacepoints
  /// @param impact the transverse impact parameter of the seed
  /// @param z_vertex the z origin of the seed
  /// @param cot_theta the cot(theta) of the middle-bottom doublet
  ///
  static TRACCC_HOST_DEVICE scalar extra_weight(const seedfilter_config& cfg,
                                                unsigned int n_compat,
                                                scalar impact, scalar z_vertex,
                                                scalar cot_theta) {
    scalar weight = 0.f;
    if (static_cast<scalar>(n_compat) > cfg.numSeedIncrement) {
      weight += cfg.seedWeightIncrement;
    }
    if (cfg.absDeltaEtaWeightFactor > 0.f &&
        impact > cfg.absDeltaEtaMinImpact) {
      const scalar etaSeed = math::asinh(cot_theta);
      const scalar etaPCA = math::asinh(z_vertex / impact);
      weight -= cfg.absDeltaEtaWeightFactor * math::fabs(etaSeed - etaPCA);
    }
    return weight;
  }

  /// The seed confirmation settings for a spacepoint at a given z
  static TRACCC_HOST_DEVICE const seed_confirmation_range& range(
      const seedfilter_config& cfg, scalar z) {
    const bool forward = z > cfg.centralSeedConfirmationRange.zMaxSeedConf ||
                         z < cfg.centralSeedConfirmationRange.zMinSeedConf;
    return forward ? cfg.forwardSeedConfirmationRange
                   : cfg.centralSeedConfirmationRange;
  }

  /// Number of compatible top spacepoints required for a spacepoint at a
  /// given radius
  static TRACCC_HOST_DEVICE unsigned int n_top_required(
      const seed_confirmation_range& rng, scalar r) {
    return r > rng.rMaxSeedConf ? rng.nTopForLargeR : rng.nTopForSmallR;
  }

  /// Minimum number of top candidates of a middle-bottom doublet for its
  /// triplets to be considered
  static TRACCC_HOST_DEVICE unsigned int min_top_candidates(
      const seedfilter_config& cfg, scalar spB_r, scalar spM_z,
      unsigned int n_quality) {
    unsigned int n = 2u;
    if (!cfg.seedConfirmation || spB_r > range(cfg, spM_z).rMaxSeedConf) {
      n = 1u;
    }
    if (cfg.seedConfirmation && n_quality > 0u) {
      ++n;
    }
    return n;
  }

  /// Seed confirmation of a triplet
  ///
  /// @param cfg the seed filter configuration
  /// @param rng the confirmation settings of the bottom spacepoint region
  /// @param spB_r the radius of the bottom spacepoint
  /// @param impact the transverse impact parameter of the seed
  /// @param z_vertex the z origin of the seed
  /// @param n_compat the number of compatible top spacepoints
  /// @param n_quality the number of quality seeds already found for the
  ///                  middle spacepoint
  /// @param[in,out] weight the seed weight
  /// @param[out] quality whether the seed is a quality seed
  /// @return false if the seed is rejected
  ///
  static TRACCC_HOST_DEVICE bool confirm(const seedfilter_config& cfg,
                                         const seed_confirmation_range& rng,
                                         scalar spB_r, scalar impact,
                                         scalar z_vertex, unsigned int n_compat,
                                         unsigned int n_quality, scalar& weight,
                                         bool& quality) {
    const unsigned int nTopSeedConf = n_top_required(rng, spB_r);
    const int deltaSeedConf =
        static_cast<int>(n_compat) + 1 - static_cast<int>(nTopSeedConf);
    if (deltaSeedConf < 0 || (n_quality != 0u && deltaSeedConf == 0)) {
      return false;
    }
    const bool seedRangeCuts = spB_r < rng.seedConfMinBottomRadius ||
                               math::fabs(z_vertex) > rng.seedConfMaxZOrigin;
    if (seedRangeCuts && deltaSeedConf == 0 && impact > rng.minImpactSeedConf) {
      return false;
    }
    weight += -(math::fabs(z_vertex) * cfg.zOriginWeightFactor) +
              cfg.compatSeedWeight;
    quality = deltaSeedConf > 0;
    return true;
  }
};

}  // namespace traccc
