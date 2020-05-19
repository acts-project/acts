// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @brief Plain structs that each define one row in a TrackML csv file

#pragma once

#include <cstdint>
#include <dfe/dfe_namedtuple.hpp>

namespace FW {

struct ParticleData {
  /// Event-unique particle identifier a.k.a barcode.
  uint64_t particle_id;
  /// Particle type number a.k.a. PDG particle number.
  int32_t particle_type;
  /// Production process type. Not available in the TrackML datasets.
  uint32_t process = 0u;
  /// Production position components in mm.
  float vx, vy, vz;
  // Production time in ns. Not available in the TrackML datasets.
  float vt = 0.0f;
  /// Momentum components in GeV.
  float px, py, pz;
  /// Mass in GeV. Not available in the TrackML datasets
  float m = 0.0f;
  /// Charge in e.
  float q;

  DFE_NAMEDTUPLE(ParticleData, particle_id, particle_type, process, vx, vy, vz,
                 vt, px, py, pz, m, q);
};

struct TruthHitData {
  /// Event-unique hit identifier. As defined for the simulated hit below and
  /// used to link back to it; same value can appear multiple times here due to
  /// shared hits in dense environments.
  uint64_t hit_id;
  /// Hit surface identifier. Not available in the TrackML datasets.
  uint64_t geometry_id = 0u;
  /// Event-unique particle identifier of the generating particle.
  uint64_t particle_id;
  /// True global hit position components in mm.
  float tx, ty, tz;
  // True global hit time in ns. Not available in the TrackML datasets.
  float tt = 0.0f;
  /// True particle momentum in GeV before interaction.
  float tpx, tpy, tpz;
  /// True particle energy in GeV before interaction.
  /// Not available in the TrackML datasets.
  float te = 0.0f;
  /// True four-momentum change in GeV due to interaction.
  /// Not available in the TrackML datasets.
  float deltapx = 0.0f;
  float deltapy = 0.0f;
  float deltapz = 0.0f;
  float deltae = 0.0f;
  // Hit index along the trajectory. Not available in the TrackML datasets.
  int32_t index = -1;

  DFE_NAMEDTUPLE(TruthHitData, hit_id, particle_id, geometry_id, tx, ty, tz, tt,
                 tpx, tpy, tpz, te, deltapx, deltapy, deltapz, deltae, index);
};

struct HitData {
  /// Event-unique hit identifier. Each value can appear at most once.
  uint64_t hit_id;
  /// Hit surface identifier. Not available in the TrackML datasets.
  uint64_t geometry_id = 0u;
  /// Partially decoded hit surface identifier components.
  uint32_t volume_id, layer_id, module_id;
  /// Global hit position components in mm.
  float x, y, z;
  /// Global hit time in ns. Not available in the TrackML datasets.
  float t = 0.0f;

  DFE_NAMEDTUPLE(HitData, hit_id, geometry_id, volume_id, layer_id, module_id,
                 x, y, z, t);
};

struct CellData {
  /// Event-unique hit identifier. As defined for the simulated hit above and
  /// used to link back to it; same value can appear multiple times for clusters
  /// with more than one active cell.
  uint64_t hit_id;
  /// Digital cell address/ channel identifier. These should have been named
  /// channel{0,1} but we cannot change it now to avoid breaking backward
  /// compatibility.
  int32_t ch0, ch1;
  /// Digital cell timestamp. Not available in the TrackML datasets.
  int32_t timestamp = 0;
  /// (Digital) measured cell value, e.g. amplitude or time-over-threshold.
  int32_t value;

  DFE_NAMEDTUPLE(CellData, hit_id, ch0, ch1, timestamp, value);
};

struct SurfaceData {
  /// Surface identifier. Not available in the TrackML datasets.
  uint64_t geometry_id;
  /// Partially decoded surface identifier components.
  uint32_t volume_id, layer_id, module_id;
  /// Center position components in mm.
  float cx, cy, cz;
  /// Rotation matrix components.
  float rot_xu, rot_xv, rot_xw;
  float rot_yu, rot_yv, rot_yw;
  float rot_zu, rot_zv, rot_zw;
  /// Limits and pitches in mm. Not always available.
  float module_t = -1;
  float module_minhu = -1;
  float module_maxhu = -1;
  float module_hv = -1;
  float pitch_u = -1;
  float pitch_v = -1;

  DFE_NAMEDTUPLE(SurfaceData, geometry_id, volume_id, layer_id, module_id, cx,
                 cy, cz, rot_xu, rot_xv, rot_xw, rot_yu, rot_yv, rot_yw, rot_zu,
                 rot_zv, rot_zw, module_t, module_minhu, module_maxhu,
                 module_hv, pitch_u, pitch_v);
};

}  // namespace FW
