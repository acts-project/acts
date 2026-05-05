// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/CaloHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimInputConverter.hpp"
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace ActsExamples {

/// Convert EDM4hep @c SimCalorimeterHit + @c CaloHitContribution collections
/// into the internal @c CaloHitContainer EDM, mirroring what
/// @c convert_calorimeter.py from ColliderML-Production does up to (but not
/// including) the per-detector energy threshold.
///
/// For every cell in the configured input collections, the converter walks
/// the cell's contributions, applies a time-of-flight correction
/// (@c corrected_time = t - (r/lightSpeed - tofOffset), where @c r is the
/// distance from the origin to the cell centre — matching pyedm4hep, which
/// also takes @c r from the parent @c SimCalorimeterHit position rather
/// than the per-contribution @c stepPosition that many simulations leave
/// unset) and a per-detector timing window, then aggregates surviving
/// contributions by source @c MCParticle:
///   - @c energy is the sum of contribution energies for the same
///     @c (cellID, particle) pair,
///   - @c time is the energy-weighted average,
///   - @c particleRow is looked up in the index map produced by
///     @c EDM4hepSimInputConverter so it matches what
///     @c ArrowParticleOutputConverter writes.
///
/// Cell-level energy thresholds are intentionally NOT applied here; the
/// downstream output converter chooses how to drop cells.
class EDM4hepCaloHitInputConverter final : public PodioInputConverter {
 public:
  /// Maps a calorimeter collection name + a hit's z coordinate to the
  /// @c uint8 detector enum stored on each @c CaloHit. The default mirrors
  /// the Python @c encode_calo_detector helper.
  using DetectorEncoder =
      std::function<std::uint8_t(std::string_view collectionName, double z)>;

  static DetectorEncoder defaultDetectorEncoder();

  struct Config {
    /// Whiteboard key for the input @c podio::Frame (typically the output
    /// of @c PodioReader).
    std::string inputFrame = "events";
    /// EDM4hep @c SimCalorimeterHit collection names to read. Each name
    /// implies a parallel @c "<name>Contributions" collection in the same
    /// frame.
    std::vector<std::string> inputCaloHitCollections;
    /// Whiteboard key of the @c EDM4hepMCParticleIndexMap that
    /// @c EDM4hepSimInputConverter publishes via its
    /// @c outputMCParticleMap option.
    std::string inputMCParticleMap;
    /// Output whiteboard key for the resulting @c CaloHitContainer.
    std::string outputCaloHits = "calohits";

    /// Per-detector contribution-level time windows applied after the
    /// time-of-flight correction, in ns. Defaults match the Python.
    double ecalTimeMin = -1.0;
    double ecalTimeMax = 10.0;
    double hcalTimeMin = -1.0;
    double hcalTimeMax = 10.0;

    /// Time-of-flight correction parameters: corrected time is
    /// @c t - (r/lightSpeed - tofOffset). Defaults match the Python
    /// (300 mm/ns, -0.1 ns).
    double lightSpeed = 300.0;
    double tofOffset = 0.1;

    /// Maps a collection name + z coordinate to the @c uint8 detector enum
    /// stored on each @c CaloHit. See @c DetectorEncoder.
    DetectorEncoder detectorEncoder = defaultDetectorEncoder();

    /// Decides which time window applies to a collection. Returns @c true
    /// for ECal-style collections and @c false for HCal-style. The default
    /// checks for an "ECal" / "HCal" prefix and throws otherwise.
    std::function<bool(std::string_view collectionName)> isEcalCollection =
        [](std::string_view name) {
          if (name.starts_with("ECal")) {
            return true;
          }
          if (name.starts_with("HCal")) {
            return false;
          }
          throw std::invalid_argument(
              std::string("EDM4hepCaloHitInputConverter: cannot infer "
                          "ECal/HCal from collection name '") +
              std::string(name) + "'");
        };
  };

  explicit EDM4hepCaloHitInputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

 private:
  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const override;

  Config m_cfg;

  ReadDataHandle<EDM4hepMCParticleIndexMap> m_inputMCParticleMap{
      this, "InputMCParticleMap"};
  WriteDataHandle<CaloHitContainer> m_outputCaloHits{this, "OutputCaloHits"};
};

}  // namespace ActsExamples
