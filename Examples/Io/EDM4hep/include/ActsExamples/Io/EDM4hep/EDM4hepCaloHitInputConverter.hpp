// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/CaloHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepSimInputConverter.hpp"
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ActsExamples {

/// Per-collection calorimeter detector codes (barrel vs endcap @c z split).
struct CaloCollectionDetectorCodes {
  bool isBarrel = true;
  std::uint8_t barrelCode = 255;
  std::uint8_t endcapNegCode = 255;
  std::uint8_t endcapPosCode = 255;

  static CaloCollectionDetectorCodes barrel(std::uint8_t code);
  static CaloCollectionDetectorCodes endcap(std::uint8_t negZ,
                                            std::uint8_t posZ);

  /// Resolve the encoded detector code given a hit's z coordinate.
  /// For barrel collections, the code is independent of z.
  /// For endcaps, the sign of z selects negative vs positive endcap.
  std::uint8_t detectorCode(double z) const {
    if (isBarrel) {
      return barrelCode;
    }
    return z < 0.0 ? endcapNegCode : endcapPosCode;
  }
};

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
  static std::function<bool(std::string_view)> defaultIsEcalCollection();

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
    /// time-of-flight correction
    double ecalTimeMin = -1.0 * Acts::UnitConstants::ns;
    double ecalTimeMax = 10.0 * Acts::UnitConstants::ns;
    double hcalTimeMin = -1.0 * Acts::UnitConstants::ns;
    double hcalTimeMax = 10.0 * Acts::UnitConstants::ns;

    /// Time-of-flight correction offset. The corrected
    /// time is @c t - (r/c - tofOffset).
    double tofOffset = 0.1 * Acts::UnitConstants::ns;

    /// Maps each @c SimCalorimeterHit collection name to detector codes.
    /// At construction, entries are aligned to @c inputCaloHitCollections
    /// order (by index); missing names get code @c 255 (barrel unknown).
    /// No per-hit @c std::function — the hot loop uses collection index only.
    std::unordered_map<std::string, CaloCollectionDetectorCodes>
        caloDetectorCodesByCollectionName;

    /// Decides which time window applies to a collection. Returns @c true
    /// for ECal-style collections and @c false for HCal-style. The default
    /// checks for an "ECal" / "HCal" prefix and throws otherwise.
    std::function<bool(std::string_view collectionName)> isEcalCollection =
        defaultIsEcalCollection();
  };

  explicit EDM4hepCaloHitInputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

  ProcessCode finalize() override;

 private:
  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const override;

  /// Resolved once in the ctor: same order and length as
  /// @c m_cfg.inputCaloHitCollections.
  static std::vector<CaloCollectionDetectorCodes>
  buildDetectorCodesByCollectionIndex(
      const std::vector<std::string>& inputCaloHitCollections,
      const std::unordered_map<std::string, CaloCollectionDetectorCodes>&
          caloDetectorCodesByCollectionName);

  Config m_cfg;
  std::vector<CaloCollectionDetectorCodes> m_detectorCodesByCollectionIndex;

  ReadDataHandle<EDM4hepMCParticleIndexMap> m_inputMCParticleMap{
      this, "InputMCParticleMap"};
  WriteDataHandle<CaloHitContainer> m_outputCaloHits{this, "OutputCaloHits"};

  /// Per-section job-lifetime averaging timers. Each event takes one
  /// sample per section; the dtor logs aggregate stats at end-of-job.
  /// Wrapped in @c std::optional so the heavy timer type doesn't need a
  /// default constructor — they're emplaced in the converter ctor with the
  /// algorithm's logger.
  mutable std::optional<Acts::AveragingScopedTimer> m_timerResolve;
  mutable std::optional<Acts::AveragingScopedTimer> m_timerAggregate;
  mutable std::optional<Acts::AveragingScopedTimer> m_timerFinalise;
};

}  // namespace ActsExamples
