// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/CaloHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/// Convert a @c CaloHitContainer into an @c arrow::Table that mirrors the
/// per-event layout produced by @c convert_calorimeter.py from
/// ColliderML-Production. The container is expected to come from an upstream
/// input converter (e.g. @c EDM4hepCaloHitInputConverter) that has already
/// applied any time-of-flight and timing-window cuts and aggregated
/// contributions per @c (cellID, particle).
///
/// The resulting table has one row per event with the following columns:
///   - @c detector : list<uint8>
///   - @c total_energy, @c x, @c y, @c z : list<float32>
///   - @c contrib_particle_ids : list<list<uint64>>
///   - @c contrib_energies, @c contrib_times : list<list<float32>>
///
/// @c event_id is added downstream by the @c ParquetWriter. Cells whose
/// @c totalEnergy is below the threshold returned by @c cellThreshold are
/// dropped at this stage.
class ACTS_ARROW_EXPORT ArrowCaloHitOutputConverter final
    : public ArrowOutputConverter {
 public:
  /// Returns the per-cell energy threshold (in the same units as
  /// @c CaloHit::totalEnergy) for a given subdetector code.
  using CellThresholdFn = std::function<double(std::uint8_t detector)>;

  /// Default threshold callback that mirrors the Python
  /// @c convert_calorimeter.py defaults: ECal codes 9..11 get
  /// @c ecalEnergyThreshold, HCal codes 12..14 get @c hcalEnergyThreshold,
  /// every other code gets 0 (i.e. nothing dropped).
  static CellThresholdFn defaultCellThreshold(double ecalEnergyThreshold,
                                              double hcalEnergyThreshold);

  struct Config {
    /// Whiteboard key of the input @c CaloHitContainer.
    std::string inputCaloHits;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "calohits";

    /// Cell-level energy threshold in ACTS energy units, evaluated per cell
    /// using its @c detector code. Defaults match the Python (50 keV ECal,
    /// 250 keV HCal).
    double ecalEnergyThreshold = 50.0 * Acts::UnitConstants::keV;
    double hcalEnergyThreshold = 250.0 * Acts::UnitConstants::keV;

    /// Custom threshold callback. If unset (default), one is constructed
    /// from @c ecalEnergyThreshold and @c hcalEnergyThreshold via
    /// @c defaultCellThreshold.
    CellThresholdFn cellThreshold;
  };

  explicit ArrowCaloHitOutputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;
  CellThresholdFn m_cellThreshold;

  ReadDataHandle<CaloHitContainer> m_inputCaloHits{this, "InputCaloHits"};
  WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
