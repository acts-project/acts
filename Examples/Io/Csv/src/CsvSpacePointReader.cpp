// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"

#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>

#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <vector>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvSpacePointReader::CsvSpacePointReader(
    const ActsExamples::CsvSpacePointReader::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(
          determineEventFilesRange(cfg.inputDir, cfg.inputStem + ".csv")),
      m_logger(Acts::getDefaultLogger("CsvSpacePointReader", lvl)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
}

std::string ActsExamples::CsvSpacePointReader::CsvSpacePointReader::name() const {
  return "CsvSpacePointReader";
}

std::pair<size_t, size_t> ActsExamples::CsvSpacePointReader::availableEvents()
    const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::CsvSpacePointReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  SimSpacePointContainer spacePointsPixel;
  SimSpacePointContainer spacePointsStrip;
  SimSpacePointContainer spacePointsStripOverlap;
  
  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".csv",
                               ctx.eventNumber);

  dfe::NamedTupleCsvReader<SpacePointData> reader(path);
  SpacePointData data;

  while (reader.read(data)) {
    std::cout << data.measurement_id << ", " << data.sp_type << ", " << data.module_idhash << ", " << 
    data.sp_x << ", " << data.sp_y << ", " << data.sp_z << ", " << data.sp_radius << ", " << 
    data.sp_covr << ", " << data.sp_covz << std::endl;
    
    Acts::Vector3 globalPos(data.sp_x, data.sp_y, data.sp_z);
    
    if (data.sp_type==0)
      spacePointsPixel.emplace_back(globalPos, data.sp_covr, data.sp_covz, data.measurement_id);
    else if (data.sp_type==1)
      spacePointsStrip.emplace_back(globalPos, data.sp_covr, data.sp_covz, data.measurement_id);
    else if (data.sp_type==2)
      spacePointsStripOverlap.emplace_back(globalPos, data.sp_covr, data.sp_covz, data.measurement_id);
    else {
      ACTS_ERROR("Invalid space point type " << data.sp_type);
      return ProcessCode::ABORT;
    }
  }

  ACTS_DEBUG("Created " << spacePointsPixel.size() << " pixel space points");
  ACTS_DEBUG("Created " << spacePointsStrip.size() << " strip space points");
  ACTS_DEBUG("Created " << spacePointsStripOverlap.size() << " overlap space points");
  
  ctx.eventStore.add("PixelSpacePoints", std::move(spacePointsPixel));
  ctx.eventStore.add("StripSpacePoints", std::move(spacePointsStrip));
  ctx.eventStore.add("OverlapSpacePoints", std::move(spacePointsStripOverlap));
  
  return ProcessCode::SUCCESS;
}
