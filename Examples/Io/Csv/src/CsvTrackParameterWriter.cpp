// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvTrackParameterWriter.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <ios>
#include <optional>
#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvTrackParameterWriter::CsvTrackParameterWriter(
    const ActsExamples::CsvTrackParameterWriter::Config& config,
    Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("CsvTrackParameterWriter", level)) {
  if (m_cfg.inputTrackParameters.empty() == m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument(
        "You have to either provide track parameters or trajectories");
  }

  m_inputTrackParameters.maybeInitialize(m_cfg.inputTrackParameters);
  m_inputTrajectories.maybeInitialize(m_cfg.inputTrajectories);
}

ActsExamples::CsvTrackParameterWriter::~CsvTrackParameterWriter() = default;

std::string ActsExamples::CsvTrackParameterWriter::name() const {
  return "CsvTrackParameterWriter";
}

ActsExamples::ProcessCode ActsExamples::CsvTrackParameterWriter::finalize() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::CsvTrackParameterWriter::write(
    const AlgorithmContext& ctx) {
  std::vector<Acts::BoundTrackParameters> inputTrackParameters;

  if (!m_cfg.inputTrackParameters.empty()) {
    const auto& tmp = m_inputTrackParameters(ctx);
    inputTrackParameters = tmp;
  } else {
    const auto& inputTrajectories = m_inputTrajectories(ctx);

    for (const auto& trajectories : inputTrajectories) {
      for (auto tip : trajectories.tips()) {
        if (!trajectories.hasTrackParameters(tip)) {
          continue;
        }
        const auto& trackParam = trajectories.trackParameters(tip);
        inputTrackParameters.push_back(trackParam);
      }
    }
  }

  std::string path =
      perEventFilepath(m_cfg.outputDir, m_cfg.outputStem, ctx.eventNumber);

  dfe::NamedTupleCsvWriter<TrackParameterData> writer(path,
                                                      m_cfg.outputPrecision);

  TrackParameterData data{};
  for (const auto& tp : inputTrackParameters) {
    const auto& params = tp.parameters();
    const auto& cov = tp.covariance().value();

    data.d0 = params[Acts::eBoundLoc0];
    data.z0 = params[Acts::eBoundLoc1];
    data.phi = params[Acts::eBoundPhi];
    data.theta = params[Acts::eBoundTheta];
    data.qop = params[Acts::eBoundQOverP];

    data.var_d0 = cov(Acts::eBoundLoc0, Acts::eBoundLoc0);
    data.var_z0 = cov(Acts::eBoundLoc1, Acts::eBoundLoc1);
    data.var_phi = cov(Acts::eBoundPhi, Acts::eBoundPhi);
    data.var_theta = cov(Acts::eBoundTheta, Acts::eBoundTheta);
    data.var_qop = cov(Acts::eBoundQOverP, Acts::eBoundQOverP);

    data.cov_d0z0 = cov(Acts::eBoundLoc0, Acts::eBoundLoc1);
    data.cov_d0phi = cov(Acts::eBoundLoc0, Acts::eBoundPhi);
    data.cov_d0theta = cov(Acts::eBoundLoc0, Acts::eBoundTheta);
    data.cov_d0qop = cov(Acts::eBoundLoc0, Acts::eBoundQOverP);

    data.cov_z0d0 = cov(Acts::eBoundLoc1, Acts::eBoundLoc0);
    data.cov_z0phi = cov(Acts::eBoundLoc1, Acts::eBoundPhi);
    data.cov_z0theta = cov(Acts::eBoundLoc1, Acts::eBoundTheta);
    data.cov_z0qop = cov(Acts::eBoundLoc1, Acts::eBoundQOverP);

    data.cov_phid0 = cov(Acts::eBoundPhi, Acts::eBoundLoc0);
    data.cov_phiz0 = cov(Acts::eBoundPhi, Acts::eBoundLoc1);
    data.cov_phitheta = cov(Acts::eBoundPhi, Acts::eBoundTheta);
    data.cov_phiqop = cov(Acts::eBoundPhi, Acts::eBoundQOverP);

    data.cov_thetad0 = cov(Acts::eBoundTheta, Acts::eBoundLoc0);
    data.cov_thetaz0 = cov(Acts::eBoundTheta, Acts::eBoundLoc1);
    data.cov_thetaphi = cov(Acts::eBoundTheta, Acts::eBoundPhi);
    data.cov_thetaqop = cov(Acts::eBoundTheta, Acts::eBoundQOverP);

    data.cov_qopd0 = cov(Acts::eBoundQOverP, Acts::eBoundLoc0);
    data.cov_qopz0 = cov(Acts::eBoundQOverP, Acts::eBoundLoc1);
    data.cov_qopphi = cov(Acts::eBoundQOverP, Acts::eBoundPhi);
    data.cov_qoptheta = cov(Acts::eBoundQOverP, Acts::eBoundTheta);

    writer.append(data);
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
