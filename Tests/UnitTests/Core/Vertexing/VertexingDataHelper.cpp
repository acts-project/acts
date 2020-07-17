// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fstream>
#include <iterator>
#include <regex>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;
using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoCtx = GeometryContext();

/// @brief Helper struct to store reference vertex related information
struct VertexInfo {
  // The position
  Vector3D position;
  // The covariance
  ActsSymMatrixD<3> covariance;
  // Number of tracks
  int nTracks;
  // Weight of first track
  double trk1Weight;
  // Vertex compatibility value of first track
  double trk1Comp;
  // Chi2 of first track
  double trk1Chi2;
};

std::tuple<Vertex<BoundParameters>, std::vector<VertexInfo>,
           std::vector<BoundParameters>>
readTracksAndVertexCSV(std::string file) {
  const std::regex comma(",");

  // Open source file
  std::ifstream mesh(file);

  // Here we will store the result
  std::vector<std::vector<std::string>> point_coordinates;

  // We want to read all lines of the file
  std::string line{};
  bool isBeamSpot = false;
  bool isTrack = false;
  bool isVertex = false;

  std::shared_ptr<PerigeeSurface> perigeeSurface;
  std::vector<BoundParameters> tracks;
  std::vector<VertexInfo> vertices;
  Vertex<BoundParameters> beamspotConstraint;

  while (mesh && getline(mesh, line)) {
    // Tokenize line and store result in vector
    std::vector<std::string> row{
        std::sregex_token_iterator(line.begin(), line.end(), comma, -1),
        std::sregex_token_iterator()};

    if (row[0] == std::string("beamspot")) {
      isBeamSpot = true;
      continue;
    }
    if (row[0] == "tracks") {
      isTrack = true;
      isBeamSpot = false;
      continue;
    }
    if (row[0] == "vertices") {
      isTrack = false;
      isBeamSpot = false;
      isVertex = true;
      continue;
    }

    if (isBeamSpot) {
      Vector3D beamspotPos;
      ActsSymMatrixD<3> beamspotCov;
      beamspotPos << std::stod(row[0]) * (1_mm), std::stod(row[1]) * (1_mm),
          std::stod(row[2]) * (1_mm);
      beamspotCov << std::stod(row[3]), 0, 0, 0, std::stod(row[4]), 0, 0, 0,
          std::stod(row[5]);
      beamspotConstraint.setPosition(beamspotPos);
      beamspotConstraint.setCovariance(beamspotCov);
      perigeeSurface = Surface::makeShared<PerigeeSurface>(beamspotPos);
    }

    if (isTrack) {
      BoundVector params;
      params << std::stod(row[0]), std::stod(row[1]), std::stod(row[2]),
          std::stod(row[3]), std::stod(row[4]) * 1. / (1_MeV),
          std::stod(row[5]);
      Covariance covMat;
      covMat << std::stod(row[6]), std::stod(row[7]), std::stod(row[8]),
          std::stod(row[9]), std::stod(row[10]) * 1. / (1_MeV),
          std::stod(row[11]), std::stod(row[12]), std::stod(row[13]),
          std::stod(row[14]), std::stod(row[15]),
          std::stod(row[16]) * 1. / (1_MeV), std::stod(row[17]),
          std::stod(row[18]), std::stod(row[19]), std::stod(row[20]),
          std::stod(row[21]), std::stod(row[22]) * 1. / (1_MeV),
          std::stod(row[23]), std::stod(row[24]), std::stod(row[25]),
          std::stod(row[26]), std::stod(row[27]),
          std::stod(row[28]) * 1. / (1_MeV), std::stod(row[29]),
          std::stod(row[30]) * 1. / (1_MeV), std::stod(row[31]) * 1. / (1_MeV),
          std::stod(row[32]) * 1. / (1_MeV), std::stod(row[33]) * 1. / (1_MeV),
          std::stod(row[34]) * 1. / (1_MeV * 1_MeV), std::stod(row[35]),
          std::stod(row[36]), std::stod(row[37]), std::stod(row[38]),
          std::stod(row[39]), std::stod(row[40]), std::stod(row[41]);

      auto boundParams =
          BoundParameters(geoCtx, std::move(covMat), params, perigeeSurface);
      tracks.push_back(boundParams);
    }

    if (isVertex) {
      Vector3D pos;
      pos << std::stod(row[0]) * (1_mm), std::stod(row[1]) * (1_mm),
          std::stod(row[2]) * (1_mm);
      ActsSymMatrixD<3> cov;
      cov << std::stod(row[3]), std::stod(row[4]), std::stod(row[5]),
          std::stod(row[6]), std::stod(row[7]), std::stod(row[8]),
          std::stod(row[9]), std::stod(row[10]), std::stod(row[11]);
      VertexInfo vertexInfo;
      vertexInfo.position = pos;
      vertexInfo.covariance = cov;
      vertexInfo.nTracks = std::stoi(row[12]);
      vertexInfo.trk1Weight = std::stod(row[13]);
      vertexInfo.trk1Comp = std::stod(row[14]);
      vertexInfo.trk1Chi2 = std::stod(row[15]);
      vertices.push_back(vertexInfo);
    }
  }

  return std::make_tuple(beamspotConstraint, vertices, tracks);
}

}  // namespace Test
}  // namespace Acts