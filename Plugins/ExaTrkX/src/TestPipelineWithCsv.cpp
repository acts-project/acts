// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// not sure where to put this

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFinding.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

int main(int argc, char** argv) {
  assert(argc >= 3);
  std::vector<std::string> args(argv, argv + argc);
  std::ifstream in_file(args.at(1));

  std::vector<float> input_values;
  std::vector<int> spacepoint_ids;

  std::string header;
  std::getline(in_file, header);

  for (std::string line; in_file; std::getline(in_file, line)) {
    std::stringstream lineStream(line);
    std::vector<std::string> result;

    if (line.empty()) {
      continue;
    }

    for (std::string cell; lineStream; std::getline(lineStream, cell, ',')) {
      if (cell.empty()) {
        continue;
      }

      result.push_back(cell);
    }

    spacepoint_ids.push_back(std::stoi(result[0]));
    input_values.push_back(std::stof(result[1]));
    input_values.push_back(std::stof(result[2]));
    input_values.push_back(std::stof(result[3]));
  }

  auto sp_it = spacepoint_ids.begin();
  auto val_it = input_values.begin();
  std::cout << "Sample of imported values:\n";
  for (auto i = 0ul; i < 5; ++i) {
    std::cout << *sp_it << "\t";
    ++sp_it;
    for (auto j = 0ul; j < 3; ++j) {
      std::cout << *val_it << "\t";
      ++val_it;
    }
    std::cout << "\n";
  }
  std::cout << "\n\n";

  Acts::ExaTrkXTrackFinding::Config cfg;
  cfg.verbose = (std::find(args.begin(), args.end(), "-v") != args.end());
  cfg.modelDir = args.at(2);
  cfg.spacepointFeatures = 3;
  cfg.embeddingDim = 8;
  cfg.rVal = 0.15;
  cfg.knnVal = 500;
  cfg.filterCut = 0.05;
  cfg.edgeCut = 0.5;
  Acts::ExaTrkXTrackFinding pipeline(cfg);

  std::vector<std::vector<int>> track_candidates;
  ExaTrkXTime timing_info{};
  pipeline.getTracks(input_values, spacepoint_ids, track_candidates,
                     timing_info);

  std::cout << "Num track candidates:" << track_candidates.size() << "\n";
}
