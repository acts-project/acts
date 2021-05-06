// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

#include <TApplication.h>

#include "boundParamResolution.C"

int main(int argc, char **argv) {
  if (argc > 1) {
    std::string farg = argv[1];
    if (farg == "-h" or farg.find("help") != std::string::npos) {
      std::cout << "*** ACTS Residual and Pull plotting " << std::endl;
      std::cout << "*** Usage: ./ActsExampleResidualAndPulls input_file "
                   "input_tree output_file <save_extension>"
                << std::endl;
    }
    if (argc >= 4) {
      TApplication theApp("ResidualAndPulls", 0, 0);
      std::string inputFile = farg;
      std::string inputTree = argv[2];
      std::string outputFile = argv[3];
      std::string saveAs = (argc > 4) ? std::string(argv[4]) : std::string("");
      boundParamResolution(inputFile, inputTree, outputFile, saveAs);
      theApp.Run();
    }
  }
  return 1;
}
