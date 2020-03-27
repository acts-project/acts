// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/BField/BFieldUtils.hpp"

#include <fstream>

#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

Acts::InterpolatedBFieldMapper<
    Acts::detail::Grid<Acts::Vector2D, Acts::detail::EquidistantAxis,
                       Acts::detail::EquidistantAxis>>
FW::BField::txt::fieldMapperRZ(
    std::function<size_t(std::array<size_t, 2> binsRZ,
                         std::array<size_t, 2> nBinsRZ)>
        localToGlobalBin,
    std::string fieldMapFile, double lengthUnit, double BFieldUnit,
    size_t nPoints, bool firstQuadrant) {
  /// [1] Read in field map file
  // Grid position points in r and z
  std::vector<double> rPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Acts::Vector2D> bField;
  // reserve estimated size
  rPos.reserve(nPoints);
  zPos.reserve(nPoints);
  bField.reserve(nPoints);
  // [1] Read in file and fill values
  std::ifstream map_file(fieldMapFile.c_str(), std::ios::in);
  std::string line;
  double r = 0., z = 0.;
  double br = 0., bz = 0.;
  while (std::getline(map_file, line)) {
    if (line.empty() || line[0] == '%' || line[0] == '#' ||
        line.find_first_not_of(' ') == std::string::npos)
      continue;

    std::istringstream tmp(line);
    tmp >> r >> z >> br >> bz;
    rPos.push_back(r);
    zPos.push_back(z);
    bField.push_back(Acts::Vector2D(br, bz));
  }
  map_file.close();
  /// [2] use helper function in core
  return Acts::fieldMapperRZ(localToGlobalBin, rPos, zPos, bField, lengthUnit,
                             BFieldUnit, firstQuadrant);
}

Acts::InterpolatedBFieldMapper<Acts::detail::Grid<
    Acts::Vector3D, Acts::detail::EquidistantAxis,
    Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>
FW::BField::txt::fieldMapperXYZ(
    std::function<size_t(std::array<size_t, 3> binsXYZ,
                         std::array<size_t, 3> nBinsXYZ)>
        localToGlobalBin,
    std::string fieldMapFile, double lengthUnit, double BFieldUnit,
    size_t nPoints, bool firstOctant) {
  /// [1] Read in field map file
  // Grid position points in x, y and z
  std::vector<double> xPos;
  std::vector<double> yPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Acts::Vector3D> bField;
  // reserve estimated size
  xPos.reserve(nPoints);
  yPos.reserve(nPoints);
  zPos.reserve(nPoints);
  bField.reserve(nPoints);
  // [1] Read in file and fill values
  std::ifstream map_file(fieldMapFile.c_str(), std::ios::in);
  std::string line;
  double x = 0., y = 0., z = 0.;
  double bx = 0., by = 0., bz = 0.;
  while (std::getline(map_file, line)) {
    if (line.empty() || line[0] == '%' || line[0] == '#' ||
        line.find_first_not_of(' ') == std::string::npos)
      continue;

    std::istringstream tmp(line);
    tmp >> x >> y >> z >> bx >> by >> bz;
    xPos.push_back(x);
    yPos.push_back(y);
    zPos.push_back(z);
    bField.push_back(Acts::Vector3D(bx, by, bz));
  }
  map_file.close();

  return Acts::fieldMapperXYZ(localToGlobalBin, xPos, yPos, zPos, bField,
                              lengthUnit, BFieldUnit, firstOctant);
}

Acts::InterpolatedBFieldMapper<
    Acts::detail::Grid<Acts::Vector2D, Acts::detail::EquidistantAxis,
                       Acts::detail::EquidistantAxis>>
FW::BField::root::fieldMapperRZ(
    std::function<size_t(std::array<size_t, 2> binsRZ,
                         std::array<size_t, 2> nBinsRZ)>
        localToGlobalBin,
    std::string fieldMapFile, std::string treeName, double lengthUnit,
    double BFieldUnit, bool firstQuadrant) {
  /// [1] Read in field map file
  // Grid position points in r and z
  std::vector<double> rPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Acts::Vector2D> bField;
  // [1] Read in file and fill values
  TFile* inputFile = TFile::Open(fieldMapFile.c_str());
  TTree* tree = (TTree*)inputFile->Get(treeName.c_str());
  Int_t entries = tree->GetEntries();

  double r, z;
  double Br, Bz;

  tree->SetBranchAddress("r", &r);
  tree->SetBranchAddress("z", &z);

  tree->SetBranchAddress("Br", &Br);
  tree->SetBranchAddress("Bz", &Bz);

  // reserve size
  rPos.reserve(entries);
  zPos.reserve(entries);
  bField.reserve(entries);

  for (int i = 0; i < entries; i++) {
    tree->GetEvent(i);
    rPos.push_back(r);
    zPos.push_back(z);
    bField.push_back(Acts::Vector2D(Br, Bz));
  }
  inputFile->Close();
  /// [2] use helper function in core
  return Acts::fieldMapperRZ(localToGlobalBin, rPos, zPos, bField, lengthUnit,
                             BFieldUnit, firstQuadrant);
}

Acts::InterpolatedBFieldMapper<Acts::detail::Grid<
    Acts::Vector3D, Acts::detail::EquidistantAxis,
    Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>
FW::BField::root::fieldMapperXYZ(
    std::function<size_t(std::array<size_t, 3> binsXYZ,
                         std::array<size_t, 3> nBinsXYZ)>
        localToGlobalBin,
    std::string fieldMapFile, std::string treeName, double lengthUnit,
    double BFieldUnit, bool firstOctant) {
  /// [1] Read in field map file
  // Grid position points in x, y and z
  std::vector<double> xPos;
  std::vector<double> yPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Acts::Vector3D> bField;
  // [1] Read in file and fill values
  TFile* inputFile = TFile::Open(fieldMapFile.c_str());
  TTree* tree = (TTree*)inputFile->Get(treeName.c_str());
  Int_t entries = tree->GetEntries();

  double x, y, z;
  double Bx, By, Bz;

  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);

  tree->SetBranchAddress("Bx", &Bx);
  tree->SetBranchAddress("By", &By);
  tree->SetBranchAddress("Bz", &Bz);

  // reserve size
  xPos.reserve(entries);
  yPos.reserve(entries);
  zPos.reserve(entries);
  bField.reserve(entries);

  for (int i = 0; i < entries; i++) {
    tree->GetEvent(i);
    xPos.push_back(x);
    yPos.push_back(y);
    zPos.push_back(z);
    bField.push_back(Acts::Vector3D(Bx, By, Bz));
  }
  inputFile->Close();

  return Acts::fieldMapperXYZ(localToGlobalBin, xPos, yPos, zPos, bField,
                              lengthUnit, BFieldUnit, firstOctant);
}
