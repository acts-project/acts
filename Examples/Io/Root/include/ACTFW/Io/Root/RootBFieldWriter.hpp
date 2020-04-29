// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Units.hpp>
#include <TFile.h>
#include <TTree.h>
#include <array>
#include <boost/optional.hpp>
#include <ios>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"

namespace FW {

/// @class RootBFieldWriter
///
/// Writes out the Acts::InterpolatedbFieldMap. Currently implemented for 'rz'
/// and 'xyz' field maps.
template <typename bfield_t>
class RootBFieldWriter {
 public:
  /// Describes the axes definition of the grid of the magnetic field map.
  enum class GridType { rz = 0, xyz = 1 };

  struct Config {
    /// The name of the output tree
    std::string treeName = "TTree";
    /// The name of the output file
    std::string fileName = "TFile.root";
    /// the file access mode (recreate by default)
    std::string fileMode = "recreate";
    /// The magnetic field to be written out
    std::shared_ptr<const bfield_t> bField = nullptr;
    /// How the magnetic field map should be written out
    GridType gridType = GridType::xyz;
    /// [optional] Setting the range to be printed out in either r (for
    /// cylinder coordinates) or x/y (in cartesian coordinates)
    /// @note setting this parameter is optional, in case no boundaries are
    /// handed over the full magnetic field map will be printed out
    boost::optional<std::array<double, 2>> rBounds;
    /// [optional] Setting the range in z to be printed out
    /// @note setting this parameter is optional, in case no boundaries are
    /// handed over the full magnetic field map will be printed out
    boost::optional<std::array<double, 2>> zBounds;
    /// Number of bins in r
    /// @note setting this parameter is optional, in case no bin numbers are
    /// handed over the full magnetic field map will be printed out
    size_t rBins = 200;
    /// Number of bins in z
    // @note setting this parameter is optional, in case no bin numbers are
    /// handed over the full magnetic field map will be printed out
    size_t zBins = 300;
    /// Number of bins in phi
    // @note setting this parameter is optional, in case no bin numbers are
    /// handed over the full magnetic field map will be printed out
    size_t phiBins = 100;
  };

  /// Write down an interpolated magnetic field map
  static void run(const Config& cfg,
                  std::unique_ptr<const Acts::Logger> p_logger =
                      Acts::getDefaultLogger("RootBFieldWriter",
                                             Acts::Logging::INFO)) {
    // Set up (local) logging
    // @todo Remove dangerous using declaration once the logger macro
    // tolerates it
    using namespace Acts;
    ACTS_LOCAL_LOGGER(std::move(p_logger))

    // Check basic configuration
    if (cfg.treeName.empty()) {
      throw std::invalid_argument("Missing tree name");
    } else if (cfg.fileName.empty()) {
      throw std::invalid_argument("Missing file name");
    } else if (!cfg.bField) {
      throw std::invalid_argument("Missing interpolated magnetic field");
    }

    // Setup ROOT I/O
    ACTS_INFO("Registering new ROOT output File : " << cfg.fileName);
    TFile* outputFile = TFile::Open(cfg.fileName.c_str(), cfg.fileMode.c_str());
    if (!outputFile) {
      throw std::ios_base::failure("Could not open '" + cfg.fileName);
    }
    outputFile->cd();
    TTree* outputTree = new TTree(cfg.treeName.c_str(), cfg.treeName.c_str());
    if (!outputTree)
      throw std::bad_alloc();

    // The position values
    double z;
    outputTree->Branch("z", &z);

    // The BField values
    double Bz;
    outputTree->Branch("Bz", &Bz);

    // Get the underlying mapper of the InterpolatedBFieldMap
    auto mapper = cfg.bField->getMapper();

    // Access the minima and maxima of all axes
    auto minima = mapper.getMin();
    auto maxima = mapper.getMax();
    auto nBins = mapper.getNBins();

    if (cfg.gridType == GridType::xyz) {
      ACTS_INFO("Map will be written out in cartesian coordinates (x,y,z).");

      // Write out the interpolated magnetic field map
      double stepX = 0., stepY = 0., stepZ = 0.;
      double minX = 0., minY = 0., minZ = 0.;
      double maxX = 0., maxY = 0., maxZ = 0.;
      size_t nBinsX = 0, nBinsY = 0, nBinsZ = 0;

      // The position values in xy
      double x;
      outputTree->Branch("x", &x);
      double y;
      outputTree->Branch("y", &y);
      // The BField values in xy
      double Bx;
      outputTree->Branch("Bx", &Bx);
      double By;
      outputTree->Branch("By", &By);

      // check if range is user defined
      if (cfg.rBounds && cfg.zBounds) {
        ACTS_INFO("User defined ranges handed over.");

        // print out map in user defined range
        minX = cfg.rBounds->at(0);
        minY = cfg.rBounds->at(0);
        minZ = cfg.zBounds->at(0);

        maxX = cfg.rBounds->at(1);
        maxY = cfg.rBounds->at(1);
        maxZ = cfg.zBounds->at(1);

        nBinsX = cfg.rBins;
        nBinsY = cfg.rBins;
        nBinsZ = cfg.zBins;

      } else {
        ACTS_INFO(
            "No user defined ranges handed over - printing out whole map.");
        // print out whole map
        // check dimension of Bfieldmap
        if (minima.size() == 3 && maxima.size() == 3) {
          minX = minima.at(0);
          minY = minima.at(1);
          minZ = minima.at(2);

          maxX = maxima.at(0);
          maxY = maxima.at(1);
          maxZ = maxima.at(2);

          nBinsX = nBins.at(0);
          nBinsY = nBins.at(1);
          nBinsZ = nBins.at(2);

        } else if (minima.size() == 2 && maxima.size() == 2) {
          minX = -maxima.at(0);
          minY = -maxima.at(0);
          minZ = minima.at(1);

          maxX = maxima.at(0);
          maxY = maxima.at(0);
          maxZ = maxima.at(1);

          nBinsX = nBins.at(0);
          nBinsY = nBins.at(0);
          nBinsZ = nBins.at(1);
        } else {
          std::ostringstream errorMsg;
          errorMsg
              << "BField has wrong dimension. The dimension needs to be "
                 "either 2 (r,z,Br,Bz) or 3(x,y,z,Bx,By,Bz) in order to be "
                 "written out by this writer.";
          throw std::invalid_argument(errorMsg.str());
        }
      }

      stepX = fabs(minX - maxX) / nBinsX;
      stepY = fabs(minY - maxY) / nBinsY;
      stepZ = fabs(minZ - maxZ) / nBinsZ;

      for (size_t i = 0; i <= nBinsX; i++) {
        double raw_x = minX + i * stepX;
        for (size_t j = 0; j <= nBinsY; j++) {
          double raw_y = minY + j * stepY;
          for (size_t k = 0; k <= nBinsZ; k++) {
            double raw_z = minZ + k * stepZ;
            Acts::Vector3D position(raw_x, raw_y, raw_z);
            if (cfg.bField->isInside(position)) {
              auto bField = cfg.bField->getField(position);

              x = raw_x / Acts::units::_mm;
              y = raw_y / Acts::units::_mm;
              z = raw_z / Acts::units::_mm;
              Bx = bField.x() / Acts::units::_T;
              By = bField.y() / Acts::units::_T;
              Bz = bField.z() / Acts::units::_T;
              outputTree->Fill();
            }
          }  // for z
        }    // for y
      }      // for x
    } else {
      ACTS_INFO("Map will be written out in cylinder coordinates (r,z).");
      // The position value in r
      double r;
      outputTree->Branch("r", &r);
      // The BField value in r
      double Br;
      outputTree->Branch("Br", &Br);

      double minR = 0, maxR = 0;
      double minZ = 0, maxZ = 0;
      size_t nBinsR = 0, nBinsZ = 0, nBinsPhi = 0;
      double stepR = 0, stepZ = 0;

      if (cfg.rBounds && cfg.zBounds) {
        ACTS_INFO("User defined ranges handed over.");

        minR = cfg.rBounds->at(0);
        minZ = cfg.zBounds->at(0);

        maxR = cfg.rBounds->at(1);
        maxZ = cfg.zBounds->at(1);

        nBinsR = cfg.rBins;
        nBinsZ = cfg.zBins;
        nBinsPhi = cfg.phiBins;
      } else {
        ACTS_INFO(
            "No user defined ranges handed over - printing out whole map.");

        if (minima.size() == 3 && maxima.size() == 3) {
          minR = 0.;
          minZ = minima.at(2);

          maxR = maxima.at(0);
          maxZ = maxima.at(2);

          nBinsR = nBins.at(0);
          nBinsZ = nBins.at(2);
          nBinsPhi = 100.;

        } else if (minima.size() == 2 || maxima.size() == 2) {
          minR = minima.at(0);
          minZ = minima.at(1);

          maxR = maxima.at(0);
          maxZ = maxima.at(1);

          nBinsR = nBins.at(0);
          nBinsZ = nBins.at(1);
          nBinsPhi = 100.;

        } else {
          std::ostringstream errorMsg;
          errorMsg
              << "BField has wrong dimension. The dimension needs to be "
                 "either 2 (r,z,Br,Bz) or 3(x,y,z,Bx,By,Bz) in order to be "
                 "written out by this writer.";
          throw std::invalid_argument(errorMsg.str());
        }
      }
      double minPhi = -M_PI;
      stepR = fabs(minR - maxR) / nBinsR;
      stepZ = fabs(minZ - maxZ) / nBinsZ;
      double stepPhi = (2 * M_PI) / nBinsPhi;

      for (size_t i = 0; i < nBinsPhi; i++) {
        double phi = minPhi + i * stepPhi;
        for (size_t k = 0; k < nBinsZ; k++) {
          double raw_z = minZ + k * stepZ;
          for (size_t j = 0; j < nBinsR; j++) {
            double raw_r = minR + j * stepR;
            Acts::Vector3D position(raw_r * cos(phi), raw_r * sin(phi), raw_z);
            if (cfg.bField->isInside(position)) {
              auto bField = cfg.bField->getField(position);
              z = raw_z / Acts::units::_mm;
              r = raw_r / Acts::units::_mm;
              Bz = bField.z() / Acts::units::_T;
              Br = VectorHelpers::perp(bField) / Acts::units::_T;
              outputTree->Fill();
            }
          }
        }
      }  // for
    }

    // Tear down ROOT I/O
    ACTS_INFO("Closing and Writing ROOT output File : " << cfg.fileName);
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
  }
};

}  // namespace FW
