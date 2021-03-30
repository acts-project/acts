// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/ContextualDetector/PrealignedDetector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/PrealignedDetectorElement.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorOptions.hpp"

#include <boost/program_options.hpp>

void PrealignedDetector::addOptions(
    boost::program_options::options_description& opt) const {
  // Add the generic geometry options
  ActsExamples::Options::addGenericGeometryOptions(opt);
  // specify the rotation setp
  opt.add_options()(
      "prealign-seed",
      boost::program_options::value<size_t>()->default_value(1324354657),
      "Seed for the decorator random numbers.")(
      "prealign-iovs",
      boost::program_options::value<size_t>()->default_value(25),
      "Number of a valid IOVs.")(
      "prealign-sigma-iplane",
      boost::program_options::value<double>()->default_value(100.),
      "Sigma of the in-plane misalignment in [um]")(
      "prealign-sigma-oplane",
      boost::program_options::value<double>()->default_value(50.),
      "Sigma of the out-of-plane misalignment in [um]")(
      "prealign-sigma-irot",
      boost::program_options::value<double>()->default_value(20.),
      "Sigma of the in-plane rotation misalignment in [mrad]")(
      "prealign-sigma-orot",
      boost::program_options::value<double>()->default_value(0.),
      "Sigma of the out-of-plane rotation misalignment in [mrad]")(
      "prealign-loglevel",
      boost::program_options::value<size_t>()->default_value(3),
      "Output log level of the alignment decorator.")(
      "prealign-firstnominal",
      boost::program_options::value<bool>()->default_value(false),
      "Keep the first iov batch nominal.");
}

auto PrealignedDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;

  auto buildLevel = vm["geo-generic-buildlevel"].template as<size_t>();
  // set geometry building logging level
  Acts::Logging::Level surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  bool buildProto =
      (vm["mat-input-type"].template as<std::string>() == "proto");

  /// return the generic detector - with aligned context decorator
  TrackingGeometryPtr aTrackingGeometry =
      ActsExamples::Generic::buildDetector<DetectorElement>(
          nominalContext, detectorStore, buildLevel, std::move(mdecorator),
          buildProto, surfaceLogLevel, layerLogLevel, volumeLogLevel);

  Acts::Logging::Level decoratorLogLevel =
      Acts::Logging::Level(vm["prealign-loglevel"].template as<size_t>());

  // Let's create a reandom number service
  ActsExamples::RandomNumbers::Config randomNumberConfig;
  randomNumberConfig.seed = vm["prealign-seed"].template as<size_t>();
  auto randomNumberSvc =
      std::make_shared<ActsExamples::RandomNumbers>(randomNumberConfig);

  // Alignment decorator service
  Decorator::Config pagcsConfig;
  pagcsConfig.detectorStore = detectorStore;
  pagcsConfig.nIovs = vm["prealign-iovs"].template as<size_t>();
  pagcsConfig.randomNumberSvc = randomNumberSvc;

  // Now create the alignment decorator
  ContextDecorators aContextDecorators = {std::make_shared<Decorator>(
      pagcsConfig,
      Acts::getDefaultLogger("PrealignedDecorator", decoratorLogLevel))};

  bool firstIovNominal = vm["prealign-firstnominal"].template as<bool>();

  // Create an local generator for the preContext
  size_t max = std::numeric_limits<size_t>::max();
  ActsExamples::WhiteBoard preStore(Acts::getDefaultLogger(
      "PreEventStore#" + std::to_string(max), decoratorLogLevel));
  ActsExamples::AlgorithmContext preContext(max, max, preStore);
  ActsExamples::RandomEngine rng = randomNumberSvc->spawnGenerator(preContext);
  std::normal_distribution<double> gauss(0., 1.);

  // The misalingments
  double sigmaIp = vm["prealign-sigma-iplane"].template as<double>();
  double sigmaOp = vm["prealign-sigma-oplane"].template as<double>();
  double sigmaIr = vm["prealign-sigma-irot"].template as<double>();
  double sigmaOr = vm["prealign-sigma-orot"].template as<double>();
  double gSigmaX = sigmaIp * Acts::UnitConstants::um;
  double gSigmaY = sigmaIp * Acts::UnitConstants::um;
  double gSigmaZ = sigmaOp * Acts::UnitConstants::um;
  double aSigmaX = sigmaOr * 0.001;  // millirad
  double aSigmaY = sigmaOr * 0.001;  // millirad
  double aSigmaZ = sigmaIr * 0.001;  // millirad

  for (unsigned int iov = 0; iov < pagcsConfig.nIovs; ++iov) {
    // Loop over the detector store and create the different alignments
    for (auto& lstore : detectorStore) {
      for (auto& ldet : lstore) {
        // Get the nominal transform
        auto& tForm = ldet->nominalTransform(preContext.geoContext);
        // Create a new transform
        auto atForm = std::make_unique<Acts::Transform3>(tForm);

        if (iov != 0 or not firstIovNominal) {
          // The shifts in x, y, z
          double tx = gSigmaX != 0 ? gSigmaX * gauss(rng) : 0.;
          double ty = gSigmaY != 0 ? gSigmaY * gauss(rng) : 0.;
          double tz = gSigmaZ != 0 ? gSigmaZ * gauss(rng) : 0.;
          // Add a translation - if there is any
          if (tx != 0. or ty != 0. or tz != 0.) {
            const auto& tMatrix = atForm->matrix();
            auto colX = tMatrix.block<3, 1>(0, 0).transpose();
            auto colY = tMatrix.block<3, 1>(0, 1).transpose();
            auto colZ = tMatrix.block<3, 1>(0, 2).transpose();
            Acts::Vector3 newCenter = tMatrix.block<3, 1>(0, 3).transpose() +
                                      tx * colX + ty * colY + tz * colZ;
            atForm->translation() = newCenter;
          }
          // Modify it - rotation around local X
          if (aSigmaX != 0.) {
            (*atForm) *=
                Acts::AngleAxis3(aSigmaX * gauss(rng), Acts::Vector3::UnitX());
          }
          if (aSigmaY != 0.) {
            (*atForm) *=
                Acts::AngleAxis3(aSigmaY * gauss(rng), Acts::Vector3::UnitY());
          }
          if (aSigmaZ != 0.) {
            (*atForm) *=
                Acts::AngleAxis3(aSigmaZ * gauss(rng), Acts::Vector3::UnitZ());
          }
        }
        // Add a new Alignment transform to the store
        ldet->addAlignedTransform(std::move(atForm), iov);
      }
    }
  }

  // return the pair of geometry and the alignment decorator(s)
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(aTrackingGeometry), std::move(aContextDecorators));
}
