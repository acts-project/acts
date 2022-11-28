#include "ActsExamples/Detector/GenericDetectorWithOptions.hpp"

#include "ActsExamples/Options/GenericDetectorOptions.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {

void GenericDetectorWithOptions::addOptions(
    boost::program_options::options_description& opt) const {
  ActsExamples::Options::addGenericGeometryOptions(opt);
}

auto GenericDetectorWithOptions::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  GenericDetector::Config cfg;

  cfg.buildLevel = vm["geo-generic-buildlevel"].template as<size_t>();
  // set geometry building logging level
  cfg.surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  cfg.layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  cfg.volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  cfg.buildProto = (vm["mat-input-type"].template as<std::string>() == "proto");

  return m_detector.finalize(cfg, std::move(mdecorator));
}

}  // namespace ActsExamples
