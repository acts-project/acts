#include "ActsExamples/Detray/ActsToDetrayDetectorAlg.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"

#include <detray/builders/detector_builder.hpp>
#include <detray/io/backend/geometry_reader.hpp>
#include <detray/io/backend/homogeneous_material_reader.hpp>
#include <detray/io/backend/material_map_reader.hpp>
#include <detray/io/backend/surface_grid_reader.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/utils/consistency_checker.hpp>

#include <vecmem/memory/host_memory_resource.hpp>
#include <detray/io/json/json.hpp>
#include <fstream>
#include <filesystem>

namespace ActsExamples {

ActsToDetrayDetectorAlg::ActsToDetrayDetectorAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ActsToDetrayDetectorAlg", std::move(logger)), m_cfg(cfg) {
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument(
        "ActsToDetrayDetectorAlg: trackingGeometry is null");
  }
  m_outputDetrayToActsMap.initialize(m_cfg.outputDetrayToActsMap);
}

ProcessCode ActsToDetrayDetectorAlg::initialize() {

  Acts::GeometryContext gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();

  // ── Convert TrackingGeometry → detray payloads ──────────────────────────
  ActsPlugins::DetrayPayloadConverter::Config convCfg;

  ACTS_INFO("Looking for beampipe volume: " << m_cfg.beampipeVolumeName);

  m_cfg.trackingGeometry->apply(
      [&convCfg, this](const Acts::TrackingVolume& volume) {
          ACTS_DEBUG("Found volume: " << volume.volumeName());
          if (volume.volumeName() == m_cfg.beampipeVolumeName) {
              convCfg.beampipeVolume = &volume;
              ACTS_INFO("Found beampipe volume: " << volume.volumeName());
          }
      });

  // Find beampipe volume
  m_cfg.trackingGeometry->apply(
      [&convCfg, this](const Acts::TrackingVolume& volume) {
        if (volume.volumeName() == m_cfg.beampipeVolumeName) {
          convCfg.beampipeVolume = &volume;
        }
      });

  if (convCfg.beampipeVolume == nullptr) {
    ACTS_WARNING("DetrayGeometryProvider: beampipe volume '"
                 << m_cfg.beampipeVolumeName << "' not found");
  }

  ActsPlugins::DetrayPayloadConverter converter(
      convCfg, logger().clone("DetrayPayloadConverter"));

  auto payloads = converter.convertTrackingGeometry(gctx, *m_cfg.trackingGeometry);

  // ── Build detray detector from payloads ──────────────────────────────────
  using detector_t = detray::detector<detray::odd_metadata<detray::array<double>>>;

  vecmem::host_memory_resource mr;
  detray::detector_builder<detector_t::metadata> detectorBuilder{};

  detray::io::geometry_reader::from_payload<detector_t>(
      detectorBuilder, *payloads.detector);

  detray::io::homogeneous_material_reader::from_payload<detector_t>(
      detectorBuilder, *payloads.homogeneousMaterial);

  detray::io::material_map_reader<std::integral_constant<std::size_t, 2>>::
      from_payload<detector_t>(detectorBuilder,
                               std::move(*payloads.materialGrids));

  detray::io::surface_grid_reader<
    typename detector_t::surface_type,
    std::integral_constant<std::size_t, 0>,
    std::integral_constant<std::size_t, 2>>::
    template from_payload<detector_t>(detectorBuilder, *payloads.surfaceGrids);

  detector_t detrayDetector(detectorBuilder.build(mr));

  detray::detail::check_consistency(detrayDetector);

  ACTS_INFO("DetrayGeometryProvider: built detray detector with "
            << detrayDetector.volumes().size() << " volumes and "
            << detrayDetector.surfaces().size() << " surfaces");

  // ── Optionally write JSON files ──────────────────────────────────────────
  if (!m_cfg.outputJsonDir.empty()) {
    auto toDetrayNameMap =
        [](const std::map<unsigned int, std::string>& src) {
          detray::name_map result;
          for (const auto& [idx, name] : src) {
            result.emplace(static_cast<detray::dindex>(idx), name);
          }
          return result;
        };

    auto detrayNames = toDetrayNameMap(payloads.names);
    ACTS_INFO("Detector name: " << (payloads.names.empty() ? "EMPTY" : payloads.names.at(0)));
    auto writer_cfg = detray::io::detector_writer_config{}
                      .format(detray::io::format::json)
                      .replace_files(true)
                      .path(m_cfg.outputJsonDir);
    // detrayNames[0] = "odd";
    // std::cerr << "detrayNames size=" << detrayNames.size() << "\n";
    // for (const auto& [k, v] : detrayNames) {
    //     std::cerr << "  [" << k << "] = " << v << "\n";
    // }

    // This one took me such a long time to find
    // Traccc uses the JSON header to decide what detector metadata type
    // to use and it's basically:
    // "Cylindrical detector from DD4hep blueprint" for ODD
    // "detray_detector" for ITk
    // anything else is default_metadata
    // detrayNames.set_detector_name("odd");

    detrayNames.set_detector_name("odd");
    detray::io::write_detector(detrayDetector, detrayNames, writer_cfg);

    namespace fs = std::filesystem;
    fs::path geoFile = fs::path(m_cfg.outputJsonDir) / "odd_geometry.json";

    // After writing:
    nlohmann::json j;
    {
        std::ifstream ifs(geoFile);
        if (!ifs.is_open()) {
            throw std::runtime_error("Cannot open " + geoFile.string());
        }
        ifs >> j;
    }

    ACTS_INFO("Before patch: " << j["header"]["common"]["detector"].get<std::string>());
    j["header"]["common"]["detector"] = "Cylindrical detector from DD4hep blueprint";
    ACTS_INFO("After patch: " << j["header"]["common"]["detector"].get<std::string>());

    {
        std::ofstream ofs(geoFile, std::ios::trunc);
        if (!ofs.is_open()) {
            throw std::runtime_error("Cannot write " + geoFile.string());
        }
        ofs << std::setw(2) << j;
        ofs.flush();
    }

    // Verify
    {
        std::ifstream verify(geoFile);
        nlohmann::json j2;
        verify >> j2;
        ACTS_INFO("Verified: " << j2["header"]["common"]["detector"].get<std::string>());
    }
  }

  // ── Build detray→Acts geometry ID map ───────────────────────────────────
  std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> detrayToActsMap;

  for (const auto& surface : detrayDetector.surfaces()) {
    // surface.source is the Acts GeometryIdentifier encoded as uint64
    const Acts::GeometryIdentifier actsId(surface.source);
    if (actsId.sensitive() == 0) {
      continue;  // skip portals and passives
    }
    detrayToActsMap[surface.identifier().value()] = actsId;
  }

  ACTS_INFO("DetrayGeometryProvider: built detray→Acts map with "
            << detrayToActsMap.size() << " sensitive surfaces");

  m_detrayToActsMap = std::move(detrayToActsMap);
  return ProcessCode::SUCCESS;
}

ProcessCode ActsToDetrayDetectorAlg::execute(const AlgorithmContext& ctx) const {
    // Just write map to whiteboard each event
    m_outputDetrayToActsMap(ctx,
        std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>(m_detrayToActsMap));
    return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples