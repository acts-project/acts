// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/EventCsvWriter.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SeedingTruth.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/TrackingGeometryLayout.hpp"

#include <limits>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace ActsPython {

void addFatrasSynthetic(py::module_& fatras) {
  using namespace ActsFatras::Synthetic;

  py::module_ m = fatras.def_submodule(
      "synthetic",
      "Synthetic space point event generation on a simplified, "
      "RZ-symmetric detector layout");

  py::enum_<SurfaceShape>(m, "SurfaceShape")
      .value("Cylinder", SurfaceShape::Cylinder)
      .value("Disc", SurfaceShape::Disc);

  py::enum_<SurfaceSide>(m, "SurfaceSide")
      .value("Negative", SurfaceSide::Negative)
      .value("Barrel", SurfaceSide::Barrel)
      .value("Positive", SurfaceSide::Positive);

  // Neither is bound by the core module, and a layout that carries its
  // detector's material has to be able to say what that material is.
  py::class_<Acts::Material>(m, "Material")
      .def_static(
          "fromMolarDensity",
          [](float x0, float l0, float ar, float z, float molarRho) {
            return Acts::Material::fromMolarDensity(x0, l0, ar, z, molarRho);
          },
          py::arg("x0"), py::arg("l0"), py::arg("ar"), py::arg("z"),
          py::arg("molarRho"))
      .def_static(
          "fromMassDensity",
          [](float x0, float l0, float ar, float z, float massRho) {
            return Acts::Material::fromMassDensity(x0, l0, ar, z, massRho);
          },
          py::arg("x0"), py::arg("l0"), py::arg("ar"), py::arg("z"),
          py::arg("massRho"))
      .def_property_readonly("X0", &Acts::Material::X0)
      .def_property_readonly("L0", &Acts::Material::L0)
      .def_property_readonly("Ar", &Acts::Material::Ar)
      .def_property_readonly("Z", &Acts::Material::Z)
      .def_property_readonly("molarDensity", &Acts::Material::molarDensity)
      .def_property_readonly("massDensity", &Acts::Material::massDensity)
      .def_property_readonly("meanExcitationEnergy",
                             &Acts::Material::meanExcitationEnergy);

  py::class_<Acts::MaterialSlab>(m, "MaterialSlab")
      .def(py::init<const Acts::Material&, float>(), py::arg("material"),
           py::arg("thickness"))
      .def_property_readonly("material", &Acts::MaterialSlab::material)
      .def_property_readonly("thickness", &Acts::MaterialSlab::thickness)
      .def_property_readonly("thicknessInX0",
                             &Acts::MaterialSlab::thicknessInX0)
      .def_property_readonly("thicknessInL0",
                             &Acts::MaterialSlab::thicknessInL0)
      .def("scaleThickness", &Acts::MaterialSlab::scaleThickness,
           py::arg("scale"))
      .def_static("combineLayers",
                  py::overload_cast<const std::vector<Acts::MaterialSlab>&>(
                      &Acts::MaterialSlab::combineLayers),
                  py::arg("layers"));

  m.def("materialSlab", &materialSlab, py::arg("x0"), py::arg("l0"),
        py::arg("ar"), py::arg("z"), py::arg("molarDensity"),
        py::arg("thickness"));

  py::class_<BandComposition>(m, "BandComposition")
      .def(py::init(
               [](float ar, float z, float molarDensityX0, float thickness) {
                 return BandComposition{ar, z, molarDensityX0, thickness};
               }),
           py::arg("ar"), py::arg("z"), py::arg("molarDensityX0"),
           py::arg("thickness"))
      .def_readwrite("ar", &BandComposition::ar)
      .def_readwrite("z", &BandComposition::z)
      .def_readwrite("molarDensityX0", &BandComposition::molarDensityX0)
      .def_readwrite("thickness", &BandComposition::thickness);

  py::class_<SurfaceMaterial>(m, "SurfaceMaterial")
      .def(py::init<>())
      .def(py::init([](std::vector<float> bounds,
                       std::vector<Acts::MaterialSlab> bands) {
             return SurfaceMaterial{std::move(bounds), std::move(bands)};
           }),
           py::arg("bounds"), py::arg("bands"))
      .def(py::init([](const BandComposition& composition,
                       std::vector<float> bounds,
                       const std::vector<float>& radiationLengths,
                       const std::vector<float>& nuclearLengths) {
             return SurfaceMaterial{composition, std::move(bounds),
                                    radiationLengths, nuclearLengths};
           }),
           py::arg("composition"), py::arg("bounds"),
           py::arg("radiationLengths"), py::arg("nuclearLengths"))
      .def(py::init(
               [](Acts::MaterialSlab slab) { return SurfaceMaterial{slab}; }),
           py::arg("slab"))
      .def("average", &SurfaceMaterial::average)
      .def_readwrite("bounds", &SurfaceMaterial::bounds)
      .def_readonly("bands", &SurfaceMaterial::bands)
      .def("at", &SurfaceMaterial::at, py::arg("along"),
           py::return_value_policy::reference_internal);

  py::enum_<EnergyLossModel>(m, "EnergyLossModel")
      .value("Mean", EnergyLossModel::Mean)
      .value("Mode", EnergyLossModel::Mode);

  py::class_<DetectorLayer>(m, "DetectorLayer")
      .def_readonly("shape", &DetectorLayer::shape)
      .def_readonly("refCoord", &DetectorLayer::refCoord)
      .def_readonly("minBound", &DetectorLayer::minBound)
      .def_readonly("maxBound", &DetectorLayer::maxBound)
      .def_readonly("side", &DetectorLayer::side)
      .def_readonly("layer", &DetectorLayer::layer)
      .def_readonly("moduleIndex", &DetectorLayer::moduleIndex);

  // what `GeneratedParticle::productionSurface` says for a particle that came
  // off no surface at all
  m.attr("kNoSurface") = kNoSurface;

  py::class_<DetectorSurface>(m, "DetectorSurface")
      .def_readonly("shape", &DetectorSurface::shape)
      .def_readonly("refCoord", &DetectorSurface::refCoord)
      .def_readonly("layers", &DetectorSurface::layers)
      .def_readwrite("material", &DetectorSurface::material)
      .def_readonly("minBound", &DetectorSurface::minBound)
      .def_readonly("maxBound", &DetectorSurface::maxBound)
      .def_readwrite("overlapProbability", &DetectorSurface::overlapProbability)
      .def_readwrite("overlapOffset", &DetectorSurface::overlapOffset)
      .def_property_readonly(
          "passive", [](const DetectorSurface& s) { return s.layers.empty(); });

  py::class_<DetectorLayout>(m, "DetectorLayout")
      .def(py::init<>())
      .def_readonly("layers", &DetectorLayout::layers)
      .def_readonly("surfaces", &DetectorLayout::surfaces)
      .def_readwrite("escapeRadius", &DetectorLayout::escapeRadius)
      .def_readwrite("escapeHalfZ", &DetectorLayout::escapeHalfZ);

  // The `add...` methods return the builder so that calls can be chained. That
  // has to be handed back as a reference to the same object: the default policy
  // for an lvalue reference is to copy, which would quietly add everything
  // after the first call to a copy that is then thrown away.
  py::class_<DetectorLayoutBuilder>(m, "DetectorLayoutBuilder")
      .def(py::init<>())
      .def("addPassiveCylinder", &DetectorLayoutBuilder::addPassiveCylinder,
           py::arg("radius"),
           py::arg("maxAbsZ") = std::numeric_limits<float>::infinity(),
           py::arg("minAbsZ") = 0.f,
           py::return_value_policy::reference_internal)
      .def("addPassiveDisc", &DetectorLayoutBuilder::addPassiveDisc,
           py::arg("side"), py::arg("absZ"), py::arg("rMin") = 0.f,
           py::arg("rMax") = std::numeric_limits<float>::infinity(),
           py::return_value_policy::reference_internal)
      .def("addCylinder", &DetectorLayoutBuilder::addCylinder,
           py::arg("radius"), py::arg("halfLengthZ"), py::arg("numModules"),
           py::return_value_policy::reference_internal)
      // the span overload cannot be bound directly: pybind11 has no caster for
      // `std::span`, so the rings are taken as a vector and converted here
      .def(
          "addDisc",
          [](DetectorLayoutBuilder& builder, SurfaceSide side, float absZ,
             const std::vector<RingBounds>& rings) -> DetectorLayoutBuilder& {
            return builder.addDisc(side, absZ, rings);
          },
          py::arg("side"), py::arg("absZ"), py::arg("rings"),
          py::return_value_policy::reference_internal)
      .def("addDisc",
           py::overload_cast<SurfaceSide, float, float, float, std::uint32_t>(
               &DetectorLayoutBuilder::addDisc),
           py::arg("side"), py::arg("absZ"), py::arg("rMin"), py::arg("rMax"),
           py::arg("numRings"), py::return_value_policy::reference_internal)
      .def("build", &DetectorLayoutBuilder::build);

  py::class_<RingBounds>(m, "RingBounds")
      .def(py::init<>())
      .def(py::init(
               [](float rMin, float rMax) { return RingBounds{rMin, rMax}; }),
           py::arg("rMin"), py::arg("rMax"))
      .def_readwrite("rMin", &RingBounds::rMin)
      .def_readwrite("rMax", &RingBounds::rMax);

  py::class_<DiscDescription>(m, "DiscDescription")
      .def(py::init<>())
      .def_readwrite("absZ", &DiscDescription::absZ)
      .def_readwrite("material", &DiscDescription::material)
      .def_readwrite("rings", &DiscDescription::rings);

  py::class_<PassiveSurfaceDescription>(m, "PassiveSurfaceDescription")
      .def(py::init<>())
      .def_readwrite("shape", &PassiveSurfaceDescription::shape)
      .def_readwrite("refCoord", &PassiveSurfaceDescription::refCoord)
      .def_readwrite("minBound", &PassiveSurfaceDescription::minBound)
      .def_readwrite("maxBound", &PassiveSurfaceDescription::maxBound)
      .def_readwrite("material", &PassiveSurfaceDescription::material);

  m.def("uniformRings", &uniformRings, py::arg("rMin"), py::arg("rMax"),
        py::arg("numRings"));
  m.def(
      "subdivideRings",
      [](const std::vector<RingBounds>& rings, std::uint32_t parts) {
        return subdivideRings(rings, parts);
      },
      py::arg("rings"), py::arg("parts"));

  py::class_<BarrelEndcapDescription>(m, "BarrelEndcapDescription")
      .def(py::init<>())
      .def_readwrite("beamPipeRadius", &BarrelEndcapDescription::beamPipeRadius)
      .def_readwrite("beamPipeMaterial",
                     &BarrelEndcapDescription::beamPipeMaterial)
      .def_readwrite("passiveSurfaces",
                     &BarrelEndcapDescription::passiveSurfaces)
      .def_readwrite("barrelRadii", &BarrelEndcapDescription::barrelRadii)
      .def_readwrite("barrelHalfLengthsZ",
                     &BarrelEndcapDescription::barrelHalfLengthsZ)
      .def_readwrite("barrelModules", &BarrelEndcapDescription::barrelModules)
      .def_readwrite("discs", &BarrelEndcapDescription::discs)
      .def_readwrite("barrelMaterials",
                     &BarrelEndcapDescription::barrelMaterials)
      .def_readwrite("barrelOverlapProbabilities",
                     &BarrelEndcapDescription::barrelOverlapProbabilities)
      .def_readwrite("discOverlapProbability",
                     &BarrelEndcapDescription::discOverlapProbability)
      .def_readwrite("barrelOverlapOffset",
                     &BarrelEndcapDescription::barrelOverlapOffset)
      .def_readwrite("discOverlapOffset",
                     &BarrelEndcapDescription::discOverlapOffset)
      .def_readwrite("escapeRadius", &BarrelEndcapDescription::escapeRadius)
      .def_readwrite("escapeHalfZ", &BarrelEndcapDescription::escapeHalfZ);

  m.def("makeLayout", &makeLayout, py::arg("description"));
  m.def("itkPixelDescription", &itkPixelDescription);
  m.def("makeItkPixelLayout", &makeItkPixelLayout);
  m.def("openDataDetectorPixelDescription", &openDataDetectorPixelDescription);
  m.def("makeOpenDataDetectorPixelLayout", &makeOpenDataDetectorPixelLayout);
  m.def("genericDetectorPixelDescription", &genericDetectorPixelDescription);
  m.def("makeGenericDetectorPixelLayout", &makeGenericDetectorPixelLayout);

  py::class_<TrackingGeometryLayoutOptions>(m, "TrackingGeometryLayoutOptions")
      .def(py::init<>())
      // The selector takes a pointer rather than a reference: `Acts::Surface`
      // is not copyable, and pybind11 would try to copy it to hand a reference
      // argument back to a Python callable.
      .def(
          "setSurfaceSelector",
          [](TrackingGeometryLayoutOptions& options,
             const std::function<bool(const Acts::Surface*)>& selector) {
            if (!selector) {
              options.surfaceSelector = nullptr;
              return;
            }
            options.surfaceSelector = [selector](const Acts::Surface& surface) {
              return selector(&surface);
            };
          },
          py::arg("selector"))
      .def_readwrite("includeMaterialSurfaces",
                     &TrackingGeometryLayoutOptions::includeMaterialSurfaces)
      .def_readwrite("materialFromGeometry",
                     &TrackingGeometryLayoutOptions::materialFromGeometry)
      .def_readwrite("materialBands",
                     &TrackingGeometryLayoutOptions::materialBands)
      .def_readwrite("materialBandTolerance",
                     &TrackingGeometryLayoutOptions::materialBandTolerance)
      .def_readwrite("passiveBeamPipeRadius",
                     &TrackingGeometryLayoutOptions::passiveBeamPipeRadius)
      .def_readwrite("modulesPerSurface",
                     &TrackingGeometryLayoutOptions::modulesPerSurface)
      .def_readwrite("cylinderRTolerance",
                     &TrackingGeometryLayoutOptions::cylinderRTolerance)
      .def_readwrite("discZTolerance",
                     &TrackingGeometryLayoutOptions::discZTolerance)
      .def_readwrite("maxRingOverlap",
                     &TrackingGeometryLayoutOptions::maxRingOverlap)
      .def_readwrite("ringRTolerance",
                     &TrackingGeometryLayoutOptions::ringRTolerance);

  m.def("makeLayoutFromTrackingGeometry", &makeLayoutFromTrackingGeometry,
        py::arg("trackingGeometry"), py::arg("geometryContext"),
        py::arg("options") = TrackingGeometryLayoutOptions{});

  py::class_<GenerationConfig>(m, "GenerationConfig")
      .def(py::init<>())
      .def_readwrite("pileup", &GenerationConfig::pileup)
      .def_readwrite("chargedPerUnitRapidity",
                     &GenerationConfig::chargedPerUnitRapidity)
      .def_readwrite("minPt", &GenerationConfig::minPt)
      .def_readwrite("ptScale", &GenerationConfig::ptScale)
      .def_readwrite("ptExponent", &GenerationConfig::ptExponent)
      .def_readwrite("minRapidity", &GenerationConfig::minRapidity)
      .def_readwrite("maxRapidity", &GenerationConfig::maxRapidity)
      .def_readwrite("rapidityEdge", &GenerationConfig::rapidityEdge)
      .def_readwrite("rapidityEdgeWidth", &GenerationConfig::rapidityEdgeWidth)
      .def_readwrite("beamspotSigmaZ", &GenerationConfig::beamspotSigmaZ)
      .def_readwrite("d0Sigma", &GenerationConfig::d0Sigma)
      .def("numPrimaries", &GenerationConfig::numPrimaries);

  py::class_<PropagationConfig>(m, "PropagationConfig")
      .def(py::init<>())
      .def_readwrite("maxTurns", &PropagationConfig::maxTurns)
      .def_readwrite("tracksPerPrimary", &PropagationConfig::tracksPerPrimary)
      .def_readwrite("hitsPerPrimary", &PropagationConfig::hitsPerPrimary)
      .def_readwrite("maxTracksPerPrimary",
                     &PropagationConfig::maxTracksPerPrimary);

  py::class_<MaterialConfig>(m, "MaterialConfig")
      .def(py::init<>())
      .def_readwrite("maxDiscPathLength", &MaterialConfig::maxDiscPathLength)
      .def_readwrite("maxCylinderPathLength",
                     &MaterialConfig::maxCylinderPathLength)
      .def_readwrite("scale", &MaterialConfig::scale)
      .def_readwrite("multipleScattering", &MaterialConfig::multipleScattering)
      .def_readwrite("energyLoss", &MaterialConfig::energyLoss)
      .def_readwrite("energyLossModel", &MaterialConfig::energyLossModel)
      .def_readwrite("maxEnergyLossFraction",
                     &MaterialConfig::maxEnergyLossFraction);

  py::class_<MeasurementConfig>(m, "MeasurementConfig")
      .def(py::init<>())
      .def_readwrite("positionSmearing", &MeasurementConfig::positionSmearing)
      .def_readwrite("overlapScale", &MeasurementConfig::overlapScale);

  py::class_<SecondarySamplingConfig>(m, "SecondarySamplingConfig")
      .def(py::init<>())
      .def_readwrite("minPt", &SecondarySamplingConfig::minPt)
      .def_readwrite("electronScale", &SecondarySamplingConfig::electronScale)
      .def_readwrite("electronExponent",
                     &SecondarySamplingConfig::electronExponent)
      .def_readwrite("electronSpread", &SecondarySamplingConfig::electronSpread)
      .def_readwrite("electronKt", &SecondarySamplingConfig::electronKt)
      .def_readwrite("momentumScale", &SecondarySamplingConfig::momentumScale)
      .def_readwrite("momentumExponent",
                     &SecondarySamplingConfig::momentumExponent)
      .def_readwrite("momentumSpread", &SecondarySamplingConfig::momentumSpread)
      .def_readwrite("kt", &SecondarySamplingConfig::kt)
      .def_readwrite("evaporationFraction",
                     &SecondarySamplingConfig::evaporationFraction)
      .def_readwrite("evaporationScale",
                     &SecondarySamplingConfig::evaporationScale);

  py::class_<SecondaryConfig>(m, "SecondaryConfig")
      .def(py::init<>())
      .def_readwrite("electronRate", &SecondaryConfig::electronRate)
      .def_readwrite("nuclearRate", &SecondaryConfig::nuclearRate)
      .def_readwrite("decayYield", &SecondaryConfig::decayYield)
      .def_readwrite("decayLength", &SecondaryConfig::decayLength)
      .def_readwrite("stubRate", &SecondaryConfig::stubRate)
      .def_readwrite("stubClusters", &SecondaryConfig::stubClusters)
      .def_readwrite("stubReach", &SecondaryConfig::stubReach)
      .def_readwrite("maxGenerations", &SecondaryConfig::maxGenerations)
      .def_readwrite("maxPerCrossing", &SecondaryConfig::maxPerCrossing)
      .def_readwrite("sampling", &SecondaryConfig::sampling);

  py::class_<SimulationConfig>(m, "SimulationConfig")
      .def(py::init<>())
      .def_readwrite("propagation", &SimulationConfig::propagation)
      .def_readwrite("material", &SimulationConfig::material)
      .def_readwrite("measurement", &SimulationConfig::measurement)
      .def_readwrite("secondaries", &SimulationConfig::secondaries);

  py::class_<EventConfig>(m, "EventConfig")
      .def(py::init<>())
      .def_readwrite("generation", &EventConfig::generation)
      .def_readwrite("simulation", &EventConfig::simulation)
      .def_readwrite("particlePdg", &EventConfig::particlePdg)
      .def_readwrite("bFieldZ", &EventConfig::bFieldZ)
      .def_readwrite("seed", &EventConfig::seed)
      .def("particleMass", &EventConfig::particleMass)
      .def_static("itkPixelTtbarPu200", &EventConfig::itkPixelTtbarPu200)
      .def_static("openDataDetectorTtbarPu200",
                  &EventConfig::openDataDetectorTtbarPu200);

  py::class_<GeneratedParticle>(m, "GeneratedParticle")
      .def_readonly("pt", &GeneratedParticle::pt)
      .def_readonly("eta", &GeneratedParticle::eta)
      .def_readonly("phi", &GeneratedParticle::phi)
      .def_readonly("d0", &GeneratedParticle::d0)
      .def_readonly("z0", &GeneratedParticle::z0)
      .def_readonly("charge", &GeneratedParticle::charge)
      .def_readonly("productionRadius", &GeneratedParticle::productionRadius)
      .def_readonly("productionZ", &GeneratedParticle::productionZ)
      .def_readonly("productionSurface", &GeneratedParticle::productionSurface)
      .def_readonly("generation", &GeneratedParticle::generation)
      .def_property_readonly("primary", &GeneratedParticle::primary)
      .def_readonly("numHits", &GeneratedParticle::numHits);

  // `Acts::SpacePointContainer` exposes only its fixed columns to Python, so
  // the two dynamic columns that carry the generator truth are pulled out here.
  // Without them the space points cannot be related back to the layout or to
  // the particles that made them.
  const auto column = [](const Event& event, const std::string& name) {
    const auto proxy = event.spacePoints.column<std::uint32_t>(name);
    std::vector<std::uint32_t> values;
    values.reserve(event.spacePoints.size());
    for (const auto sp : event.spacePoints) {
      values.push_back(sp.extra(proxy));
    }
    return values;
  };

  py::class_<Event>(m, "Event")
      .def(py::init<>())
      .def_readonly("spacePoints", &Event::spacePoints)
      .def_readonly("particles", &Event::particles)
      .def_property_readonly(
          "layerIds",
          [column](const Event& event) { return column(event, "layerId"); },
          "Index into `DetectorLayout.layers` of each space point's layer")
      .def_property_readonly(
          "particleIds",
          [column](const Event& event) { return column(event, "particleId"); },
          "Index into `Event.particles` of each space point's particle");

  py::class_<EventGenerator>(m, "EventGenerator")
      .def(py::init<const DetectorLayout&, const EventConfig&>(),
           py::arg("layout"), py::arg("config"),
           // the generator only keeps a pointer to the layout
           py::keep_alive<1, 2>())
      .def("generate", py::overload_cast<>(&EventGenerator::generate))
      .def("reset", &EventGenerator::reset, py::arg("seed"));

  m.def("generateEvent", &generateEvent, py::arg("layout"), py::arg("config"));

  py::class_<EventSummary>(m, "EventSummary")
      .def_readonly("spacePoints", &EventSummary::spacePoints)
      .def_readonly("primaries", &EventSummary::primaries)
      .def_readonly("secondaries", &EventSummary::secondaries)
      .def_readonly("seedablePrimaries", &EventSummary::seedablePrimaries)
      .def_readonly("primaryHits", &EventSummary::primaryHits)
      .def_readonly("secondaryHits", &EventSummary::secondaryHits);

  m.def("summarize", &summarize, py::arg("event"), py::arg("ptThreshold"));

  m.def("writeEventCsv", &writeEventCsv, py::arg("event"), py::arg("layout"),
        py::arg("prefix"));
}

}  // namespace ActsPython
