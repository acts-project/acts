// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/CosmicSeedCuts.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Utilities/StripModulePairing.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

/// Create space point representations from measurements.
///
/// This implements simple space point construction, where each surface-based
/// measurement translates into one space point using the surface
/// local-to-global transform.
///
/// The algorithm takes both the source links and measurements container as
/// input. The source link container is geometry-sorted and each element is
/// small compared to a measurement. The geometry selection is therefore much
/// easier to perform on the source links than on the unsorted measurements.
///
/// There are no explicit requirements on the content of the input measurements.
/// If no local positions are measured, the transformed global positions will
/// always be the position of the module origin.
///
/// Modo cósmico (opcional, todo apagado por defecto): con
/// cosmicStripSpacePoints los space points de strips se forman sin vértice, y
/// outputSpacePointClasses / outputSpacePointSurfaces /
/// outputSpacePointParticlesMap dan las entradas que usan
/// CosmicTripletSeedingAlgorithm y la medición de eficiencia. Es la
/// contraparte con geometría ACTS (p. ej. el ODD) de SimHitSpacePointMaker.
class SpacePointMaker final : public IAlgorithm {
 public:
  struct Config {
    /// Input measurement subset (initial full subset from digitization for pass
    /// 1; filtered subset from MeasurementFilterAlgorithm for subsequent
    /// passes).
    std::string inputMeasurements;
    /// Output space points collection.
    std::string outputSpacePoints;
    /// Tracking geometry for transformation lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// For which part of the detector geometry should space points be created.
    ///
    /// Only volumes and layers can be set. Zero values can be used as wildcards
    /// to select larger parts of the hierarchy, i.e. setting only the volume
    /// selects all measurements within that volume. Adding a single identifier
    /// with all components set to zero selects all available measurements. The
    /// selection must not have duplicates.
    std::vector<Acts::GeometryIdentifier> geometrySelection;
    /// Geometry selection for strip modules
    std::vector<Acts::GeometryIdentifier> stripGeometrySelection;
    /// Assumed vertex position for the strip space point formation.
    Acts::Vector3 stripVertex = Acts::Vector3::Zero();
    /// Tolerance scaling factor on the strip detector element length.
    double stripLengthTolerance = 0.01;
    /// Tolerance on the gap between the two strip detector elements, in mm.
    double stripLengthGapTolerance = 0.01;

    /// Cósmicos: forma los space points de strips con la máxima aproximación
    /// entre las dos strips (StripSpacePointBuilder::computeCosmicSpacePoint)
    /// en lugar de la formación restringida al vértice. Un cósmico no viene
    /// de stripVertex, así que stripVertex y stripLengthGapTolerance no se
    /// usan (stripLengthTolerance sí). Sin vértice tampoco hay una dirección
    /// con la que elegir el mejor compañero estéreo: cada par de medidas de
    /// los dos módulos cuyo cruce cae dentro de ambas strips da un space
    /// point.
    bool cosmicStripSpacePoints = false;
    /// sin^2 mínimo del ángulo entre las dos strips en el modo cósmico
    /// (StripSpacePointBuilder::CosmicOptions::tolerance).
    double cosmicStripTolerance = 1e-6;
    /// Modo cósmico: guarda todos los pairs estéreo, aunque la máxima
    /// aproximación caiga fuera de las strips (sin el corte de
    /// stripLengthTolerance), con el space point en la máxima aproximación
    /// llevada al extremo más cercano de la strip. Sin la dirección de la
    /// traza no se puede saber dónde la cruza; esa posición es sólo un punto
    /// de partida que CosmicTripletSeedFinder recalcula con la dirección de
    /// cada triplete (CosmicSeedCuts::useStripCalibration), y es ahí donde se
    /// descartan los pairs que no son de una traza.
    bool cosmicStripKeepAllPairs = false;

    /// Opcional: si no está vacío, escribe un
    /// std::vector<Acts::CosmicSpacePointClass> con la clase de cada space
    /// point (mismo índice que el space point). Es la entrada
    /// inputSpacePointClasses de CosmicTripletSeedingAlgorithm.
    std::string outputSpacePointClasses;
    /// Clase por geometryId().volume() de la (primera) medida del space
    /// point; have falta al menos una lista con outputSpacePointClasses. Un
    /// mismo volumen no puede estar en dos listas, y una medida de un volumen
    /// que no está en ninguna aborta el evento (listando los que faltan).
    std::vector<std::uint32_t> pixelBarrelVolumes;
    std::vector<std::uint32_t> pixelEndcapVolumes;
    std::vector<std::uint32_t> stripBarrelVolumes;
    std::vector<std::uint32_t> stripEndcapVolumes;

    /// Opcional: si no está vacío, escribe un
    /// std::vector<Acts::GeometryIdentifier> con la superficie de la primera
    /// medida de cada space point (mismo índice que el space point), como
    /// sur() de Athena, que es la del primer cluster. Es la entrada
    /// inputSpacePointSurfaces de CosmicTripletSeedingAlgorithm. Dentro de
    /// un evento, todos los space points de un par de módulos estéreo tienen
    /// la primera medida en el mismo módulo del par, así que el id identifica
    /// al módulo de dos lados.
    std::string outputSpacePointSurfaces;

    /// Opcional, para medir eficiencia/pureza de los seeds: con los dos
    /// configurados escribe un IndexMultimap<SimBarcode> de índice de space
    /// point a las partículas de sus medidas (sin repetidas; un space point
    /// de strips con dos partículas distintas aparece con las dos).
    std::string inputMeasurementParticlesMap;
    std::string outputSpacePointParticlesMap;
  };

  /// Construct the space point maker.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  explicit SpacePointMaker(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the space point construction.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  StripModulePairMap m_stripModulePairMap;

  std::optional<IndexSourceLink::SurfaceAccessor> m_slSurfaceAccessor;

  ReadDataHandle<MeasurementSubset> m_inputMeasurements{this,
                                                        "InputMeasurements"};

  WriteDataHandle<SpacePointContainer> m_outputSpacePoints{this,
                                                           "OutputSpacePoints"};

  WriteDataHandle<std::vector<Acts::CosmicSpacePointClass>>
      m_outputSpacePointClasses{this, "OutputSpacePointClasses"};
  WriteDataHandle<std::vector<Acts::GeometryIdentifier>>
      m_outputSpacePointSurfaces{this, "OutputSpacePointSurfaces"};

  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  WriteDataHandle<IndexMultimap<SimBarcode>> m_outputSpacePointParticlesMap{
      this, "OutputSpacePointParticlesMap"};
};

}  // namespace ActsExamples
