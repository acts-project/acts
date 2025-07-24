// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Module/Entries.hpp"
#include "ActsPython/Utilities/Context.hpp"

/// This adds the python module entries for the Examples
namespace ActsPython {
    // These will be added to `acts.examples` module
    void addFramework(Context& ctx);
    void addDetector(Context& ctx);
    void addGenerators(Context& ctx);
    void addAmbiguityResolution(Context& ctx);
    void addMagneticFieldMaps(Context& ctx);
    void addDigitizationAlgorithms(Context& ctx);
    void addMaterialMappingAlgorithms(Context& ctx);
    void addPropagationAlgorithms(Context& ctx);
    void addTrackFindingAlgorithms(Context& ctx);
    void addTrackFittingAlgorithms(Context& ctx);
    void addTruthTrackingAlgorithms(Context& ctx);
    void addVertexingAlgorithms(Context& ctx);
    void addPrinterAlgorithms(Context& ctx);
    void addUtilityAlgorithms(Context& ctx);
    void addInput(Context& ctx);
    void addOutput(Context& ctx);
    // Plugin dependent components
    void addGeoModelDetector(Context& ctx);
    void addTGeoDetector(Context& ctx);
    void addRootInput(Context& ctx);
    void addRootOutput(Context& ctx);
    void addJsonInputOutput(Context& ctx);
    void addSvgOutput(Context& ctx);
    void addFatrasAlgorithms(Context& ctx);
    void addTracccAlgorithms(Context& ctx);
    void addPythia8Generator(Context& ctx);
    void addTruthJetAlgorithms(Context& ctx);
    void addOnnxAlgorithms(Context& ctx);

}  // namespace ActsPython

void ActsPython::addExamplesModule(Context& ctx) {
    addFramework(ctx);
    addDetector(ctx);
    addGenerators(ctx);
    addAmbiguityResolution(ctx);
    addDigitizationAlgorithms(ctx);
    addPropagationAlgorithms(ctx);
    addMagneticFieldMaps(ctx);
    addMaterialMappingAlgorithms(ctx);
    addTrackFindingAlgorithms(ctx);
    addTrackFittingAlgorithms(ctx);
    addTruthTrackingAlgorithms(ctx);
    addVertexingAlgorithms(ctx);
    addPrinterAlgorithms(ctx);
    addUtilityAlgorithms(ctx);
    addInput(ctx);
    addOutput(ctx);
    // Plugin dependent components
    addGeoModelDetector(ctx);
    addTGeoDetector(ctx);
    addRootInput(ctx);
    addRootOutput(ctx);
    addJsonInputOutput(ctx);
    addSvgOutput(ctx);
    addFatrasAlgorithms(ctx);
    addTracccAlgorithms(ctx);
    addPythia8Generator(ctx);
    addTruthJetAlgorithms(ctx);
    addOnnxAlgorithms(ctx);
}