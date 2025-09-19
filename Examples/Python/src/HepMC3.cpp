// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3InputConverter.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3OutputConverter.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {
void addHepMC3(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  auto hepmc3 = mex.def_submodule("_hepmc3");

  ACTS_PYTHON_DECLARE_WRITER(HepMC3Writer, hepmc3, "HepMC3Writer", outputPath,
                             perEvent, inputEvent, compression,
                             maxEventsPending, writeEventsInOrder);

  ACTS_PYTHON_DECLARE_READER(HepMC3Reader, hepmc3, "HepMC3Reader", inputPath,
                             inputPaths, outputEvent, printListing, numEvents,
                             checkEventNumber, maxEventBufferSize,
                             vertexGenerator, randomNumbers);

  ACTS_PYTHON_DECLARE_ALGORITHM(HepMC3OutputConverter, hepmc3,
                                "HepMC3OutputConverter", inputParticles,
                                inputVertices, outputEvent);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      HepMC3InputConverter, hepmc3, "HepMC3InputConverter", inputEvent,
      outputParticles, outputVertices, printListing, checkConsistency,
      mergePrimaries, primaryVertexSpatialThreshold, vertexSpatialThreshold,
      mergeSecondaries);

  {
    using enum HepMC3Util::Compression;
    py::enum_<HepMC3Util::Compression>(hepmc3, "Compression")
        .value("none", none)
        .value("zlib", zlib)
        .value("lzma", lzma)
        .value("bzip2", bzip2)
        .value("zstd", zstd);
  }

  hepmc3.def("availableCompressionModes", []() {
    auto modes = HepMC3Util::availableCompressionModes();
    return std::vector(modes.begin(), modes.end());
  });

  hepmc3.def("compressionExtension", &HepMC3Util::compressionExtension);
}
}  // namespace ActsPython
