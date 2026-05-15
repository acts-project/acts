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
#include "ActsExamples/Utilities/MultiplicityGenerators.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsHepMC3, hepmc3) {
  ACTS_PYTHON_DECLARE_WRITER(HepMC3Writer, hepmc3, "HepMC3Writer", outputPath,
                             inputEvent, compression, maxEventsPending,
                             writeEventsInOrder);

  // Declare the HepMC3Reader class first
  auto reader =
      py::class_<HepMC3Reader, IReader, std::shared_ptr<HepMC3Reader>>(
          hepmc3, "HepMC3Reader")
          .def(py::init<const HepMC3Reader::Config&, Acts::Logging::Level>(),
               py::arg("config"), py::arg("level"))
          .def_property_readonly("config", &HepMC3Reader::config);

  // Expose Input struct as a nested class of HepMC3Reader
  py::class_<HepMC3Reader::Input>(reader, "Input")
      .def(py::init<>())
      .def(py::init([](const std::filesystem::path& path,
                       const std::shared_ptr<const MultiplicityGenerator>&
                           multiplicityGenerator) {
             HepMC3Reader::Input inp;
             inp.path = path;
             inp.multiplicityGenerator = multiplicityGenerator;
             return inp;
           }),
           py::arg("path"), py::arg("multiplicityGenerator"))
      .def_readwrite("path", &HepMC3Reader::Input::path)
      .def_readwrite("multiplicityGenerator",
                     &HepMC3Reader::Input::multiplicityGenerator)
      // Factory methods for convenience
      .def_static(
          "Fixed",
          [](const std::filesystem::path& path, std::size_t n) {
            HepMC3Reader::Input inp;
            inp.path = path;
            inp.multiplicityGenerator =
                std::make_shared<FixedMultiplicityGenerator>(n);
            return inp;
          },
          py::arg("path"), py::arg("n") = 1,
          "Create Input with FixedMultiplicityGenerator")
      .def_static(
          "Poisson",
          [](const std::filesystem::path& path, double mean) {
            HepMC3Reader::Input inp;
            inp.path = path;
            inp.multiplicityGenerator =
                std::make_shared<PoissonMultiplicityGenerator>(mean);
            return inp;
          },
          py::arg("path"), py::arg("mean"),
          "Create Input with PoissonMultiplicityGenerator");

  auto config = py::class_<HepMC3Reader::Config>(reader, "Config")
                    .def(py::init<>(), "Default constructor");
  // Now configure the HepMC3Reader itself
  ACTS_PYTHON_STRUCT(config, inputs, inputPath, outputEvent, printListing,
                     numEvents, checkEventNumber, maxEventBufferSize,
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

  {
    using enum HepMC3Util::Format;
    py::enum_<HepMC3Util::Format>(hepmc3, "Format")
        .value("ascii", ascii)
        .value("root", root);
  }

  hepmc3.def("availableCompressionModes", []() {
    auto modes = HepMC3Util::availableCompressionModes();
    return std::vector(modes.begin(), modes.end());
  });

  hepmc3.def("availableFormats", []() {
    auto formats = HepMC3Util::availableFormats();
    return std::vector(formats.begin(), formats.end());
  });

  hepmc3.def("compressionExtension", &HepMC3Util::compressionExtension);
  hepmc3.def("compressionFromFilename", &HepMC3Util::compressionFromFilename,
             py::arg("filename"));
  hepmc3.def("formatFromFilename", &HepMC3Util::formatFromFilename,
             py::arg("filename"));

  // HepMC3 normalize function and result
  {
    auto result =
        py::class_<HepMC3Util::NormalizeResult>(hepmc3, "NormalizeResult")
            .def(py::init<>())
            .def_readonly("numEvents", &HepMC3Util::NormalizeResult::numEvents)
            .def_readonly("outputFiles",
                          &HepMC3Util::NormalizeResult::outputFiles)
            .def_readonly("totalInputSize",
                          &HepMC3Util::NormalizeResult::totalInputSize)
            .def_readonly("totalOutputSize",
                          &HepMC3Util::NormalizeResult::totalOutputSize)
            .def_readonly("totalReadTime",
                          &HepMC3Util::NormalizeResult::totalReadTime)
            .def_readonly("totalWriteTime",
                          &HepMC3Util::NormalizeResult::totalWriteTime);
  }

  hepmc3.def("normalizeFiles", &HepMC3Util::normalizeFiles,
             py::arg("inputFiles"), py::arg("singleOutputPath") = std::nullopt,
             py::arg("outputDir") = ".", py::arg("outputPrefix") = "events",
             py::arg("eventsPerFile") = 10000,
             py::arg("maxEvents") = std::nullopt,
             py::arg("format") = HepMC3Util::Format::ascii,
             py::arg("compression") = HepMC3Util::Compression::none,
             py::arg("compressionLevel") = 6, py::arg("verbose") = false);
}
