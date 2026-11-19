// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <fstream>

#include <gtest/gtest.h>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/detectors/io.hpp"
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/framework/types.hpp"

using algebra = detray::test::algebra;
using metadata = detray::registry::telescope_metadata<algebra>;

TEST(detray_registered_io, round_trip_and_explicit_override) {
  vecmem::host_memory_resource mr;
  detray::tel_det_config<algebra> config{20.f, 20.f};
  auto [concrete, names] =
      detray::build_telescope_detector<algebra>(mr, config);
  auto n_surfaces = concrete.surfaces().size();
  auto n_volumes = concrete.volumes().size();
  detray::host_detector detector;
  detector.set<metadata>(std::move(concrete));
  detray::io::detector_writer_config writer;
  writer.replace_files(true);
  detray::io::write_registered_detector<algebra>(detector, names, writer);
  const auto filename = detector.as<metadata>().name(names) + "_geometry.json";
  nlohmann::ordered_json json;
  {
    std::ifstream stream{filename};
    stream >> json;
  }
  EXPECT_EQ(json["header"]["common"]["metadata"], "telescope");
  detray::io::detector_reader_config reader;
  reader.add_file(filename);
  auto [restored, restored_names] =
      detray::io::read_registered_detector_json<algebra>(mr, reader);
  ASSERT_TRUE(restored.is<metadata>());
  EXPECT_EQ(restored.as<metadata>().surfaces().size(), n_surfaces);
  EXPECT_EQ(restored.as<metadata>().volumes().size(), n_volumes);
  EXPECT_EQ(restored.as<metadata>().name(restored_names),
            detector.as<metadata>().name(names));

  json["header"]["common"]["metadata"] = "not_registered";
  {
    std::ofstream stream{filename};
    stream << json;
  }
  EXPECT_THROW(detray::io::read_registered_detector_json<algebra>(mr, reader),
               std::invalid_argument);
  auto overridden = detray::io::read_registered_detector_json<algebra>(
      mr, reader, "telescope");
  EXPECT_TRUE(overridden.first.is<metadata>());
}

TEST(detray_registered_io, metadata_identity) {
  EXPECT_EQ(detray::registry::metadata_id<metadata>(), "telescope");
  EXPECT_EQ(detray::registry::resolve_metadata_id("odd_metadata"), "odd");
  EXPECT_THROW(detray::registry::resolve_metadata_id("unknown"),
               std::invalid_argument);
  EXPECT_THROW(detray::registry::resolve_metadata_id(""),
               std::invalid_argument);
}
