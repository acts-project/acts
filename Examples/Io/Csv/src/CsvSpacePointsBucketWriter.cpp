// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSpacePointsBucketWriter.hpp"

#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <cstddef>

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvSpacePointsBucketWriter::CsvSpacePointsBucketWriter(
    const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputBuckets, "CsvSpacePointsBucketWriter", level),
      m_cfg(config) {}

CsvSpacePointsBucketWriter::~CsvSpacePointsBucketWriter() = default;

ProcessCode CsvSpacePointsBucketWriter::finalize() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ProcessCode CsvSpacePointsBucketWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<SimSpacePointContainer>& buckets) {
  // Open per-event file for all components
  std::string pathBucket =
      perEventFilepath(m_cfg.outputDir, "buckets.csv", ctx.eventNumber);

  NamedTupleCsvWriter<SpacePointBucketData> writerBucket(pathBucket,
                                                         m_cfg.outputPrecision);

  SpacePointBucketData bucketData{};

  int bucketIdx = 0;
  for (const auto& bucket : buckets) {
    if (bucket.empty()) {
      continue;
    }
    // Split the bucket into lines of 20 space points to manage variable sizes
    auto numLines =
        static_cast<int>((bucket.size() / 20) +
                         static_cast<std::size_t>(bucket.size() % 20 != 0));
    for (int nLines = 0; nLines < numLines; nLines++) {
      bucketData.bucketIdx = bucketIdx;
      bucketData.bucketSize = bucket.size();
      int maxSPIdx =
          std::min(20, static_cast<int>(bucket.size() - nLines * 20));
      for (int SPIdx = 0; SPIdx < maxSPIdx; SPIdx++) {
        bucketData.measurement_id[SPIdx] = (bucket[nLines * 20 + SPIdx])
                                               .sourceLinks()[0]
                                               .get<IndexSourceLink>()
                                               .index();
      }
      writerBucket.append(bucketData);
    }
    bucketIdx++;
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
