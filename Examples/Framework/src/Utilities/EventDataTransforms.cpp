#include "ActsExamples/Utilities/EventDataTransforms.hpp"

ActsExamples::ProtoTrack ActsExamples::seedToPrototrack(const ActsExamples::SimSeed &seed) {
  ProtoTrack track;
  track.reserve(seed.sp().size());
  for (auto spacePointPtr : seed.sp()) {
    for (const auto& slink : spacePointPtr->sourceLinks()) {
      const auto& islink = slink.get<IndexSourceLink>();
      track.emplace_back(islink.index());
    }
  }
  return track;
}

ActsExamples::ProtoTrackContainer ActsExamples::seedsToPrototracks(const ActsExamples::SimSeedContainer &seeds) {
  ProtoTrackContainer tracks;
  tracks.reserve(seeds.size());

  for (const auto& seed : seeds) {
    tracks.push_back(seedToPrototrack(seed));
  }

  return tracks;
}
