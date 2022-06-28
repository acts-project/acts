#include "ActsExamples/TrackFinding/TrajectoriesToPrototracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

namespace ActsExamples {
  
ProcessCode TrajectoriesToPrototracks::execute(
    const AlgorithmContext& ctx) const {
  const auto trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  ProtoTrackContainer tracks;

  for (const auto& trajectory : trajectories) {
    for (const auto tip : trajectory.tips()) {
      ProtoTrack track;

      trajectory.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
        if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        const auto& source_link =
            static_cast<const IndexSourceLink&>(state.uncalibrated());
        track.push_back(source_link.index());

        return true;
      });

      tracks.push_back(track);
    }
  }

  ctx.eventStore.add(m_cfg.outputPrototracks, std::move(tracks));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
