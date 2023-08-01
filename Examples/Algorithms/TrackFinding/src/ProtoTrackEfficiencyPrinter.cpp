#include "ActsExamples/TrackFindingX/ProtoTrackEfficiencyPrinter.hpp"

#include <Acts/Utilities/Zip.hpp>

#include <boost/histogram.hpp>
#include <boost/histogram/ostream.hpp>

#include <list>

namespace {
  template<typename T>
  std::list<T> vecToList(const std::vector<T> &v) {
    return std::list<T>(v.begin(), v.end());
  }
}

ActsExamples::ProcessCode ActsExamples::ProtoTrackEfficiencyPrinter::execute(
    const ActsExamples::AlgorithmContext &context) const {
  auto testTracks = vecToList(m_testProtoTracks(context));
  auto refTracks = m_refProtoTracks(context);

  std::for_each(testTracks.begin(), testTracks.end(),
                [](auto &t) { std::sort(t.begin(), t.end()); });
  std::for_each(refTracks.begin(), refTracks.end(),
                [](auto &t) { std::sort(t.begin(), t.end()); });

  ACTS_INFO("Receiving " << refTracks.size() << " reference tracks");
  refTracks.erase(std::remove_if(refTracks.begin(), refTracks.end(),
                                 [](const auto &t) { return t.size() < 3; }),
                  refTracks.end());
  ACTS_INFO(" -> " << refTracks.size() << " tracks with size >= 3");
  ACTS_INFO("Receiving " << testTracks.size() << " test tracks");

  std::vector<double> effs(refTracks.size(), 0.0);

  for (auto [refTrack, eff] : Acts::zip(refTracks, effs)) {
    ProtoTrack intersection;
    std::optional<typename decltype(testTracks)::iterator> toDelete;

    for (auto testTrackIt = testTracks.begin(); testTrackIt != testTracks.end(); ++testTrackIt) {
      const auto &testTrack = *testTrackIt;
      intersection.resize(std::max(testTrack.size(), refTrack.size()));
      const auto it = std::set_intersection(refTrack.begin(), refTrack.end(), testTrack.begin(),
                            testTrack.end(), intersection.begin());
      const auto size = static_cast<double>(std::distance(intersection.begin(), it));
      intersection.clear();

      eff = std::max(eff, size / refTrack.size());

      if( eff > 0.50001) {
        toDelete = testTrackIt;
        break;
      }
    }

    if( toDelete ) {
      testTracks.erase(*toDelete);
    }
  }

  namespace bh = boost::histogram;
  auto h = bh::make_histogram(bh::axis::regular<>(11, 0.0, 1.1));
  h.fill(effs);

  ACTS_INFO("Prototrack efficiency histogram:\n" << h);

  return {};
}
