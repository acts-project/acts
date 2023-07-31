#include "ActsExamples/TrackFindingX/ProtoTrackEfficiencyPrinter.hpp"

#include <Acts/Utilities/Zip.hpp>

#include <boost/histogram.hpp>
#include <boost/histogram/ostream.hpp>

ActsExamples::ProcessCode ActsExamples::ProtoTrackEfficiencyPrinter::execute(
    const ActsExamples::AlgorithmContext &context) const {
  auto testTracks = m_testProtoTracks(context);
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
    for (const auto &testTrack : testTracks) {
      std::set_intersection(refTrack.begin(), refTrack.end(), testTrack.begin(),
                            testTrack.end(), intersection.begin());
      eff = std::max(
          eff, static_cast<double>(intersection.size()) / refTrack.size());
      intersection.clear();
    }
  }

  for (auto eff : effs) {
    std::cout << eff << " ";
  }
  std::cout << std::endl;

  namespace bh = boost::histogram;
  auto h = bh::make_histogram(bh::axis::regular<>(10, 0.0, 1.0));
  h.fill(effs);

  ACTS_INFO("Prototrack efficiency histogram:\n" << h);

  return {};
}
