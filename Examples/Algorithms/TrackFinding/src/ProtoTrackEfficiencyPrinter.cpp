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
      intersection.resize(std::max(testTrack.size(), refTrack.size()));
      // std::cout << "\nref: ";
      // std::copy(refTrack.begin(), refTrack.end(), std::ostream_iterator<std::size_t>(std::cout, ","));
      // std::cout << "\ntest: ";
      // std::copy(testTrack.begin(), testTrack.end(), std::ostream_iterator<std::size_t>(std::cout, ","));
      // std::cout << "\n-----------------------" << std::endl;

      const auto it = std::set_intersection(refTrack.begin(), refTrack.end(), testTrack.begin(),
                            testTrack.end(), intersection.begin());
      const auto size = static_cast<double>(std::distance(intersection.begin(), it));

      eff = std::max(eff, size / refTrack.size());
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
