#include "ActsExamples/TrackFindingX/ProtoTrackEfficiencyPrinter.hpp"

#include <Acts/Utilities/Zip.hpp>

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

  std::sort(effs.begin(), effs.end());

  const static std::vector<double> thresholds = {0., .1, .2, .3, .4,
                                                 .5, .6, .7, .8, .9};

  auto it = effs.begin();
  std::vector<std::size_t> hist;
  for (double threshold : thresholds) {
    auto endIt = std::find_if(it, effs.end(),
                              [&](double eff) { return eff > threshold; });
    hist.push_back(std::distance(it, endIt));
    it = endIt;
  }

  const auto max = *std::max_element(hist.begin(), hist.end());
  const int colMax = 40;

  ACTS_INFO("Prototrack efficiency histogram:");
  for (const auto &[v, th] : Acts::zip(hist, thresholds)) {
    auto rel = v / static_cast<float>(max);
    auto l = std::round(rel * colMax);
    std::string str(l, '#');

    auto label = std::to_string(th);
    label.resize(7, ' ');

    ACTS_INFO(">" << label << " | " << str);
  }

  return {};
}
