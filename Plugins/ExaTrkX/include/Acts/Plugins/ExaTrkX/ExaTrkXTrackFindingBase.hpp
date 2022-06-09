#pragma once

#include <string>
#include <vector>

#include "ExaTrkXTiming.hpp"

namespace Acts {

class ExaTrkXTrackFindingBase {
 public:
  ExaTrkXTrackFindingBase(std::string name, bool verbose)
      : m_name(std::move(name)), m_verbose(verbose) {}
  virtual ~ExaTrkXTrackFindingBase() {}

  ExaTrkXTrackFindingBase() = delete;
  ExaTrkXTrackFindingBase(const ExaTrkXTrackFindingBase&) = delete;
  ExaTrkXTrackFindingBase& operator=(const ExaTrkXTrackFindingBase&) = delete;

  void getTracks(std::vector<float>& inputValues,
                 std::vector<int>& spacepointIDs,
                 std::vector<std::vector<int> >& trackCandidates) const;

  /// implement the algorithm with timing information collected.
  virtual void getTracks(std::vector<float>& inputValues,
                         std::vector<int>& spacepointIDs,
                         std::vector<std::vector<int> >& trackCandidates,
                         ExaTrkXTime& timeInfo) const = 0;

  std::string name() const { return m_name; }

 protected:
  std::string m_name;
  bool m_verbose = false;
};

}  // namespace Acts
