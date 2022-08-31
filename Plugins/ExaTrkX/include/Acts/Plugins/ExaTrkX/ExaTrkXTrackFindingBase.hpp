#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <string>
#include <vector>

namespace Acts {

struct ExaTrkXTime;

class ExaTrkXTrackFindingBase {
 public:
  ExaTrkXTrackFindingBase(const std::string &name)
      : m_name(name) {}
      
  virtual ~ExaTrkXTrackFindingBase() {}

  ExaTrkXTrackFindingBase() = delete;
  ExaTrkXTrackFindingBase(const ExaTrkXTrackFindingBase&) = delete;
  ExaTrkXTrackFindingBase& operator=(const ExaTrkXTrackFindingBase&) = delete;

  void getTracks(std::vector<float>& inputValues,
                 std::vector<int>& spacepointIDs,
                 std::vector<std::vector<int> >& trackCandidates,
                 Acts::LoggerWrapper logger = Acts::getDummyLogger()) const;

  /// implement the algorithm with timing information collected.
  virtual void getTracks(
      std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
      std::vector<std::vector<int> >& trackCandidates, ExaTrkXTime& timeInfo,
      Acts::LoggerWrapper logger = Acts::getDummyLogger()) const = 0;

  const std::string &name() const { return m_name; }

 protected:
  std::string m_name;
};

}  // namespace Acts
