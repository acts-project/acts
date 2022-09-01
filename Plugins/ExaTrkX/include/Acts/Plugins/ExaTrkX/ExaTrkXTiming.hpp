#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <chrono>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

namespace Acts {

struct ExaTrkXTime {
  float embedding = 0.0;
  float building = 0.0;
  float filtering = 0.0;
  float gnn = 0.0;
  float labeling = 0.0;
  float total = 0.0;

  void summary(LoggerWrapper& logger) const {
    ACTS_VERBOSE("1) embedding: " << embedding);
    ACTS_VERBOSE("2) building: " << building);
    ACTS_VERBOSE("3) filtering: " << filtering);
    ACTS_VERBOSE("4) gnn: " << gnn);
    ACTS_VERBOSE("5) labeling: " << labeling);
    ACTS_VERBOSE("6) total: " << total);
  }
};

struct ExaTrkXTimeList {
  std::vector<float> embedding;
  std::vector<float> building;
  std::vector<float> filtering;
  std::vector<float> gnn;
  std::vector<float> labeling;
  std::vector<float> total;

  void add(const ExaTrkXTime& time) {
    embedding.push_back(time.embedding);
    building.push_back(time.building);
    filtering.push_back(time.filtering);
    gnn.push_back(time.gnn);
    labeling.push_back(time.labeling);
    total.push_back(time.total);
  }

  ExaTrkXTime get(size_t evtid) {
    if (evtid >= embedding.size()) {
      throw std::runtime_error("Error: event id is out of range.");
    }
    ExaTrkXTime timing{embedding[evtid], building[evtid], filtering[evtid],
                       gnn[evtid],       labeling[evtid], total[evtid]};
    return timing;
  }

  void summary(Acts::LoggerWrapper logger, size_t start = 0) {
    size_t num = embedding.size();
    if (num <= start) {
      throw std::runtime_error("Not enough data events");
    }
    num -= start;
    float tot_embedding = 0;
    float tot_building = 0;
    float tot_filtering = 0;
    float tot_gnn = 0;
    float tot_labeling = 0;
    float tot_total = 0;

    if (num > 0) {
      tot_embedding =
          std::accumulate(embedding.begin() + start, embedding.end(), 0.0f);
      tot_building =
          std::accumulate(building.begin() + start, building.end(), 0.0f);
      tot_filtering =
          std::accumulate(filtering.begin() + start, filtering.end(), 0.0f);
      tot_gnn = std::accumulate(gnn.begin() + start, gnn.end(), 0.0f);
      tot_labeling =
          std::accumulate(labeling.begin() + start, labeling.end(), 0.0f);
      tot_total = std::accumulate(total.begin() + start, total.end(), 0.0f);
    }

    ACTS_VERBOSE("1) embedding: " << tot_embedding / num);
    ACTS_VERBOSE("2) building: " << tot_building / num);
    ACTS_VERBOSE("3) filtering: " << tot_filtering / num);
    ACTS_VERBOSE("4) gnn: " << tot_gnn / num);
    ACTS_VERBOSE("5) labeling: " << tot_labeling / num);
    ACTS_VERBOSE("6) total: " << tot_total / num);
  }

  void summaryOneEvent(Acts::LoggerWrapper logger, int evtid) {
    get(evtid).summary(logger);
  }

  std::size_t numEvts() { return embedding.size(); }

  void save(const std::string& filename) {
    std::ofstream file(filename);
    if (not file.is_open()) {
      throw std::runtime_error("Cannot write to file " + filename);
    }

    file << "embedding,building,filtering,gnn,labeling,total\n";

    for (auto i = 0ul; i < embedding.size(); i++) {
      file << embedding[i] << building[i] << filtering[i] << gnn[i]
           << labeling[i] << total[i];
    }
  }
};

class ExaTrkXTimer {
 public:
  void start() {
    m_start = std::chrono::high_resolution_clock::now();
    m_running = true;
  }
  void stop() {
    m_end = std::chrono::high_resolution_clock::now();
    m_running = false;
  }
  double stopAndGetElapsedTime() {
    stop();
    return elapsedSeconds();
  }
  double elapsed() {
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    if (m_running) {
      end = std::chrono::high_resolution_clock::now();
    } else {
      end = m_end;
    }
    return std::chrono::duration<double, std::milli>(end - m_start).count();
  }
  double elapsedSeconds() { return elapsed() / 1000.0; }

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
  bool m_running = false;
};

}  // namespace Acts
