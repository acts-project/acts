// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/Pythia8ProcessGenerator.hpp"

#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <algorithm>
#include <iterator>
#include <ostream>
#include <random>
#include <utility>

#include <HepMC3/WriterAscii.h>
#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC3.h>

namespace ActsExamples {

struct Pythia8RandomEngineWrapper : public Pythia8::RndmEngine {
  RandomEngine* rng{nullptr};

  struct {
    std::size_t numUniformRandomNumbers = 0;
    double first = std::numeric_limits<double>::quiet_NaN();
    double last = std::numeric_limits<double>::quiet_NaN();
  } statistics;

  Pythia8RandomEngineWrapper() = default;

  double flat() override {
    if (rng == nullptr) {
      throw std::runtime_error(
          "Pythia8RandomEngineWrapper: no random engine set");
    }

    double value = std::uniform_real_distribution<double>(0.0, 1.0)(*rng);
    if (statistics.numUniformRandomNumbers == 0) {
      statistics.first = value;
    }
    statistics.last = value;
    statistics.numUniformRandomNumbers++;
    return value;
  }

  void setRandomEngine(RandomEngine& rng_) { rng = &rng_; }
  void clearRandomEngine() { rng = nullptr; }
};

struct Pythia8GeneratorImpl {
  std::unique_ptr<HepMC3::Writer> m_hepMC3Writer;
  std::unique_ptr<HepMC3::Pythia8ToHepMC3> m_hepMC3Converter;
  std::shared_ptr<Pythia8RandomEngineWrapper> m_pythia8RndmEngine;
};

Pythia8Generator::Pythia8Generator(const Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("Pythia8Generator", lvl)),
      m_pythia8(std::make_unique<Pythia8::Pythia>("", false)) {
  ACTS_INFO("Pythia8Generator: init");
  // disable all output by default but allow re-enable via config
  m_pythia8->settings.flag("Print:quiet", true);
  for (const auto& setting : m_cfg.settings) {
    ACTS_VERBOSE("use Pythia8 setting '" << setting << "'");
    m_pythia8->readString(setting.c_str());
  }
  m_pythia8->settings.mode("Beams:idA", m_cfg.pdgBeam0);
  m_pythia8->settings.mode("Beams:idB", m_cfg.pdgBeam1);
  m_pythia8->settings.mode("Beams:frameType", 1);
  m_pythia8->settings.parm("Beams:eCM",
                           m_cfg.cmsEnergy / Acts::UnitConstants::GeV);

  m_impl = std::make_unique<Pythia8GeneratorImpl>();

  m_impl->m_pythia8RndmEngine = std::make_shared<Pythia8RandomEngineWrapper>();

#if PYTHIA_VERSION_INTEGER >= 8310
  m_pythia8->setRndmEnginePtr(m_impl->m_pythia8RndmEngine);
#else
  m_pythia8->setRndmEnginePtr(m_impl->m_pythia8RndmEngine.get());
#endif

  RandomEngine rng{m_cfg.initializationSeed};
  m_impl->m_pythia8RndmEngine->setRandomEngine(rng);
  m_pythia8->init();
  m_impl->m_pythia8RndmEngine->clearRandomEngine();

  m_impl->m_hepMC3Converter = std::make_unique<HepMC3::Pythia8ToHepMC3>();
  if (m_cfg.writeHepMC3.has_value()) {
    ACTS_DEBUG("Initializing HepMC3 output to: " << m_cfg.writeHepMC3.value());
    m_impl->m_hepMC3Writer =
        std::make_unique<HepMC3::WriterAscii>(m_cfg.writeHepMC3.value());
  }
}

// needed to allow unique_ptr of forward-declared Pythia class
Pythia8Generator::~Pythia8Generator() {
  if (m_impl->m_hepMC3Writer) {
    m_impl->m_hepMC3Writer->close();
  }

  ACTS_DEBUG("Pythia8Generator produced "
             << m_impl->m_pythia8RndmEngine->statistics.numUniformRandomNumbers
             << " uniform random numbers");
  ACTS_DEBUG("                 first = "
             << m_impl->m_pythia8RndmEngine->statistics.first);
  ACTS_DEBUG("                  last = "
             << m_impl->m_pythia8RndmEngine->statistics.last);
}

std::shared_ptr<HepMC3::GenEvent> Pythia8Generator::operator()(
    RandomEngine& rng) {
  using namespace Acts::UnitLiterals;

  // pythia8 is not thread safe and generation needs to be protected
  std::lock_guard<std::mutex> lock(m_pythia8Mutex);
  // use per-thread random engine also in pythia

  m_impl->m_pythia8RndmEngine->setRandomEngine(rng);

  {
    Acts::FpeMonitor mon{0};  // disable all FPEs while we're in Pythia8
    m_pythia8->next();
  }

  if (m_cfg.printShortEventListing) {
    m_pythia8->process.list();
  }
  if (m_cfg.printLongEventListing) {
    m_pythia8->event.list();
  }

  auto genEvent = std::make_shared<HepMC3::GenEvent>();
  genEvent->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

  assert(m_impl->m_hepMC3Converter != nullptr);
  m_impl->m_hepMC3Converter->fill_next_event(*m_pythia8, genEvent.get(),
                                             genEvent->event_number());

  if (m_impl->m_hepMC3Converter && m_impl->m_hepMC3Writer) {
    m_impl->m_hepMC3Writer->write_event(*genEvent);
  }

  m_impl->m_pythia8RndmEngine->clearRandomEngine();

  return genEvent;
}

}  // namespace ActsExamples
