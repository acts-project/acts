// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "NativeSimulationSession.hpp"

#include <stdexcept>

#include <G4Types.hh>

#ifdef G4MULTITHREADED
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Geant4/Geant4Manager.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <mutex>
#include <optional>
#include <utility>

#include <G4Event.hh>
#include <G4FieldManager.hh>
#include <G4LogicalVolume.hh>
#include <G4MTRunManager.hh>
#include <G4Navigator.hh>
#include <G4PropagatorInField.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4RunManagerKernel.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>
#include <G4TransportationManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4VUserActionInitialization.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserEventInformation.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <Randomize.hh>

#include "NativeEventQueue.hpp"
#include "SimulationActions.hpp"

namespace ActsExamples::Geant4 {
namespace {

// This object is owned by Geant4. Drop the binding before publishing the
// result: a waiting ACTS task may immediately destroy the input collections.
struct NativeEvent final : G4VUserEventInformation {
  explicit NativeEvent(std::shared_ptr<NativeEventRequest> job)
      : request(std::move(job)) {}
  void Print() const override {}

  struct Binding {
    Binding(EventContext& context, EventStore& store)
        : guard(context.bind(store)) {}
    decltype(std::declval<EventContext&>().bind(
        std::declval<EventStore&>())) guard;
  };
  std::shared_ptr<NativeEventRequest> request;
  std::optional<Binding> binding;
  std::exception_ptr error;
};

class StreamingRunManager final : public G4MTRunManager {
 public:
  explicit StreamingRunManager(NativeEventQueue& queue) : m_queue(queue) {
    SetEventModulo(1);
  }

  ~StreamingRunManager() override {
    // Geant4's process-wide volume stores otherwise clean up at program exit,
    // on an unrelated thread lacking the master's split geometry data.
    // Join workers before removing the geometry they share.
    TerminateWorkers();
    if (geometryInitialized) {
      // Geant4 11.4 removes named region roots after deleting logical
      // volumes. Removing a root from a multi-root region writes to that
      // volume. Detach the roots while all volumes are still alive.
      for (auto* region : *G4RegionStore::GetInstance()) {
        if (region->GetName() == "DefaultRegionForTheWorld") {
          continue;
        }
        while (region->GetNumberOfRootVolumes() != 0) {
          region->RemoveRootLogicalVolume(
              *region->GetRootLogicalVolumeIterator(), false);
        }
      }
      ReinitializeGeometry(true, false);
    }
  }

  // This virtual is called by each worker, before creating primaries. The
  // sequencer can skip events before transport and does not know the count.
  G4bool SetUpAnEvent(G4Event* event, G4long& s1, G4long& s2, G4long& s3,
                      G4bool /*reseedRequired*/) override {
    if (fakeRun) {
      return false;
    }
    auto request = m_queue.take();
    if (!request) {
      return false;
    }
    // Geant4 initializes the worker engine before GeneratePrimaries. That
    // callback applies the same full ACTS seed used by the serial backend.
    s1 = 1 + static_cast<long>(request->seed % 900000000);
    s2 = 1 + static_cast<long>((request->seed >> 32) % 900000000);
    s3 = 0;
    event->SetEventID(request->id);
    event->SetUserInformation(new NativeEvent(std::move(request)));
    return true;
  }

 protected:
  G4bool InitializeSeeds(G4int /*events*/) override { return true; }

 private:
  NativeEventQueue& m_queue;
};

// Geant4 owns a distinct generator and action list on every worker thread.
class Worker final : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Worker(NativeEventQueue& queue) : m_queue(queue) {}
  std::shared_ptr<EventContext> context = std::make_shared<EventContext>();
  std::unique_ptr<SimulationActions> actions;
  NativeEvent* current = nullptr;

  void GeneratePrimaries(G4Event* event) override {
    current = static_cast<NativeEvent*>(event->GetUserInformation());
    try {
      current->binding.emplace(*context, current->request->store);
      G4Random::setTheSeed(current->request->seed);
      actions->generator->GeneratePrimaries(event);
      current->error = current->request->store.primaryError;
    } catch (...) {
      current->error = std::current_exception();
      // GeneratePrimaries precedes EventProc. BeginOfEventAction performs
      // AbortEvent once Geant4 has established the event.
    }
  }

  // Never unwind an ACTS callback through a Geant4 worker thread.
  template <typename F>
  void duringEvent(F&& callback) {
    if (current != nullptr && current->error) {
      return;
    }
    try {
      callback();
    } catch (...) {
      if (current == nullptr) {
        m_queue.fail(std::current_exception());
      } else {
        current->error = std::current_exception();
        G4RunManager::GetRunManager()->AbortEvent();
      }
    }
  }

 private:
  NativeEventQueue& m_queue;
};

class EventAction final : public G4UserEventAction {
 public:
  explicit EventAction(Worker& worker) : m_worker(worker) {}
  void BeginOfEventAction(const G4Event* /*event*/) override {
    if (m_worker.current->error) {
      G4RunManager::GetRunManager()->AbortEvent();
    }
  }
  void EndOfEventAction(const G4Event* event) override {
    auto* data = std::exchange(m_worker.current, nullptr);
    data->binding.reset();
    auto request = std::move(data->request);
    if (data->error) {
      request->response.set_exception(data->error);
    } else if (event->IsAborted()) {
      request->response.set_exception(
          std::make_exception_ptr(std::runtime_error("Geant4 event aborted")));
    } else {
      request->response.set_value(std::move(request->store.result));
    }
  }

 private:
  Worker& m_worker;
};

class Tracking final : public G4UserTrackingAction {
 public:
  explicit Tracking(Worker& worker) : m_worker(worker) {}
  void SetTrackingManagerPointer(G4TrackingManager* manager) override {
    G4UserTrackingAction::SetTrackingManagerPointer(manager);
    if (m_worker.actions) {
      m_worker.actions->tracking->SetTrackingManagerPointer(manager);
    }
  }
  void PreUserTrackingAction(const G4Track* track) override {
    m_worker.duringEvent(
        [&] { m_worker.actions->tracking->PreUserTrackingAction(track); });
  }
  void PostUserTrackingAction(const G4Track* track) override {
    m_worker.duringEvent(
        [&] { m_worker.actions->tracking->PostUserTrackingAction(track); });
  }

 private:
  Worker& m_worker;
};

class Stepping final : public G4UserSteppingAction {
 public:
  explicit Stepping(Worker& worker) : m_worker(worker) {}
  void SetSteppingManagerPointer(G4SteppingManager* manager) override {
    G4UserSteppingAction::SetSteppingManagerPointer(manager);
    if (m_worker.actions) {
      m_worker.actions->stepping->SetSteppingManagerPointer(manager);
    }
  }
  void UserSteppingAction(const G4Step* step) override {
    m_worker.duringEvent(
        [&] { m_worker.actions->stepping->UserSteppingAction(step); });
  }

 private:
  Worker& m_worker;
};

class WorkerField final : public G4MagneticField {
 public:
  WorkerField(Worker& worker,
              std::shared_ptr<const Acts::MagneticFieldProvider> provider)
      : m_worker(worker),
        m_field(MagneticFieldWrapper::Config{std::move(provider)}) {}
  void GetFieldValue(const G4double point[4], G4double* output) const override {
    std::fill_n(output, 3, 0.);
    m_worker.duringEvent([&] { m_field.GetFieldValue(point, output); });
  }

 private:
  Worker& m_worker;
  MagneticFieldWrapper m_field;
};

class WorkerRun final : public G4UserRunAction {
 public:
  WorkerRun(const Geant4Simulation::Config& cfg, Worker& worker)
      : m_cfg(cfg), m_worker(worker) {}
  void BeginOfRunAction(const G4Run* /*run*/) override {
    m_worker.duringEvent([&] {
      auto* transport = G4TransportationManager::GetTransportationManager();
      if (std::isfinite(m_cfg.propagatorLargestAcceptableStep)) {
        transport->GetPropagatorInField()->SetLargestAcceptableStep(
            m_cfg.propagatorLargestAcceptableStep / Acts::UnitConstants::mm *
            CLHEP::mm);
      }
      if (m_cfg.magneticField) {
        m_field = std::make_unique<WorkerField>(m_worker, m_cfg.magneticField);
        auto* manager = transport->GetFieldManager();
        manager->SetDetectorField(m_field.get());
        manager->CreateChordFinder(m_field.get());
        transport->GetNavigatorForTracking()
            ->GetWorldVolume()
            ->GetLogicalVolume()
            ->SetFieldManager(manager, true);
      }
    });
  }
  void EndOfRunAction(const G4Run* /*run*/) override {
    if (m_field) {
      G4TransportationManager::GetTransportationManager()
          ->GetFieldManager()
          ->SetDetectorField(nullptr);
    }
  }

 private:
  const Geant4Simulation::Config& m_cfg;
  Worker& m_worker;
  std::unique_ptr<WorkerField> m_field;
};

class ActionInitialization final : public G4VUserActionInitialization {
 public:
  ActionInitialization(const Geant4Simulation::Config& cfg,
                       const Acts::Logger& logger,
                       SensitiveSurfaceMapper::VolumeToSurfAssocMap_t mapping,
                       NativeEventQueue& queue)
      : m_cfg(cfg),
        m_logger(logger),
        m_mapping(std::make_shared<
                  const SensitiveSurfaceMapper::VolumeToSurfAssocMap_t>(
            std::move(mapping))),
        m_queue(queue) {}
  void Build() const override {
    auto* worker = new Worker(m_queue);
    SetUserAction(worker);
    try {
      worker->actions =
          std::make_unique<SimulationActions>(m_cfg, worker->context, m_logger);
      worker->actions->sensitive->shareSurfaceMapping(m_mapping);
      // Invoke potentially stateful factories sequentially during worker setup.
      std::lock_guard lock(m_factoryMutex);
      appendDiagnostics(*worker->actions->stepping, m_cfg, worker->context);
    } catch (...) {
      m_queue.fail(std::current_exception());
      return;
    }
    SetUserAction(new EventAction(*worker));
    SetUserAction(new Tracking(*worker));
    SetUserAction(new Stepping(*worker));
    SetUserAction(new WorkerRun(m_cfg, *worker));
  }

 private:
  const Geant4Simulation::Config& m_cfg;
  const Acts::Logger& m_logger;
  const std::shared_ptr<const SensitiveSurfaceMapper::VolumeToSurfAssocMap_t>
      m_mapping;
  NativeEventQueue& m_queue;
  mutable std::mutex m_factoryMutex;
};
}  // namespace

struct NativeSimulationSession::Impl {
  Impl(const Geant4Simulation::Config& config, const Acts::Logger& log)
      : cfg(config), logger(log.clone()), queue(cfg.numGeant4Workers) {}
  const Geant4Simulation::Config cfg;
  std::unique_ptr<const Acts::Logger> logger;
  NativeEventQueue queue;
  std::future<void> master;

  void run() {
    ActsPlugins::FpeMonitor monitor{0};
    // Keep the master alive until waiting workers and callers have been woken
    // on failure, before its destructor joins the Geant4 threads.
    std::optional<StreamingRunManager> runManager;
    try {
      auto& factory = Geant4Manager::instance();
      factory.reserveRunManager();
      auto& manager = runManager.emplace(queue);
      const int workers = static_cast<int>(cfg.numGeant4Workers);
      manager.SetNumberOfThreads(workers);
      if (manager.GetNumberOfThreads() != workers) {
        throw std::runtime_error(
            "G4FORCENUMBEROFTHREADS conflicts with numGeant4Workers");
      }
      manager.SetUserInitialization(
          factory.createPhysicsList(cfg.physicsList).release());
      manager.SetUserInitialization(
          cfg.detector->buildGeant4DetectorConstruction(cfg.constructionOptions)
              .release());
      manager.InitializeGeometry();
      auto* world =
          G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld();
      manager.SetUserInitialization(new ActionInitialization(
          cfg, *logger, mapSensitiveVolumes(cfg, *world, *logger), queue));
      Geant4Manager::tweakLogging(
          manager, logger->level() == Acts::Logging::VERBOSE ? 2 : 0);
      manager.Initialize();
      queue.checkError();
      // BeamOn starts Geant4's workers and retains all run barriers. The
      // master's own loop is bounded by the pool size; workers obtain real
      // events via SetUpAnEvent until endRun closes the queue.
      manager.BeamOn(workers);
      queue.checkFinished();
    } catch (...) {
      auto error = std::current_exception();
      queue.fail(error);
      throw;
    }
  }
};

NativeSimulationSession::NativeSimulationSession(
    const Geant4Simulation::Config& cfg, const Acts::Logger& logger)
    : m_impl(std::make_unique<Impl>(cfg, logger)) {
  if (cfg.numGeant4Workers == 0 ||
      cfg.numGeant4Workers >
          static_cast<unsigned int>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("Invalid native Geant4 worker count");
  }
}

NativeSimulationSession::~NativeSimulationSession() {
  m_impl->queue.close();
  if (m_impl->master.valid()) {
    m_impl->master.wait();
  }
}

void NativeSimulationSession::beginRun() {
  if (m_impl->master.valid()) {
    throw std::logic_error("Native Geant4 run already initialized");
  }
  m_impl->master = std::async(std::launch::async, &Impl::run, m_impl.get());
  try {
    m_impl->queue.waitUntilOpen();
  } catch (...) {
    // Initialization failed. Wait for teardown while config, geometry and
    // logging are still alive, and consume the master's completion state.
    m_impl->master.get();
    throw;
  }
}

void NativeSimulationSession::endRun() {
  if (G4Threading::IsWorkerThread()) {
    throw std::logic_error(
        "Cannot finalize the Geant4 run from a transport callback");
  }
  m_impl->queue.close();
  if (m_impl->master.valid()) {
    m_impl->master.get();
  }
}

EventResult NativeSimulationSession::simulate(
    const SimParticleContainer& particles,
    const Acts::GeometryContext& geometry, std::uint64_t seed,
    std::size_t eventNumber, const SimParticleDecayGraph* decayGraph) {
  if (G4Threading::IsWorkerThread()) {
    throw std::logic_error(
        "Cannot submit a Geant4 event from a transport callback");
  }
  if (eventNumber > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error(
        "ACTS event number exceeds the Geant4 event ID range");
  }
  return m_impl->queue.submit(particles, geometry, decayGraph, seed,
                              static_cast<int>(eventNumber));
}

}  // namespace ActsExamples::Geant4

#else
namespace ActsExamples::Geant4 {
struct NativeSimulationSession::Impl {};
NativeSimulationSession::NativeSimulationSession(
    const Geant4Simulation::Config& /*cfg*/, const Acts::Logger& /*logger*/) {
  throw std::runtime_error(
      "Native Geant4 workers require a multithreaded Geant4 build");
}
NativeSimulationSession::~NativeSimulationSession() = default;
void NativeSimulationSession::beginRun() {}
void NativeSimulationSession::endRun() {}
EventResult NativeSimulationSession::simulate(
    const SimParticleContainer& /*particles*/,
    const Acts::GeometryContext& /*geometry*/, std::uint64_t /*seed*/,
    std::size_t /*eventNumber*/, const SimParticleDecayGraph* /*decayGraph*/) {
  throw std::logic_error("Native Geant4 workers unavailable");
}
}  // namespace ActsExamples::Geant4
#endif
