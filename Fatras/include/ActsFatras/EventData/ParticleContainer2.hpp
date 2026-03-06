// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/EnumBitwiseOperators.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "Acts/Utilities/detail/ContainerSubset.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/GenerationProcess.hpp"
#include "ActsFatras/EventData/SimulationOutcome.hpp"

#include <cstdint>
#include <span>

namespace Acts {
class Surface;
}

namespace ActsFatras {

using ParticleIndex2 = std::uint32_t;
using ParticleIndexSubset2 = std::span<const ParticleIndex2>;

static constexpr ParticleIndex2 InvalidParticleIndex2 =
    std::numeric_limits<ParticleIndex2>::max();

using HitIndex2 = std::uint32_t;

class ParticleContainer2;

template <bool>
class ParticleProxy2;
using MutableParticleProxy2 = ParticleProxy2<false>;
using ConstParticleProxy2 = ParticleProxy2<true>;

namespace detail {
struct GenerationStateAccessor;
struct SimulationStateAccessor;
}  // namespace detail

template <typename, bool>
class ParticleStateProxy2;
using MutableGeneratorParticleProxy2 =
    ParticleStateProxy2<detail::GenerationStateAccessor, false>;
using ConstGeneratorParticleProxy2 =
    ParticleStateProxy2<detail::GenerationStateAccessor, true>;
using MutableSimulationParticleProxy2 =
    ParticleStateProxy2<detail::SimulationStateAccessor, false>;
using ConstSimulationParticleProxy2 =
    ParticleStateProxy2<detail::SimulationStateAccessor, true>;

template <typename T, bool>
class ParticleColumnProxy2;
template <typename T>
using MutableParticleColumnProxy2 = ParticleColumnProxy2<T, false>;
template <typename T>
using ConstParticleColumnProxy2 = ParticleColumnProxy2<T, true>;

enum class ParticleColumns2 : std::uint32_t {
  None = 0,

  // General particle properties
  Parents = 1 << 0,
  Barcode = 1 << 1,
  Pdg = 1 << 2,
  Charge = 1 << 3,
  Mass = 1 << 4,
  GenerationProcess = 1 << 5,

  // Initial kinematic properties at production vertex
  InitialReferenceSurface = 1 << 6,
  InitialFourPosition = 1 << 7,
  InitialAbsoluteMomentum = 1 << 8,
  InitialDirection = 1 << 9,

  // Simulation-specific properties
  CurrentReferenceSurface = 1 << 10,
  CurrentFourPosition = 1 << 11,
  CurrentAbsoluteMomentum = 1 << 12,
  CurrentDirection = 1 << 13,
  ProperTime = 1 << 14,
  PathInX0 = 1 << 15,
  PathInL0 = 1 << 16,
  NumberOfHits = 1 << 17,
  SimulationOutcome = 1 << 18,

  Generated = Parents | Barcode | Pdg | Charge | Mass | GenerationProcess |
              InitialReferenceSurface | InitialFourPosition |
              InitialAbsoluteMomentum | InitialDirection,
  Simulated = Generated | CurrentReferenceSurface | CurrentFourPosition |
              CurrentAbsoluteMomentum | CurrentDirection | ProperTime |
              PathInX0 | PathInL0 | NumberOfHits | SimulationOutcome,
  All = Parents | Barcode | Pdg | Charge | Mass | GenerationProcess |
        InitialReferenceSurface | InitialFourPosition |
        InitialAbsoluteMomentum | InitialDirection | CurrentReferenceSurface |
        CurrentFourPosition | CurrentAbsoluteMomentum | CurrentDirection |
        ProperTime | PathInX0 | PathInL0 | NumberOfHits | SimulationOutcome,
};

/// Enable bitwise operators for ParticleColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(ParticleColumns2);

namespace detail {

class ColumnHolderBase {
 public:
  virtual ~ColumnHolderBase() = default;

  virtual std::unique_ptr<ColumnHolderBase> copy() const = 0;

  virtual std::size_t size() const = 0;
  virtual void reserve(std::size_t size) = 0;
  virtual void resize(std::size_t size) = 0;
  virtual void clear() = 0;
  virtual void emplace_back() = 0;
};

template <typename T>
class ColumnHolder final : public ColumnHolderBase {
 public:
  using Value = T;
  using Container = std::vector<Value>;
  using MutableProxy = MutableParticleColumnProxy2<Value>;
  using ConstProxy = ConstParticleColumnProxy2<Value>;

  ColumnHolder() = default;
  explicit ColumnHolder(Value defaultValue)
      : m_default(std::move(defaultValue)) {}

  MutableProxy proxy(ParticleContainer2 &container) {
    return MutableProxy(container, m_data);
  }
  ConstProxy proxy(const ParticleContainer2 &container) const {
    return ConstProxy(container, m_data);
  }

  std::unique_ptr<ColumnHolderBase> copy() const override {
    return std::make_unique<ColumnHolder<T>>(*this);
  }

  std::size_t size() const override { return m_data.size(); }
  void reserve(std::size_t size) override { m_data.reserve(size); }
  void clear() override { m_data.clear(); }
  void resize(std::size_t size) override { m_data.resize(size, m_default); }
  void emplace_back() override { m_data.emplace_back(m_default); }

 private:
  Value m_default{};
  Container m_data;
};

}  // namespace detail

class ParticleContainer2 final {
 public:
  /// Type alias for particle index in container
  using Index = ParticleIndex2;
  /// Type alias for subset of particle indices
  using IndexSubset = ParticleIndexSubset2;
  /// Type alias for mutable particle proxy
  using MutableProxy = MutableParticleProxy2;
  /// Type alias for const particle proxy
  using ConstProxy = ConstParticleProxy2;

  explicit ParticleContainer2(
      ParticleColumns2 columns = ParticleColumns2::None) noexcept;

  /// Constructs a copy of the given space point container.
  /// @param other The space point container to copy.
  ParticleContainer2(const ParticleContainer2 &other) noexcept;

  /// Move constructs a space point container.
  /// @param other The space point container to move.
  ParticleContainer2(ParticleContainer2 &&other) noexcept;

  /// Detructs the space point container.
  ~ParticleContainer2() noexcept = default;

  /// Assignment operator for copying a space point container.
  /// @param other The space point container to copy.
  /// @return A reference to this space point container.
  ParticleContainer2 &operator=(const ParticleContainer2 &other) noexcept;

  /// Move assignment operator for a space point container.
  /// @param other The space point container to move.
  /// @return A reference to this space point container.
  ParticleContainer2 &operator=(ParticleContainer2 &&other) noexcept;

  /// Returns the number of particles in the container.
  /// @return The number of particles in the container.
  [[nodiscard]] std::uint32_t size() const noexcept { return m_size; }
  /// Checks if the container is empty.
  /// @return True if the container is empty, false otherwise.
  [[nodiscard]] bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of particles.
  /// @param size The number of particles to reserve space for.
  /// @param averageParentIndices The average number of parent indices per particle.
  void reserve(std::uint32_t size, float averageParentIndices = 1) noexcept;

  /// Clears the container, removing all particles and columns.
  void clear() noexcept;

  /// Creates a new particle at the end of the container.
  /// @return A mutable proxy to the newly created particle.
  MutableProxy createParticle() noexcept;

  /// Creates additional columns. This will create the columns if they do not
  /// already exist.
  /// @param columns The columns to create.
  void createColumns(ParticleColumns2 columns) noexcept;

  /// Drops the specified columns from the container.
  /// This will only drop columns if they exist.
  /// @param columns The columns to drop.
  void dropColumns(ParticleColumns2 columns) noexcept;

  /// Checks if the container has the given Columns.
  /// @param columns The Columns to check for.
  /// @return True if the container has all the specified Columns, false
  ///         otherwise.
  bool hasColumns(ParticleColumns2 columns) const noexcept {
    return (m_knownColumns & columns) == columns;
  }

  /// Creates a new column with the given name.
  /// If a column with the same name already exists, an exception is thrown.
  /// @param name The name of the column.
  /// @return A reference to the newly created column.
  /// @throws std::runtime_error if a column with the same name already exists.
  /// @throws std::runtime_error if the column name is reserved.
  template <typename T>
  MutableParticleColumnProxy2<T> createColumn(const std::string &name) {
    return createColumnImpl<ColumnHolder<T>>(name);
  }

  /// Drops the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @throws std::runtime_error if the column does not exist.
  /// @throws std::runtime_error if the column name is reserved.
  void dropColumn(const std::string &name);

  /// Checks if an Column with the given name exists.
  /// @param name The name of the column.
  /// @return True if the column exists, false otherwise.
  bool hasColumn(const std::string &name) const noexcept {
    return m_allColumns.contains(name);
  }

  /// Returns a mutable reference to the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A mutable reference to the column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  MutableParticleColumnProxy2<T> column(const std::string &name) {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a const reference to the column with the given name.
  /// If the column does not exist, an exception is thrown.
  /// @param name The name of the column.
  /// @return A const reference to the column.
  /// @throws std::runtime_error if the column does not exist.
  template <typename T>
  ConstParticleColumnProxy2<T> column(const std::string &name) const {
    return columnImpl<ColumnHolder<T>>(name);
  }

  /// Returns a mutable proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the space point at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A mutable proxy to the space point at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the space point at the given index.
  /// @param index The index of the space point to access.
  /// @return A const proxy to the space point at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Type alias for template iterator over space points in container
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      ParticleContainer2,
      std::conditional_t<read_only, ConstProxy, MutableProxy>, Index,
      read_only>;

  /// Type alias for mutable iterator over space points
  using iterator = Iterator<false>;
  /// Type alias for const iterator over space points
  using const_iterator = Iterator<true>;

  /// @brief Returns mutable iterator to the beginning of the container
  /// @return Mutable iterator pointing to the first space point
  iterator begin() noexcept { return iterator(*this, 0); }
  /// @brief Returns mutable iterator to the end of the container
  /// @return Mutable iterator pointing past the last space point
  iterator end() noexcept { return iterator(*this, size()); }

  /// @brief Returns const iterator to the beginning of the container
  /// @return Const iterator pointing to the first space point
  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  /// @brief Returns const iterator to the end of the container
  /// @return Const iterator pointing past the last space point
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

  /// Subset facade over arbitrary index sets.
  template <bool read_only>
  class Subset final
      : public Acts::detail::ContainerSubset<
            Subset<read_only>, Subset<true>, ParticleContainer2,
            std::conditional_t<read_only, ConstProxy, MutableProxy>,
            std::span<const Index>, read_only> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<
        Subset<read_only>, Subset<true>, ParticleContainer2,
        std::conditional_t<read_only, ConstProxy, MutableProxy>,
        std::span<const Index>, read_only>;

    using Base::Base;
  };

  /// Type alias for mutable subset of space points
  using MutableSubset = Subset<false>;
  /// Type alias for const subset of space points
  using ConstSubset = Subset<true>;

  /// Creates a mutable subset of space points from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A mutable subset of space points.
  MutableSubset subset(const IndexSubset &subset) noexcept {
    return MutableSubset(*this, subset);
  }
  /// Creates a const subset of space points from the given index subset.
  /// @param subset The index subset to create the subset from.
  /// @return A const subset of space points.
  ConstSubset subset(const IndexSubset &subset) const noexcept {
    return ConstSubset(*this, subset);
  }

 public:
  using ColumnHolderBase = detail::ColumnHolderBase;
  template <typename T>
  using ColumnHolder = detail::ColumnHolder<T>;

  template <bool>
  friend class ParticleProxy2;
  template <typename, bool>
  friend class ParticleStateProxy2;

  std::size_t m_size{0};

  std::unordered_map<std::string, ColumnHolderBase *> m_allColumns;
  ParticleColumns2 m_knownColumns{ParticleColumns2::None};
  std::unordered_map<std::string, std::unique_ptr<ColumnHolderBase>>
      m_dynamicColumns;

  std::vector<ParticleIndex2> m_parentIndices;

  std::optional<ColumnHolder<std::uint32_t>> m_parentIndicesOffsetColumn;
  std::optional<ColumnHolder<std::uint8_t>> m_parentIndicesCountColumn;
  std::optional<ColumnHolder<Barcode>> m_barcodeColumn;
  std::optional<ColumnHolder<Acts::PdgParticle>> m_pdgColumn;
  std::optional<ColumnHolder<double>> m_chargeColumn;
  std::optional<ColumnHolder<double>> m_massColumn;
  std::optional<ColumnHolder<GenerationProcess>> m_generationProcessColumn;

  std::optional<ColumnHolder<const Acts::Surface *>>
      m_initialReferenceSurfaceColumn;
  std::optional<ColumnHolder<Acts::Vector4>> m_initialFourPositionColumn;
  std::optional<ColumnHolder<double>> m_initialAbsoluteMomentumColumn;
  std::optional<ColumnHolder<Acts::Vector3>> m_initialDirectionColumn;

  std::optional<ColumnHolder<const Acts::Surface *>>
      m_currentReferenceSurfaceColumn;
  std::optional<ColumnHolder<Acts::Vector4>> m_currentFourPositionColumn;
  std::optional<ColumnHolder<double>> m_currentAbsoluteMomentumColumn;
  std::optional<ColumnHolder<Acts::Vector3>> m_currentDirectionColumn;
  std::optional<ColumnHolder<double>> m_properTimeColumn;
  std::optional<ColumnHolder<double>> m_pathInX0Column;
  std::optional<ColumnHolder<double>> m_pathInL0Column;
  std::optional<ColumnHolder<std::uint32_t>> m_numberOfHitsColumn;
  std::optional<ColumnHolder<SimulationOutcome>> m_simulationOutcomeColumn;

  static auto knownColumnMasks() noexcept {
    using enum ParticleColumns2;
    return std::tuple(
        Parents, Parents, Barcode, Pdg, Charge, Mass, GenerationProcess,
        InitialReferenceSurface, InitialFourPosition, InitialAbsoluteMomentum,
        InitialDirection, CurrentReferenceSurface, CurrentFourPosition,
        CurrentAbsoluteMomentum, CurrentDirection, ProperTime, PathInX0,
        PathInL0, NumberOfHits, SimulationOutcome);
  }

  static auto knownColumnNames() noexcept {
    return std::tuple(
        "parentsOffset", "parentsCount", "barcode", "pdg", "charge", "mass",
        "generationProcess", "initialReferenceSurface", "initialFourPosition",
        "initialAbsoluteMomentum", "initialDirection",
        "currentReferenceSurface", "currentFourPosition",
        "currentAbsoluteMomentum", "currentDirection", "properTime", "pathInX0",
        "pathInL0", "numberOfHits", "simulationOutcome");
  }

  static auto knownColumnDefaults() noexcept {
    return std::tuple(std::uint32_t{0}, std::uint8_t{0}, Barcode{},
                      Acts::PdgParticle::eInvalid, double{0}, double{0},
                      GenerationProcess::eUndefined,
                      static_cast<const Acts::Surface *>(nullptr),
                      Acts::Vector4(Acts::Vector4::Zero()), double{0},
                      Acts::Vector3(Acts::Vector3::Zero()),
                      static_cast<const Acts::Surface *>(nullptr),
                      Acts::Vector4(Acts::Vector4::Zero()), double{0},
                      Acts::Vector3(Acts::Vector3::Zero()), double{0},
                      double{0}, double{0}, std::uint32_t{0},
                      SimulationOutcome::Alive);
  }

  template <typename Self>
  static auto knownColumns(Self &&self) noexcept {
    return std::tie(
        self.m_parentIndicesOffsetColumn, self.m_parentIndicesCountColumn,
        self.m_barcodeColumn, self.m_pdgColumn, self.m_chargeColumn,
        self.m_massColumn, self.m_generationProcessColumn,
        self.m_initialReferenceSurfaceColumn, self.m_initialFourPositionColumn,
        self.m_initialAbsoluteMomentumColumn, self.m_initialDirectionColumn,
        self.m_currentReferenceSurfaceColumn, self.m_currentFourPositionColumn,
        self.m_currentAbsoluteMomentumColumn, self.m_currentDirectionColumn,
        self.m_properTimeColumn, self.m_pathInX0Column, self.m_pathInL0Column,
        self.m_numberOfHitsColumn, self.m_simulationOutcomeColumn);
  }
  auto knownColumns() & noexcept { return knownColumns(*this); }
  auto knownColumns() const & noexcept { return knownColumns(*this); }
  auto knownColumns() && noexcept { return knownColumns(*this); }

  void copyColumns(const ParticleContainer2 &other);
  void moveColumns(ParticleContainer2 &other) noexcept;

  static bool reservedColumn(const std::string &name) noexcept;

  template <typename Holder>
  MutableParticleColumnProxy2<typename Holder::Value> createColumnImpl(
      const std::string &name) {
    if (reservedColumn(name)) {
      throw std::runtime_error("Column name is reserved: " + name);
    }
    if (hasColumn(name)) {
      throw std::runtime_error("Column already exists: " + name);
    }
    auto holder = std::make_unique<Holder>();
    holder->resize(size());
    auto proxy = holder->proxy(*this);
    m_allColumns.try_emplace(name, holder.get());
    m_dynamicColumns.try_emplace(name, std::move(holder));
    return proxy;
  }

  template <typename Holder, typename Self>
  static auto columnImpl(Self &&self, const std::string &name) {
    const auto it = self.m_allColumns.find(name);
    if (it == self.m_allColumns.end()) {
      throw std::runtime_error("Column not found: " + name);
    }
    auto &holder = dynamic_cast<Holder &>(*it->second);
    return holder.proxy(self);
  }

  template <typename Holder>
  MutableParticleColumnProxy2<typename Holder::Value> columnImpl(
      const std::string &name) {
    return columnImpl<Holder>(*this, name);
  }

  template <typename Holder>
  ConstParticleColumnProxy2<typename Holder::Value> columnImpl(
      const std::string &name) const {
    return columnImpl<Holder>(*this, name);
  }
};

using ConstParticleSubset = ParticleContainer2::ConstSubset;
using MutableParticleSubset = ParticleContainer2::MutableSubset;

/// A proxy class for accessing individual particles.
template <bool read_only>
class ParticleProxy2 final {
 public:
  /// Indicates whether this particle proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Type alias for particle index type
  using Index = ParticleIndex2;

  /// Type alias for container type (const if read-only)
  using Container = Acts::const_if_t<ReadOnly, ParticleContainer2>;

  /// Constructs a particle proxy for the given container and index.
  /// @param container The container holding the particle.
  /// @param index The index of the particle in the container.
  ParticleProxy2(Container &container, Index index) noexcept
      : m_container(&container), m_index(index) {}

  /// Copy construct a particle proxy.
  /// @param other The particle proxy to copy.
  ParticleProxy2(const ParticleProxy2 &other) noexcept = default;

  /// Copy construct a mutable particle proxy.
  /// @param other The mutable particle proxy to copy.
  explicit ParticleProxy2(const ParticleProxy2<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  /// Copy assign a particle proxy.
  /// @param other The particle proxy to copy.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy2 &operator=(const ParticleProxy2 &other) noexcept = default;

  /// Copy assign a mutable particle proxy.
  /// @param other The mutable particle proxy to copy.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy2 &operator=(const ParticleProxy2<false> &other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Move assign a particle proxy.
  /// @param other The particle proxy to move.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy2 &operator=(ParticleProxy2 &&other) noexcept = default;

  /// Move assign a mutable particle proxy.
  /// @param other The mutable particle proxy to move.
  /// @return Reference to this particle proxy after assignment.
  ParticleProxy2 &operator=(ParticleProxy2<false> &&other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Returns a const proxy of the particle.
  /// @return A const proxy of the particle.
  ParticleProxy2<true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_index};
  }

  /// Gets the container holding the particle.
  /// @return A reference to the container holding the particle.
  ParticleContainer2 &container() const noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the particle.
  /// @return A const reference to the container holding the particle.
  const ParticleContainer2 &container() const noexcept { return *m_container; }
  /// Gets the index of the particle in the container.
  /// @return The index of the particle in the container.
  Index index() const noexcept { return m_index; }

  MutableGeneratorParticleProxy2 generationState() noexcept
    requires(!ReadOnly);
  MutableSimulationParticleProxy2 simulationState() noexcept
    requires(!ReadOnly);

  ConstGeneratorParticleProxy2 generationState() const noexcept;
  ConstSimulationParticleProxy2 simulationState() const noexcept;

  void assignParentIndices(std::span<const ParticleIndex2> parentIndices)
    requires(!ReadOnly)
  {
    if (!m_container->m_parentIndicesOffsetColumn.has_value() ||
        !m_container->m_parentIndicesCountColumn.has_value()) {
      throw std::logic_error("No parent indices column available");
    }
    if (accessImpl(m_container->m_parentIndicesCountColumn) != 0) {
      throw std::logic_error("Parent indices already assigned to the particle");
    }

    accessImpl(m_container->m_parentIndicesOffsetColumn) =
        static_cast<std::uint32_t>(m_container->m_parentIndices.size());
    accessImpl(m_container->m_parentIndicesCountColumn) =
        static_cast<std::uint8_t>(parentIndices.size());
    m_container->m_parentIndices.insert(m_container->m_parentIndices.end(),
                                        parentIndices.begin(),
                                        parentIndices.end());
  }

  std::span<ParticleIndex2> parentIndices() noexcept
    requires(!ReadOnly)
  {
    return {m_container->m_parentIndices.data() +
                accessImpl(m_container->m_parentIndicesOffsetColumn),
            accessImpl(m_container->m_parentIndicesCountColumn)};
  }

  Barcode &barcode() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_barcodeColumn);
  }

  Acts::PdgParticle &pdg() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_pdgColumn);
  }

  double &charge() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_chargeColumn);
  }

  double &mass() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_massColumn);
  }

  GenerationProcess &generationProcess() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_generationProcessColumn);
  }

  Acts::Vector4 &initialFourPosition() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_initialFourPositionColumn);
  }

  double &initialAbsoluteMomentum() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_initialAbsoluteMomentumColumn);
  }

  Acts::Vector3 &initialDirection() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_initialDirectionColumn);
  }

  const Acts::Surface *&currentReferenceSurface() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_currentReferenceSurfaceColumn);
  }

  Acts::Vector4 &currentFourPosition() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_currentFourPositionColumn);
  }

  double &currentAbsoluteMomentum() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_currentAbsoluteMomentumColumn);
  }

  Acts::Vector3 &currentDirection() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_currentDirectionColumn);
  }

  double &properTime() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_properTimeColumn);
  }

  double &pathInX0() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_pathInX0Column);
  }

  double &pathInL0() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_pathInL0Column);
  }

  std::uint32_t &numberOfHits() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_numberOfHitsColumn);
  }

  SimulationOutcome &simulationOutcome() noexcept
    requires(!ReadOnly)
  {
    return accessImpl(m_container->m_simulationOutcomeColumn);
  }

  std::span<const ParticleIndex2> parentIndices() const noexcept {
    return {m_container->m_parentIndices.data() +
                accessImpl(m_container->m_parentIndicesOffsetColumn),
            accessImpl(m_container->m_parentIndicesCountColumn)};
  }

  const Barcode &barcode() const noexcept {
    return accessImpl(m_container->m_barcodeColumn);
  }

  Acts::PdgParticle pdg() const noexcept {
    return accessImpl(m_container->m_pdgColumn);
  }

  double charge() const noexcept {
    return accessImpl(m_container->m_chargeColumn);
  }

  double mass() const noexcept { return accessImpl(m_container->m_massColumn); }

  GenerationProcess generationProcess() const noexcept {
    return accessImpl(m_container->m_generationProcessColumn);
  }

  const Acts::Vector4 &initialFourPosition() const noexcept {
    return accessImpl(m_container->m_initialFourPositionColumn);
  }

  double initialAbsoluteMomentum() const noexcept {
    return accessImpl(m_container->m_initialAbsoluteMomentumColumn);
  }

  const Acts::Vector3 &initialDirection() const noexcept {
    return accessImpl(m_container->m_initialDirectionColumn);
  }

  const Acts::Surface *currentReferenceSurface() const noexcept {
    return accessImpl(m_container->m_currentReferenceSurfaceColumn);
  }

  const Acts::Vector4 &currentFourPosition() const noexcept {
    return accessImpl(m_container->m_currentFourPositionColumn);
  }

  double currentAbsoluteMomentum() const noexcept {
    return accessImpl(m_container->m_currentAbsoluteMomentumColumn);
  }

  const Acts::Vector3 &currentDirection() const noexcept {
    return accessImpl(m_container->m_currentDirectionColumn);
  }

  double properTime() const noexcept {
    return accessImpl(m_container->m_properTimeColumn);
  }

  double pathInX0() const noexcept {
    return accessImpl(m_container->m_pathInX0Column);
  }

  double pathInL0() const noexcept {
    return accessImpl(m_container->m_pathInL0Column);
  }

  std::uint32_t numberOfHits() const noexcept {
    return accessImpl(m_container->m_numberOfHitsColumn);
  }

  const SimulationOutcome &simulationOutcome() const noexcept {
    return accessImpl(m_container->m_simulationOutcomeColumn);
  }

  /// Const access to an extra column of data for the space point.
  /// @param column The extra column proxy to access.
  /// @return A const reference to the value in the extra column for the space point.
  template <typename T>
  const T &extra(const ConstParticleColumnProxy2<T> &column) const noexcept {
    return column[m_index];
  }

  /// Check if this is a secondary particle.
  /// @return True if particle is a secondary (has non-zero vertex secondary, generation, or sub-particle), false otherwise
  bool isSecondary() const {
    return barcode().vertexSecondary() != 0 || barcode().generation() != 0 ||
           barcode().subParticle() != 0;
  }

  Acts::PdgParticle absolutePdg() const {
    return Acts::makeAbsolutePdgParticle(pdg());
  }

  /// Particle absolute charge.
  /// @return The absolute particle charge (positive value)
  double absoluteCharge() const { return std::abs(charge()); }

  /// Particle hypothesis.
  /// @return Particle hypothesis containing PDG, mass, and charge information
  Acts::ParticleHypothesis hypothesis() const {
    return Acts::ParticleHypothesis(
        absolutePdg(), static_cast<float>(mass()),
        Acts::AnyCharge{static_cast<float>(absoluteCharge())});
  }

 private:
  Container *m_container{nullptr};
  Index m_index{0};

  template <typename OptColumn>
  auto &accessImpl(OptColumn &&column) const {
    assert(column.has_value() && "Column does not exist");
    assert(m_index < column->size() && "Index out of bounds");
    return column->proxy(*m_container)[m_index];
  }
};

template <typename state_accessor, bool read_only>
class ParticleStateProxy2 final {
 public:
  ///
  using Index = ParticleIndex2;

  ///
  using StateAccessor = state_accessor;

  /// Indicates whether this particle proxy is read-only or data can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using Particle = ParticleProxy2<ReadOnly>;

  /// Type alias for container type (const if read-only)
  using Container = Acts::const_if_t<ReadOnly, ParticleContainer2>;

  ///
  explicit ParticleStateProxy2(Particle particle) noexcept
      : m_particle(std::move(particle)) {}

  ///
  explicit ParticleStateProxy2(ParticleProxy2<false> particle) noexcept
    requires(ReadOnly)
      : m_particle(std::move(particle)) {}

  /// Copy construct a particle proxy.
  /// @param other The particle proxy to copy.
  ParticleStateProxy2(const ParticleStateProxy2 &other) noexcept = default;

  /// Copy construct a mutable particle proxy.
  /// @param other The mutable particle proxy to copy.
  explicit ParticleStateProxy2(
      const ParticleStateProxy2<StateAccessor, false> &other) noexcept
    requires ReadOnly
      : m_particle(other.m_particle) {}

  /// Copy assign a particle proxy.
  /// @param other The particle proxy to copy.
  /// @return Reference to this particle proxy after assignment.
  ParticleStateProxy2 &operator=(const ParticleStateProxy2 &other) noexcept =
      default;

  /// Copy assign a mutable particle proxy.
  /// @param other The mutable particle proxy to copy.
  /// @return Reference to this particle proxy after assignment.
  ParticleStateProxy2 &operator=(
      const ParticleStateProxy2<StateAccessor, false> &other) noexcept
    requires ReadOnly
  {
    m_particle = other.m_particle;
    return *this;
  }

  /// Move assign a particle proxy.
  /// @param other The particle proxy to move.
  /// @return Reference to this particle proxy after assignment.
  ParticleStateProxy2 &operator=(ParticleStateProxy2 &&other) noexcept =
      default;

  /// Move assign a mutable particle proxy.
  /// @param other The mutable particle proxy to move.
  /// @return Reference to this particle proxy after assignment.
  ParticleStateProxy2 &operator=(
      ParticleStateProxy2<StateAccessor, false> &&other) noexcept
    requires ReadOnly
  {
    m_particle = other.m_particle;
    return *this;
  }

  /// Returns a const proxy of the particle.
  /// @return A const proxy of the particle.
  ParticleStateProxy2<StateAccessor, true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return ParticleStateProxy2<StateAccessor, true>(m_particle.asConst());
  }

  ///
  Particle particle() const noexcept { return m_particle; }

  std::span<const ParticleIndex2> parentIndices() const noexcept {
    return particle().parentIndices();
  }

  const Barcode &barcode() const noexcept { return particle().barcode(); }

  Acts::PdgParticle pdg() const noexcept { return particle().pdg(); }

  double charge() const noexcept { return particle().charge(); }

  double mass() const noexcept { return particle().mass(); }

  GenerationProcess generationProcess() const noexcept {
    return particle().generationProcess();
  }

  const Acts::Surface *&referenceSurface() noexcept
    requires(!ReadOnly)
  {
    return StateAccessor::referenceSurface(particle());
  }

  Acts::Vector4 &fourPosition() noexcept
    requires(!ReadOnly)
  {
    return StateAccessor::fourPosition(particle());
  }

  double &absoluteMomentum() noexcept
    requires(!ReadOnly)
  {
    return StateAccessor::absoluteMomentum(particle());
  }

  Acts::Vector3 &direction() noexcept
    requires(!ReadOnly)
  {
    return StateAccessor::direction(particle());
  }

  const Acts::Surface *referenceSurface() const noexcept {
    return StateAccessor::referenceSurface(particle());
  }

  const Acts::Vector4 &fourPosition() const noexcept {
    return StateAccessor::fourPosition(particle());
  }

  double absoluteMomentum() const noexcept {
    return StateAccessor::absoluteMomentum(particle());
  }

  const Acts::Vector3 &direction() const noexcept {
    return StateAccessor::direction(particle());
  }

  /// Total energy, i.e. norm of the four-momentum.
  /// @return The total energy calculated from mass and momentum
  double energy() const {
    return Acts::fastHypot(particle().mass(), absoluteMomentum());
  }

  /// Change the energy by the given amount.
  ///
  /// Energy loss corresponds to a negative change. If the updated energy
  /// would result in an unphysical value, the particle is put to rest, i.e.
  /// its absolute momentum is set to zero.
  /// @param delta Energy change (negative for energy loss)
  /// @return Reference to this particle for method chaining
  void applyEnergyDelta(double delta)
    requires(!ReadOnly)
  {
    const auto newEnergy = energy() + delta;
    if (particle().mass() >= newEnergy) {
      absoluteMomentum() = 0;
    } else {
      absoluteMomentum() = Acts::fastCathetus(newEnergy, particle().mass());
    }
  }
  /// Particle qOverP.
  /// @return The charge over momentum ratio
  double qOverP() const {
    return particle().hypothesis().qOverP(absoluteMomentum(),
                                          particle().charge());
  }

  /// Three-position, i.e. spatial coordinates without the time.
  /// @return Three-dimensional position vector (x, y, z)
  auto position() const {
    return fourPosition().template segment<3>(Acts::ePos0);
  }
  /// Time coordinate.
  /// @return The time coordinate value
  double time() const { return fourPosition()[Acts::eTime]; }
  /// Energy-momentum four-vector.
  /// @return Four-dimensional momentum vector (px, py, pz, E)
  Acts::Vector4 fourMomentum() const {
    Acts::Vector4 mom4;
    // stored direction is always normalized
    mom4[Acts::eMom0] = absoluteMomentum() * direction()[Acts::ePos0];
    mom4[Acts::eMom1] = absoluteMomentum() * direction()[Acts::ePos1];
    mom4[Acts::eMom2] = absoluteMomentum() * direction()[Acts::ePos2];
    mom4[Acts::eEnergy] = energy();
    return mom4;
  }
  /// Polar angle.
  /// @return The polar angle (theta) in radians
  double theta() const { return Acts::VectorHelpers::theta(direction()); }
  /// Azimuthal angle.
  /// @return The azimuthal angle (phi) in radians
  double phi() const { return Acts::VectorHelpers::phi(direction()); }
  /// Absolute momentum in the x-y plane.
  /// @return The transverse momentum magnitude
  double transverseMomentum() const {
    return absoluteMomentum() *
           direction().template segment<2>(Acts::eMom0).norm();
  }
  /// Absolute momentum.
  /// @return Three-dimensional momentum vector
  Acts::Vector3 momentum() const { return absoluteMomentum() * direction(); }

  /// Check if the particle has a reference surface.
  /// @return True if reference surface is set, false otherwise
  bool hasReferenceSurface() const { return referenceSurface() != nullptr; }

  /// Bound track parameters.
  /// @param gctx Geometry context for coordinate transformations
  /// @return Result containing bound track parameters or error if no reference surface
  Acts::Result<Acts::BoundTrackParameters> boundParameters(
      const Acts::GeometryContext &gctx) const {
    if (!hasReferenceSurface()) {
      return Acts::Result<Acts::BoundTrackParameters>::failure(
          std::error_code());
    }
    Acts::Result<Acts::Vector2> localResult =
        referenceSurface()->globalToLocal(gctx, position(), direction());
    if (!localResult.ok()) {
      return localResult.error();
    }
    Acts::BoundVector params;
    params << localResult.value(), phi(), theta(), qOverP(), time();
    return Acts::BoundTrackParameters(referenceSurface()->getSharedPtr(),
                                      params, std::nullopt,
                                      particle().hypothesis());
  }

  /// @return Curvilinear track parameters representation
  Acts::BoundTrackParameters curvilinearParameters() const {
    return Acts::BoundTrackParameters::createCurvilinear(
        fourPosition(), direction(), qOverP(), std::nullopt,
        particle().hypothesis());
  }

  /// Check if the particle is alive, i.e. is not at rest.
  /// @return True if particle has non-zero momentum, false otherwise
  bool isAlive() const { return absoluteMomentum() > 0; }

 private:
  Particle m_particle;
};

namespace detail {

struct GenerationStateAccessor final {
  static const Acts::Surface *&referenceSurface(
      MutableParticleProxy2 particle) {
    static_cast<void>(particle);
    throw std::logic_error(
        "Generation state does not have a reference surface");
  }
  static Acts::Vector4 &fourPosition(MutableParticleProxy2 particle) {
    return particle.initialFourPosition();
  }
  static double &absoluteMomentum(MutableParticleProxy2 particle) {
    return particle.initialAbsoluteMomentum();
  }
  static Acts::Vector3 &direction(MutableParticleProxy2 particle) {
    return particle.initialDirection();
  }
  static const Acts::Surface *referenceSurface(ConstParticleProxy2 particle) {
    static_cast<void>(particle);
    return nullptr;
  }
  static const Acts::Vector4 &fourPosition(ConstParticleProxy2 particle) {
    return particle.initialFourPosition();
  }
  static double absoluteMomentum(ConstParticleProxy2 particle) {
    return particle.initialAbsoluteMomentum();
  }
  static const Acts::Vector3 &direction(ConstParticleProxy2 particle) {
    return particle.initialDirection();
  }
};

struct SimulationStateAccessor final {
  static const Acts::Surface *&referenceSurface(
      MutableParticleProxy2 particle) {
    return particle.currentReferenceSurface();
  }
  static Acts::Vector4 &fourPosition(MutableParticleProxy2 particle) {
    return particle.currentFourPosition();
  }
  static double &absoluteMomentum(MutableParticleProxy2 particle) {
    return particle.currentAbsoluteMomentum();
  }
  static Acts::Vector3 &direction(MutableParticleProxy2 particle) {
    return particle.currentDirection();
  }
  static const Acts::Surface *referenceSurface(ConstParticleProxy2 particle) {
    return particle.currentReferenceSurface();
  }
  static const Acts::Vector4 &fourPosition(ConstParticleProxy2 particle) {
    return particle.currentFourPosition();
  }
  static double absoluteMomentum(ConstParticleProxy2 particle) {
    return particle.currentAbsoluteMomentum();
  }
  static const Acts::Vector3 &direction(ConstParticleProxy2 particle) {
    return particle.currentDirection();
  }
};

}  // namespace detail

/// Additional column of data that can be added to the space point container.
/// The column is indexed by the space point index.
template <typename T, bool read_only>
class ParticleColumnProxy2 final {
 public:
  /// Flag indicating whether this particle column proxy is read-only
  constexpr static bool ReadOnly = read_only;
  /// Type alias for particle index type
  using Index = ParticleIndex2;
  /// Type alias for particle index subset type
  using IndexSubset = ParticleIndexSubset2;
  /// Type alias for column value type
  using Value = T;
  /// Type alias for container type (const if read-only)
  using Container = Acts::const_if_t<ReadOnly, ParticleContainer2>;
  /// Type alias for column container type (const if read-only)
  using Column = Acts::const_if_t<ReadOnly, std::vector<Value>>;

  /// Constructs a particle column proxy for the given container and column.
  /// @param container The container holding the particle.
  /// @param column The column of data to access.
  ParticleColumnProxy2(Container &container, Column &column)
      : m_container{&container}, m_column(&column) {}

  /// Copy construct a particle column proxy.
  /// @param other The particle column proxy to copy.
  ParticleColumnProxy2(const ParticleColumnProxy2 &other) noexcept = default;

  /// Copy construct a mutable particle column proxy.
  /// @param other The mutable particle column proxy to copy.
  explicit ParticleColumnProxy2(
      const ParticleColumnProxy2<T, false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_column(&other.column()) {}

  /// Returns the number of entries in the particle column.
  /// @return The size of the particle column.
  std::uint32_t size() const noexcept { return column().size(); }

  /// Returns a const proxy of the particle column.
  /// @return A const proxy of the particle column.
  ParticleColumnProxy2<T, true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, *m_column};
  }

  /// Gets the container holding the particle.
  /// @return A reference to the container holding the particle.
  ParticleContainer2 &container() noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the particle.
  /// @return A const reference to the container holding the particle.
  const ParticleContainer2 &container() const noexcept { return *m_container; }

  /// Returns a const reference to the column container.
  /// @return A const reference to the column container.
  const std::vector<Value> &column() const noexcept { return *m_column; }

  /// Returns a mutable span to the column data.
  /// @return A mutable span to the column data.
  std::span<Value> data() noexcept
    requires(!ReadOnly)
  {
    return std::span<Value>(column().data(), column().size());
  }
  /// Returns a const span to the column data.
  /// @return A const span to the column data.
  std::span<const Value> data() const noexcept {
    return std::span<const Value>(column().data(), column().size());
  }

  /// Returns a mutable reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A mutable reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  Value &at(Index index)
    requires(!ReadOnly)
  {
    if (index >= column().size()) {
      throw std::out_of_range(
          "Index out of range in ParticleContainer: " + std::to_string(index) +
          " >= " + std::to_string(size()));
    }
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the particle to access.
  /// @return A const reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  const Value &at(Index index) const {
    if (index >= column().size()) {
      throw std::out_of_range(
          "Index out of range in ParticleContainer: " + std::to_string(index) +
          " >= " + std::to_string(size()));
    }
    return data()[index];
  }

  /// Returns a mutable reference to the column entry at the given index.
  /// @param index The index of the particle to access.
  /// @return A mutable reference to the column entry at the given index.
  Value &operator[](Index index) noexcept
    requires(!ReadOnly)
  {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// @param index The index of the particle to access.
  /// @return A const reference to the column entry at the given index.
  const Value &operator[](Index index) const noexcept {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }

  /// Subset view over selected column entries.
  class Subset final
      : public Acts::detail::ContainerSubset<Subset, Subset, Column, Value,
                                             Index, ReadOnly> {
   public:
    /// Base class type
    using Base = Acts::detail::ContainerSubset<Subset, Subset, Column, Value,
                                               Index, ReadOnly>;

    using Base::Base;
  };

  /// Creates a subset view of this particle column based on provided
  /// indices.
  ///
  /// This method creates a subset proxy that provides access to only the
  /// particles at the indices specified in the IndexSubset. The subset
  /// maintains a reference to the original column data without copying,
  /// enabling efficient access to selected particles for filtering, clustering,
  /// or other operations.
  ///
  /// @param subset The index subset specifying which particles to include
  /// @return A subset proxy providing access to the selected particles
  ///
  /// @note The returned subset shares data with the original column
  /// @note The subset remains valid only as long as the original column exists
  /// @note This operation does not copy data, providing efficient subset access
  Subset subset(const IndexSubset &subset) const noexcept {
    return Subset(*m_column, subset);
  }

 private:
  Container *m_container{};
  Column *m_column{};

  std::vector<Value> &column() noexcept
    requires(!ReadOnly)
  {
    return *m_column;
  }
};

template <bool read_only>
MutableGeneratorParticleProxy2
ParticleProxy2<read_only>::generationState() noexcept
  requires(!ReadOnly)
{
  return MutableGeneratorParticleProxy2(*this);
}

template <bool read_only>
MutableSimulationParticleProxy2
ParticleProxy2<read_only>::simulationState() noexcept
  requires(!ReadOnly)
{
  return MutableSimulationParticleProxy2(*this);
}

template <bool read_only>
ConstGeneratorParticleProxy2 ParticleProxy2<read_only>::generationState()
    const noexcept {
  return ConstGeneratorParticleProxy2(*this);
}

template <bool read_only>
ConstSimulationParticleProxy2 ParticleProxy2<read_only>::simulationState()
    const noexcept {
  return ConstSimulationParticleProxy2(*this);
}

inline MutableParticleProxy2 ParticleContainer2::at(Index index) {
  if (index >= size()) {
    throw std::out_of_range(
        "Index out of range in ParticleContainer: " + std::to_string(index) +
        " >= " + std::to_string(size()));
  }
  return MutableProxy(*this, index);
}

inline ConstParticleProxy2 ParticleContainer2::at(Index index) const {
  if (index >= size()) {
    throw std::out_of_range(
        "Index out of range in ParticleContainer: " + std::to_string(index) +
        " >= " + std::to_string(size()));
  }
  return ConstProxy(*this, index);
}

inline MutableParticleProxy2 ParticleContainer2::operator[](
    Index index) noexcept {
  return MutableProxy(*this, index);
}

inline ConstParticleProxy2 ParticleContainer2::operator[](
    Index index) const noexcept {
  return ConstProxy(*this, index);
}

}  // namespace ActsFatras

inline std::ostream &operator<<(std::ostream &os,
                                ActsFatras::ConstParticleProxy2 particle) {
  // compact format w/ only identity information but no kinematics
  os << "id=" << particle.index();
  os << "|barcode=" << "(" << particle.barcode() << ")";
  os << "|pdg=" << particle.pdg();
  os << "|q=" << particle.charge();
  os << "|m=" << particle.mass();
  os << "|p=" << particle.initialAbsoluteMomentum();
  return os;
}

inline std::ostream &operator<<(std::ostream &os,
                                ActsFatras::MutableParticleProxy2 particle) {
  return os << particle.asConst();
}
