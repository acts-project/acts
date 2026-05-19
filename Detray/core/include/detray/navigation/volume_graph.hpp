// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <queue>
#include <sstream>
#include <string>
#include <utility>

namespace detray {

/// @brief Placeholder struct for the implementation of an inspection function
/// that will be executed, when a node is visited.
template <typename node_t>
struct void_node_inspector {
  bool operator()(const node_t & /*n*/) const { return true; }
};

/// @brief Placeholder struct for an action while walking through the graph.
template <typename node_t>
struct void_actor {
  void operator()(const node_t & /*n*/,
                  const dindex_range & /*edge_range*/) const {
    /*Do nothing*/
  }
};

/// @brief Uses the geometry implementations to walk through their graph-like
/// structure breadth first.
///
/// @tparam detector_t the type of geometry we want to walk along.
/// @tparam node_inspector the type of inspection to perform when a node is
///         visited
///
/// @note The detector has to expose the volume/portal interface.
template <typename detector_t,
          typename node_inspector =
              void_node_inspector<typename detector_t::volume_type>,
          template <typename...> class vector_t = dvector>
class volume_graph {
 public:
  using geo_obj_ids = typename detector_t::geo_obj_ids;
  using volume_container_t = vector_t<typename detector_t::volume_type>;
  using mask_container_t = typename detector_t::mask_container;

  /// @brief Builds a graph node from the detector collections on the fly.
  ///
  /// The node collection builds a graph node from a volume and the mask
  /// links in its surfaces, which are half edges of the graph. Every node
  /// has its index in the detector volume container and a vector of half
  /// edges.
  struct node_generator
      : public detray::ranges::view_interface<node_generator> {
    /// A functor to retrieve the half-edges of a node
    struct node_builder {
      inline void operator()(
          const typename detector_t::surface_type &sf,
          vector_t<typename detector_t::surface_type::mask_link> &half_edges)
          const {
        half_edges.push_back(sf.mask());
      }
    };

    /// A node in the graph has an index (volume index) and a collection of
    /// edges that belong to it (mask link of every surface in the volume).
    /// One mask link is a half edge in the graph.
    struct node {
      /// Constructor from a detectors volume and surface collections
      explicit node(const tracking_volume<detector_t> &volume)
          : m_idx(volume.index()) {
        // @TODO: Remove duplicates from multiple placements of surfaces
        volume.template visit_surfaces<surface_id::e_all, node_builder>(
            m_half_edges);
      }

      /// @returns volume index of the node
      dindex index() const { return m_idx; }

      /// @returns edges of the node
      const auto &half_edges() const { return m_half_edges; }

      /// Node(volume) index
      dindex m_idx;
      /// Vector of half edges towards other volumes
      vector_t<typename detector_t::surface_type::mask_link> m_half_edges;
    };

    /// @brief Iterator over the graph nodes.
    struct iterator {
      // Iterator type defs
      using volume_iter = detray::ranges::const_iterator_t<volume_container_t>;

      using difference_type = std::iter_difference_t<volume_iter>;
      using value_type = node;
      using pointer = typename std::iterator_traits<volume_iter>::pointer;
      using reference = std::iter_reference_t<volume_iter>;
      using iterator_category =
          typename std::iterator_traits<volume_iter>::iterator_category;

      /// No default construction because of reference member
      iterator() = delete;

      /// Constructor from an iterator on the detector volume container
      /// and a reference to its surface container.
      iterator(volume_iter &&vol_itr, const detector_t &det)
          : m_vol_itr(std::move(vol_itr)), m_det(det) {}

      /// Equality operator
      bool operator==(const iterator &rhs) const {
        return m_vol_itr == rhs.m_vol_itr;
      }

      /// Dereference operator @returns a graph node
      node operator*() { return node({m_det, m_vol_itr->index()}); }

      /// Prefix increment. No postfix increment implemented
      iterator &operator++() {
        ++m_vol_itr;
        return *this;
      }

      /// Prefix decrement. No postfix decrement implemented
      iterator &operator--() {
        --m_vol_itr;
        return *this;
      }

      /// Advances iterator by @param j
      constexpr auto operator+=(const difference_type j) -> iterator & {
        m_vol_itr += j;
        return *this;
      }

      /// Advances iterator by - @param j
      constexpr auto operator-=(const difference_type j) -> iterator & {
        return *this += -j;
      }

     protected:
      /// @returns an iterator that has been advanced by @param j
      friend constexpr auto operator+(const difference_type j,
                                      const iterator &itr) -> iterator {
        return {itr.m_vol_itr + j, itr.m_det};
      }

      /// @returns an iterator that has been advanced by - @param j
      friend constexpr auto operator-(const difference_type j,
                                      const iterator &itr) -> iterator {
        return itr + -j;
      }

      /// @returns distance between two iterators
      friend constexpr auto operator-(const iterator &lhs, const iterator &rhs)
          -> difference_type {
        return lhs.m_vol_itr - rhs.m_vol_itr;
      }

     private:
      /// Iterator over the detector volume container.
      volume_iter m_vol_itr;
      /// Access to detector surfaces
      const detector_t &m_det;
    };

    /// Node iterator type
    using iterator_t = iterator;

    /// No default constructor (holds member reference)
    node_generator() = delete;

    /// Constructor from a detector containers
    node_generator(const volume_container_t &volumes, const detector_t &det)
        : m_volumes(volumes), m_det(det) {}

    /// @returns beginning of collection
    iterator begin() const { return iterator(m_volumes.begin(), m_det); }

    /// @returns end of collection
    iterator end() const { return iterator(m_volumes.end(), m_det); }

    /// @returns the number of nodes
    dindex size() const { return static_cast<dindex>(m_volumes.size()); }

    const volume_container_t &m_volumes;
    const detector_t &m_det;
  };

  /// @brief Builds graph edges from the detector mask collection on the fly.
  ///
  /// Iterates through the detectors mask store given a volume id (graph
  /// node) The graph edges are constructed for the node on the fly from the
  /// half edges in the @c _edges collection (detector masks).
  struct edge_generator
      : public detray::ranges::view_interface<edge_generator> {
    /// The link to the next volume
    using mask_edge_t = typename detector_t::surface_type::navigation_link;
    /// The link to the surfaces mask
    using mask_link_t = typename detector_t::surface_type::mask_link;

    /// An edge in the graph connects two nodes (volumes). It is constructed
    /// from a surface and its mask.
    struct edge {
      edge(const dindex volume_id, const dindex link)
          : _from(volume_id), _to(link) {}

      dindex from() const { return _from; }
      dindex to() const { return _to; }

      dindex _from;
      dindex _to;
    };

    /// Nested functor that fills the edges from a mask container
    struct edges_builder {
      /// @brief Builds the collection of graph edges for a given
      /// node.
      ///
      /// From the volume index and a mask link owned by one of the
      /// volumes surfaces a vector of graph edges is built.
      /// The mask container calls this functor and provides the correct
      /// mask group from which the surface masks can be obtained from
      /// the surfaces mask range.
      ///
      /// @param mask_group the group of masks in the mask container
      ///                   of the detector
      /// @param volume_id the index of the volume/node for which to
      ///                  build edges
      /// @param mask_range the range of masks in the group for which
      ///                   to build edges
      template <typename mask_group_t, typename mask_range_t>
      inline vector_t<edge> operator()(const mask_group_t &mask_group,
                                       const mask_range_t &mask_range,
                                       const dindex volume_id) {
        vector_t<edge> edges{};
        for (const auto &mask :
             detray::ranges::subrange(mask_group, mask_range)) {
          edges.emplace_back(volume_id, mask.volume_link());
        }
        return edges;
      }
    };

    /// No default constructor (holds member reference)
    edge_generator() = delete;

    /// Constructor from the detector masks store.
    explicit edge_generator(const typename detector_t::mask_container &masks,
                            const dindex volume_id = 0,
                            const mask_link_t mask_link = {})
        : _masks(masks), _volume_id(volume_id), _mask_link{mask_link} {
      _edges = _masks.template visit<edges_builder>(_mask_link, _volume_id);
    }

    /// @returns begging of graph edges container
    auto begin() { return _edges.begin(); }

    /// @returns end of graph edges container
    auto end() { return _edges.end(); }

    /// @return updated version of @c edge_generator
    edge_generator &operator()(dindex volume_id, const mask_link_t mask_link) {
      _volume_id = volume_id;
      _mask_link = mask_link;
      _edges.clear();
      _edges = _masks.template visit<edges_builder>(_mask_link, _volume_id);

      return *this;
    }

    const mask_container_t &_masks;
    vector_t<edge> _edges;
    dindex _volume_id;
    mask_link_t _mask_link;
  };

  // Graph nodes
  using node_type = typename node_generator::node;
  // Graph edges
  using edge_type = typename edge_generator::edge;

  /// Default constructor
  volume_graph() = delete;

  /// @brief Build from existing nodes and edges, which are provide by the
  /// detector.
  ///
  /// @param det provides: geometry volumes that become the graph nodes,
  /// surfaces which are needed to index the correct masks and the
  /// masks that link to volumes and become graph edges.
  explicit volume_graph(const detector_t &det)
      : _nodes(det.volumes(), det), _edges(det.mask_store()), _adj_matrix{0} {
    build();
  }

  /// Default destructor: we don't own anything.
  ~volume_graph() = default;

  /// @return number of nodes
  dindex n_nodes() const { return static_cast<dindex>(_nodes.size()); }

  /// @return node collection - const access. */
  const auto &nodes() const { return _nodes; }

  /// @return number of surfaces/portals in the geometry */
  dindex n_edges(dindex volume_id) const {
    const node_type n = _nodes[volume_id];
    // Number of half_edges is equal to number of edges for volume
    return static_cast<dindex>(n.half_edges().size());
  }

  /// @return edges collection - const access.
  const auto &edges() const { return _edges; }

  /// @return graph adjacency - const access.
  const auto &adjacency_matrix() const { return _adj_matrix; }

  /// Walks breadth first through the geometry objects.
  /*template <typename action_t = void_actor<node_type>>
  void bfs(action_t actor = {}) const {
      // Do node inspection
      const auto inspector = node_inspector();

      node_type const *current = nullptr;
      vector_t<bool> visited(_nodes.size(), false);

      // Nodes to be visited. Start at first node
      std::queue<node_type const *> node_queue;
      node first_node = _nodes.front();
      node_queue.push(&(first_node));

      // Visit adjacent nodes and check current one
      while (!node_queue.empty()) {
          // Inspect
          current = node_queue.front();
          if (visited[current->index()]) {
              node_queue.pop();
              continue;
          }
          inspector(*current);

          // Visit neighbors and perform action
          actor(*current, current->template
  range<static_cast<geo_obj_ids>(0)>());

          // Add neighbors to queue
          for (const auto &edg_link : current->edges()) {
              for (const auto &edg :
  edge_generator::iterator(current->index(), edg_link, _edges)) { dindex nbr
  = edg.to();
                  // If not leaving world and if not visited, enqueue the node
                  if ((nbr != dindex_invalid && nbr > 0) && not
  visited[nbr]) { node_queue.push(&(_nodes[nbr]));
                  }
              }
          }

          // This node has been visited
          visited[current->index()] = true;
          node_queue.pop();
      }
  }*/

  /// @returns the linking description as a string.
  inline std::string to_string() const {
    std::stringstream stream;
    dindex dim = n_nodes() + 1u;
    for (const auto &n : _nodes) {
      stream << "[>>] Node with index " << n.index() << std::endl;
      stream << " -> edges: " << std::endl;
      for (unsigned int i = 0u; i < dim; ++i) {
        const dindex degr = _adj_matrix[dim * n.index() + i];
        if (degr == 0) {
          continue;
        }
        std::string n_occur =
            degr > 1 ? "\t\t\t\t(" + std::to_string(degr) + "x)" : "";

        // Edge that leads out of the detector world
        if (i == dim - 1u && degr != 0u) {
          stream << "    -> leaving world " + n_occur << std::endl;
        } else {
          stream << "    -> " << std::to_string(i) + "\t" + n_occur
                 << std::endl;
        }
      }
    }
    return stream.str();
  }

  /// @returns the linking description as a string in DOT syntax.
  inline std::string to_dot_string() const {
    std::stringstream stream;
    dindex dim = n_nodes() + 1u;

    stream << "strict graph {" << std::endl;
    stream << "    layout=neato;" << std::endl;
    stream << "    node [margin=0,shape=circle,style=filled];" << std::endl;
    stream << "    overlap=scale;" << std::endl;
    stream << "    splines=true;" << std::endl;
    stream << "    mode=KK;" << std::endl;
    stream << std::endl;
    stream << R"(    exit [label="OOB",fillcolor="firebrick1"];)" << std::endl;

    for (const auto &n : _nodes) {
      stream << "    v" << n.index() << ";" << std::endl;
    }

    stream << std::endl;

    for (const auto &n : _nodes) {
      for (unsigned int i = 0u; i < dim; ++i) {
        const dindex degr = _adj_matrix[dim * n.index() + i];
        if (degr == 0) {
          continue;
        }

        bool to_oob = i == dim - 1u;
        bool to_self = n.index() == i;

        std::string label = degr > 1 ? std::to_string(degr) : "";
        std::string dest = to_oob ? "exit" : ("v" + std::to_string(i));
        std::string dir = (to_oob || to_self) ? "forward" : "both";

        stream << "    v" << n.index() << " -- " << dest << " [label=\""
               << label << "\", dir=" << dir << "];" << std::endl;
      }
    }

    stream << "}" << std::endl;

    return stream.str();
  }

 private:
  /// @brief Go through the nodes and fill adjacency matrix.
  ///
  /// Root node is always at zero.
  void build() {
    // Leave space for the world volume (links to dindex_invalid)
    const dindex dim{n_nodes() + 1u};
    _adj_matrix.resize(dim * dim);

    for (const auto &n : _nodes) {
      // Only works for non batched geometries
      for (const auto &edg_link : n.half_edges()) {
        // Build an edge
        for (const auto edg : _edges(n.index(), edg_link)) {
          const dindex elem =
              edg.to() < detail::invalid_value<
                             typename edge_generator::mask_edge_t>()
                  ? dim * n.index() + edg.to()
                  : dim * n.index() + dim - 1u;
          _adj_matrix[elem]++;
        }
      }
    }
  }

  /// Graph nodes
  node_generator _nodes;

  /// Graph edges
  edge_generator _edges;

  /// Adjacency matrix
  vector_t<dindex> _adj_matrix;
};

}  // namespace detray
