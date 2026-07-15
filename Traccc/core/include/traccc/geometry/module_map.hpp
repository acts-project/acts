/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

namespace traccc {
/**
 * @brief A fast and GPU-friendly map for consecutive module IDs to values.
 *
 * Using std::map to model the relationship between geometry identifiers and
 * transformation matrices works only on the CPU. std::map is not supported at
 * all on CUDA, and either doesn't work or is inefficient on other platforms.
 *
 * This class is a more efficient and more GPU-friendly implementation of the
 * geometry ID to transformation mapping concept.
 *
 * The implementation of this class relies on the idea that module IDs exist in
 * contiguous sub-ranges, and we use this fact to get more efficient accesses.
 * First, we reduce the range of module IDs to a set of contiguous sub-ranges.
 * In the TrackML detector, there are roughly 18,000 modules, but there are
 * only 71 contiguous sub-arrays. Because we can index a contiguous range in
 * constant time, this presents possible speed-up.
 *
 * Thus, the map is essentially two-layered. First, we determine which
 * contigous sub-range a given value must fall in, and then we get the value
 * from that sub-range in constant time.
 *
 * The first layer is implemented as a simple balanced binary search tree,
 * which we implement in a GPU-friendly fashion. This means that the cost of
 * finding the relevant sub-range is O(log_2(n)), with n the number of
 * sub-ranges. For the TrackML detector, this means we can find our value in
 * only 7 steps!
 *
 * @note In theory, this class can be used for any mapping of key to value, but
 * it relies quite heavily (for performance) on the contiguous nature of the
 * keys inserted into it.
 *
 * @note This map supports only construction. You can't update it in any way.
 * To build one, construct a std::map first, and then convert it.
 */
template <typename K = geometry_id, typename V = transform3>
class module_map {
    public:
    // Default constructor
    module_map() = default;

    /**
     * @brief Construct a module map from one created by the `io` library.
     *
     * This method constructs a module map from an existing module map
     * represented as a std::map.
     *
     * @param[in] input The existing module map to convert.
     */
    module_map(const std::map<K, V>& input) {
        /*
         * First, construct a list of nodes and values, which are linked
         * together.
         */
        std::vector<module_map_node> nodes;
        std::tie(m_values, nodes) = node_list(input);

        /*
         * Lay out the nodes in memory in an efficient way.
         */
        create_tree(nodes);
    }

    /**
     * @brief Find a given key in the map.
     *
     * @param[in] i The key to look-up.
     *
     * @return The value associated with the given key.
     *
     * @warning This method does no bounds checking, and will result in
     * undefined behaviour if the key does not exist in the map.
     */
    const V& operator[](const K& i) const { return *at_helper(i, 0); }

    /**
     * @brief Find a given key in the map, with bounds checking.
     *
     * @param[in] i The key to look-up.
     *
     * @return The value associated with the given key.
     */
    const V& at(const K& i) const {
        const V* r = at_helper(i, 0);

        if (r == nullptr) {
            throw std::out_of_range("Index not found in module map!");
        }

        return *r;
    }

    /**
     * @brief Get the total number of modules in the module map.
     *
     * This iterates over all of the nodes in the map and sums up their sizes.
     *
     * @return The total number of modules in this module map.
     */
    std::size_t size(void) const {
        std::size_t c = 0;

        for (const module_map_node& n : m_nodes) {
            c += n.size;
        }

        return c;
    }

    bool contains(const K& i) const { return at_helper(i, 0) != nullptr; }

    bool empty(void) const { return m_nodes.empty(); }

    private:
    /**
     * @brief The internal representation of nodes in our binary search tree.
     *
     * These objects carry three pieces of data. Firstly, there is the starting
     * ID. Then, there is the size. Since the node represents a stretch of
     * consecutive IDs, we know that the node ends at `start + size`. Finally,
     * there is the index in the value array. We keep indices instead of
     * pointers to make it easier to port this code to other devices.
     */
    struct module_map_node {
        module_map_node() = default;

        module_map_node(K s, std::size_t n, std::size_t i)
            : start(s), size(n), index(i) {}

        K start = 0;
        std::size_t size = 0;
        std::size_t index = 0;
    };

    /**
     * @brief Lay out a set of nodes in a binary tree format.
     *
     * This method updates the binary tree contained in this object in an
     * in-place fashion.
     *
     * @param[in] nodes The nodes to create a tree from.
     */
    void create_tree(const std::vector<module_map_node>& nodes) {
        /*
         * First, we resize the number of nodes in the tree to the number of
         * nodes in the input, but rounded up such that we have we end up with
         * 2^n - 1 spaces, where n is the ceiling of log2 of the size of the
         * input, plus one.
         */
        m_nodes.resize(round_up(nodes.size()));

        /*
         * This vector is our work queue for the pre-order insertion that we
         * will be doing. Each of these values contains a position where to
         * insert the next node, and the subset of nodes to insert in that
         * sub-tree.
         */
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> z;

        /*
         * To kick-start the algorithm, we want to insert all nodes (from 0 to
         * the last one) at the first spot in the array.
         */
        z.emplace_back(0, 0, nodes.size());

        /*
         * We are done when the work queue is empty, but we will keep adding
         * things to it for a while at the start.
         */
        while (!z.empty()) {
            /*
             * Pop the current index and node range off the stack.
             */
            std::size_t i, s, e;
            std::tie(i, s, e) = z.back();
            z.pop_back();

            /*
             * Pick a node roughly in the middle of the node range. This will
             * be our next child.
             */
            std::size_t m = (s + e) / 2;

            /*
             * Insert the node we selected as our center in the given position.
             */
            m_nodes[i] = nodes[m];

            /*
             * If we have any nodes on the left side, we will put an element in
             * the work queue to place the left nodes at index 2i + 1.
             */
            if (m - s > 0) {
                z.emplace_back(2 * i + 1, s, m);
            }

            /*
             * Same thing for the right side, but we insert the right side
             * starting at 2i + 2 instead.
             */
            if (e - (m + 1) > 0) {
                z.emplace_back(2 * i + 2, m + 1, e);
            }
        }
    }

    /**
     * @brief Special rounding function for binary trees.
     *
     * This function computes 2^n - 1, where n is the ceiling of log2 of the
     * input plus one.
     *
     * This sounds like an odd operation, but consider that a binary tree (when
     * complete) of depth k, can contain at most 2^k - 1 nodes. That is the
     * number of spaces we need. However, our bottom level may have some empty
     * slots, because it is rather unlikely that we will have exactly 2^k - 1
     * nodes. Thus, we need to round up k from log2 of the input. However, we
     * need to also keep in mind that 2^l (for some integer l) should be
     * rounded up to 2^(l+1), not 2^l! Since we are computing l, this would
     * mean we would be short one space in some edge cases.
     *
     * @param[in] i The number to round up.
     *
     * @return The specially rounded number.
     */
    static std::size_t round_up(std::size_t i) {
        std::size_t n = 1;

        while (n < i + 1) {
            n *= 2;
        }

        return n - 1;
    }

    /**
     * @brief Get a list of values and nodes from an existing std::map.
     *
     * Given an existing map, this function will attempt to find all the nodes,
     * that is to say contiguous sub-ranges. It will also put the values into
     * one big, contiguous array.
     *
     * @param[in] input An existing std::map module map.
     *
     * @return A pair of vectors, one has the values, and the other has the
     * nodes.
     */
    static std::pair<std::vector<V>, std::vector<module_map_node>> node_list(
        const std::map<K, V>& input) {
        /*
         * We start out by grabbing a list of all the keys in the map.
         */
        std::vector<K> keys;
        keys.reserve(input.size());
        for (const std::pair<const K, V>& i : input) {
            keys.push_back(i.first);
        }

        /*
         * Next, we declare a vector of values, which is one of the two things
         * that we will return in the end.
         */
        std::vector<V> values;
        values.reserve(keys.size());

        /*
         * To find contiguous sub-ranges, the keys must be sorted. In general,
         * it is likely that they will be pre-sorted.
         */
        std::sort(keys.begin(), keys.end());

        /*
         * This vector will hold our nodes. We'll update it as we go, and we
         * start out with an incomplete first node, which represents the single
         * value at the start of the key vector. We will add to this later!
         */
        std::vector<module_map_node> nodes;
        nodes.push_back(module_map_node(keys.front(), 1, values.size()));
        values.push_back(input.at(keys.front()));

        /*
         * Starting from the second key, we check if they are contiuous with
         * the last key, and store them if so.
         */
        for (std::size_t i = 1; i < keys.size(); ++i) {
            /*
             * If the key is not adjacent to the last one, we push a new
             * partial node onto the node vector.
             */
            if (keys.at(i) != keys.at(i - 1) + 1) {
                nodes.push_back(module_map_node(keys.at(i), 0, values.size()));
            }

            /*
             * Push the current value into the value array, and register one
             * additional element in the current node.
             */
            ++nodes.back().size;
            values.push_back(input.at(keys.at(i)));
        }

        /*
         * Return both the generated values.
         */
        return {values, nodes};
    }

    /**
     * @brief Helper function to retrieve a value from the map.
     *
     * This function searches recursively for a value in a given sub-tree. Due
     * to the formalism used for node indices, we can identify a subtree simply
     * by the node index of its root node.
     *
     * @param[in] i The value to look for.
     * @param[in] n The index of the subtree's root node.
     */
    const V* at_helper(const K& i, std::size_t n) const {
        /*
         * For memory safety, if we are out of bounds we will exit.
         */
        if (n >= m_nodes.size()) {
            return nullptr;
        }

        /*
         * Retrieve the current root node.
         */
        const module_map_node& node = m_nodes[n];

        /*
         * If the size is zero, it is essentially an invalid node (i.e. the
         * node does not exist).
         */
        if (node.size == 0) {
            return nullptr;
        }

        /*
         * If the value we are looking for is past the start of the current
         * node, there are three possibilities. Firstly, the value might be in
         * the current node. Secondly, the value might be in the right child of
         * the current node. Thirdly, the value might not be in the map at all.
         */
        if (i >= node.start) {
            /*
             * Next, we check if the value is within the range represented by
             * the current node.
             */
            if (i < node.start + node.size) {
                /*
                 * Found it! Return a pointer to the value within the
                 * contiguous range.
                 */
                return &m_values[node.index + (i - node.start)];
            } else {
                /*
                 * Two possibilties remain, we need to check the right subtree.
                 */
                return at_helper(i, 2 * n + 2);
            }
        }
        /*
         * If the value we want to find is less then the start of this node,
         * there are only two possibilities. Firstly, the value might be in the
         * left subtree, or the value might not be in the map at all.
         */
        else {
            return at_helper(i, 2 * n + 1);
        }
    }

    /**
     * @brief The internal storage of the nodes in our binary search tree.
     *
     * This follows the well-known formalism where the root node resides at
     * index 0, while for any node at position n, the left child is at index 2n
     * + 1, and the right child is at index 2n + 2.
     */
    std::vector<module_map_node> m_nodes;

    /**
     * @brief This vector stores the values in a contiguous manner. Our nodes
     * keep indices in this array instead of pointers.
     */
    std::vector<V> m_values;
};
}  // namespace traccc
