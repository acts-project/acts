// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <vector>

#include <boost/pending/disjoint_sets.hpp>


namespace Concepts {
template <typename...> using void_t = void;

template <typename, typename T = void>
struct cell_type_has_required_functions : std::false_type {};

template <typename T>
struct cell_type_has_required_functions<T, void_t<
    decltype(get_cell_row(std::declval<T>())),
    decltype(get_cell_column(std::declval<T>())),
    decltype(get_cell_activation(std::declval<T>()))>>
    : std::true_type {};
} // namespace Concepts

template <typename T>
constexpr bool CellConcept =
    Concepts::cell_type_has_required_functions<T>();

#define CHECK_CELL_TYPE(T) static_assert(CellConcept<T>, "Cell type should have the following functions:'int get_cell_row(const Cell&)', 'int get_cell_column(const Cell&)', 'float get_cell_activation(const Cell&)'")

// TODO namespace
// TODO internal namespace

// TODO descr
template <typename Cell>
struct Compare {
    bool operator()(const Cell &c0, const Cell &c1) const {
	int row0 = get_cell_row(c0);
	int row1 = get_cell_row(c1);
	int col0 = get_cell_column(c0);
	int col1 = get_cell_column(c1);
	return (col0 == col1)? row0 < row1 : col0 < col1;
    }
};

// Definition of the Cell api for the LabeledCell<Cell> wrapper

template <typename Cell>
int get_cell_row(const LabeledCell<Cell>& cell)
{
    return get_cell_row(*cell.ptr);
}

template <typename Cell>
int get_cell_column(const LabeledCell<Cell>& cell)
{
    return get_cell_column(*cell.ptr);
}

template <typename Cell>
Label& get_cell_label(LabeledCell<Cell>& lcell)
{
    return lcell.lbl;
}

// Simple wrapper around boost::disjoint_sets. In theory, could use
// boost::vector_property_map and use boost::disjoint_sets without
// wrapping, but it's way slower
class DisjointSets {
public:
    DisjointSets(size_t initial_size = 128) :
	m_global_id(1),
	m_size(initial_size),
	m_rank(m_size),
	m_parent(m_size),
	m_ds(&m_rank[0], &m_parent[0])
	{}

    Label make_set() {
	// Empirically, m_size = 128 seems to be good default. If we
	// exceed this, take a performance hit and do the right thing.
	while (m_global_id >= m_size) {
	    m_size *= 2;
	    m_rank.resize(m_size);
	    m_parent.resize(m_size);
	    m_ds = boost::disjoint_sets<size_t*, size_t*>(&m_rank[0], &m_parent[0]);
	}
	m_ds.make_set(m_global_id);
	return static_cast<Label>(m_global_id++);
    }

    void union_set(size_t x, size_t y) {m_ds.union_set(x, y);}
    Label find_set(size_t x) {return static_cast<Label>(m_ds.find_set(x));}
    
private:
    size_t m_global_id;
    size_t m_size;
    std::vector<size_t> m_rank;
    std::vector<size_t> m_parent;
    boost::disjoint_sets<size_t*, size_t*> m_ds;
};


// TODO descr
template <typename Cell>
int get_connexions(typename std::vector<Cell>::iterator it, std::vector<Cell>& set, bool commonCorner, std::array<Label, 4>& seen)
{
    int curr_row = get_cell_row(*it);
    int curr_col = get_cell_column(*it);
    int nconn = 0;
    seen[0] = seen[1] = seen[2] = seen[3] = Label::None;
    typename std::vector<Cell>::iterator it_2 {it};

    while (it_2 != set.begin()) {
	it_2 = std::prev(it_2);
	int delta_row = std::abs(curr_row - get_cell_row(*it_2));
	int delta_col = std::abs(curr_col - get_cell_column(*it_2));

	// Iteration is column-wise, so if too far in column, can
	// safely stop
	if (delta_col > 1)
	    break;
	// For same reason, if too far in row we know the pixel is not
	// connected, but need to keep iterating
	if (delta_row > 1)
	    continue;

	// Decide whether or not cluster is connected based on 4- or
	// 8-connectivity
	if ((delta_row + delta_col) == 1 ||
	    (commonCorner and delta_row == 1 and delta_col == 1)) {
	    seen[nconn++] = get_cell_label(*it_2);
	    if (nconn == (commonCorner? 4 : 2))
		break;
	}
    }
	
    return nconn;
}

// TODO descr
// TODO should be exposed at top-level? i.e. in header
template <typename Cell, typename CellCollection>
void label_clusters_inplace(CellCollection& lcells, bool commonCorner)
{
    DisjointSets ds{};
    std::array<Label, 4> seen = {
	Label::None, Label::None, Label::None, Label::None
    };

    // Sort cells by position to enable in-order scan
    std::sort(lcells.begin(), lcells.end(), Compare<Cell>());

    for (auto it = lcells.begin(); it != lcells.end(); ++it) {
	    int nconn = get_connexions<Cell>(it, lcells, commonCorner, seen);
	    if (nconn == 0) {
		// Allocate new label
		get_cell_label(*it) = ds.make_set();
	    } else {
		// Sanity check: first element should always have
		// label if nconn > 0
		if (seen[0] == Label::None)
		    throw std::logic_error("nconn > 0 but seen[0] == Label::None");

		// Record equivalences
		for (size_t i = 1; i < (commonCorner? 4 : 2); i++) {
		    if (seen[i] != Label::None and seen[0] != seen[i]) {
			ds.union_set(seen[0], seen[i]);
		    }
		}
		// Set label for current cell
		get_cell_label(*it) = seen[0];
	    }
    }

    // Second pass: Merge labels
    for (auto& cell : lcells) {
	Label& lbl = get_cell_label(cell);
	lbl = ds.find_set(lbl);
    }
}

// TODO descr
template <typename Cell, typename InputIt>
LabeledCellCollection<Cell> labelClusters(InputIt begin, InputIt end, bool commonCorner, float threshold)
{
    CHECK_CELL_TYPE(Cell);
    LabeledCellCollection<Cell> lcells;
    for (InputIt it = begin; it != end; ++it) {
	if (get_cell_activation(*it) >= threshold)
	    lcells.push_back(LabeledCell<Cell>(*it));
    }
    label_clusters_inplace<LabeledCell<Cell>>(lcells, commonCorner);
    return lcells;
}

// TODO descr
template <typename Cell, typename ClusterT, typename InputIt, typename OutputIt>
void createClusters(InputIt begin, InputIt end, OutputIt out, bool commonCorner, float threshold)
{
    CHECK_CELL_TYPE(Cell);
    LabeledCellCollection<Cell> lcells {
	labelClusters<Cell>(begin, end, commonCorner, threshold)
    };
    if (lcells.empty())
	return;

    // Sort the cells by their cluster label
    std::sort(lcells.begin(), lcells.end(),
	      [](const LabeledCell<Cell>& lhs, const LabeledCell<Cell>& rhs)
		  {return lhs.lbl < rhs.lbl;});

    // Accumulate clusters into the output collections
    ClusterT cl;
    size_t lbl = lcells.front().lbl;
    for (auto& cell : lcells) {
	if (cell.lbl != lbl) {
	    // New cluster, save previous one
	    *out++ = std::move(cl);
	    cl = ClusterT();
	    lbl = cell.lbl;
	}
	cluster_add_cell(cl, *cell.ptr);
    }
    // Get the last cluster as well
    *out++ = std::move(cl);
}
