#pragma once

#include "traccc/device/global_index.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


TRACCC_HOST_DEVICE
inline unsigned short count_bits(unsigned long long x) {
	unsigned short nOnes = 0;
	for (int i(0); i<64; ++i) {
	   if ((x & (1ULL << i)) != 0ULL) ++nOnes;
	}
	return nOnes;
}


namespace traccc::device {

TRACCC_HOST_DEVICE
int quant(scalar min, scalar max, unsigned nSteps, scalar val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}

TRACCC_HOST_DEVICE
scalar transform(const scalar qopt_value,
                 const scalar r,
                 const scalar phi,
		 const scalar magnetic_field_z) {
    const scalar d0 = 0;
    scalar phi_value = asin(0.50f * r * magnetic_field_z * qopt_value - d0 / r) + phi;
    ::remainderf(phi_value, 2.00f * std::numbers::pi_v<scalar>);
    return phi_value;
}

template<typename T,
	unsigned int N_COLS,
	unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void ranges(const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device& container,
	    const scalar magnetic_field_z,
            const scalar r,
            const scalar phi,
            unsigned int* min_ranges,
            unsigned int* max_ranges) {
    for (unsigned int i(0); i<N_COLS; ++i) {
        scalar min_value_phi = transform(container.edges_x()[i], r, phi, magnetic_field_z);
        scalar max_value_phi = transform(container.edges_x()[i + 1], r, phi, magnetic_field_z);
        if (min_value_phi > max_value_phi) std::swap(min_value_phi, max_value_phi);
        
        // out of bounds
        if (min_value_phi < container.edges_y()[0] ||
            max_value_phi > container.edges_y()[N_ROWS]) {
            min_ranges[i] = 0;
            max_ranges[i] = 0;
            return;
        }
        
        unsigned int phi_bin_min = quant(container.edges_y()[0], 
                                         container.edges_y()[N_ROWS],
                                         N_ROWS, 
                                         min_value_phi);
        unsigned int phi_bin_max = quant(container.edges_y()[0], 
                                         container.edges_y()[N_ROWS],
                                         N_ROWS, 
                                         max_value_phi) + 1;
        
        // Clamp bins
        if (phi_bin_max > device::N_ROWS) phi_bin_max = N_ROWS;
        min_ranges[i] = phi_bin_min;
        max_ranges[i] = phi_bin_max;
    }
}


template <typename bfield_t>
TRACCC_HOST_DEVICE
void fill_hough_plane_by_column(const global_index_t globalIndex,
                const bfield_t& bfield,
                const edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
                const vecmem::data::vector_view<scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::view houghplanes_view) {
    edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);
    if (globalIndex >= device::N_COLS * houghplanes.size()) return;
    
    unsigned long long n_col = globalIndex / houghplanes.size();
    unsigned long long n_wedge = globalIndex % houghplanes.size();

    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
    const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
    const vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
    const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
    const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);

    unsigned long long column_values[device::N_ROWS] {};
    unsigned long long wedge_bit = (1ULL << n_wedge);

    for (unsigned int i(0); i<spacepoints.size(); ++i) {
	    // space point is not contributing to this wedge
	    if ((spacepoint_wedge_association[i] & wedge_bit) == 0ULL) continue;
	    const auto covfie_field_at_sp = bfield.at(spacepoints.x()[i], spacepoints.y()[i], spacepoints.z()[i]);

	    // Get the range
	    scalar min_value_phi = transform(houghplanes.edges_x()[n_col], spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);
            scalar max_value_phi = transform(houghplanes.edges_x()[n_col + 1], spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);
            if (min_value_phi > max_value_phi) std::swap(min_value_phi, max_value_phi);

            // out of bounds
            if (min_value_phi < houghplanes.edges_y()[0] ||
                max_value_phi > houghplanes.edges_y()[N_ROWS]) continue;
	    
            unsigned int phi_bin_min = quant(houghplanes.edges_y()[0],
                                         houghplanes.edges_y()[device::N_ROWS],
                                         device::N_ROWS,
                                         min_value_phi);
            unsigned int phi_bin_max = quant(houghplanes.edges_y()[0],
                                         houghplanes.edges_y()[device::N_ROWS],
                                         device::N_ROWS,
                                         max_value_phi) + 1;

            // Clamp bins
            if (phi_bin_max > device::N_ROWS) phi_bin_max = N_ROWS;
            const unsigned int diff_phi = phi_bin_max - phi_bin_min;

	    if (diff_phi < device::N_ROWS / 2) {
		    for (unsigned int idx(phi_bin_min); idx < phi_bin_max; ++idx) column_values[idx] |= wedge_bit;
	    } else {
		    for (unsigned int idx(0); idx < phi_bin_min; ++idx) column_values[idx] |= wedge_bit;
		    for (unsigned int idx(phi_bin_max); idx < device::N_ROWS; ++idx) column_values[idx] |= wedge_bit;
	    }
    } // loop on space points

    static const std::array<unsigned short, 25> min_contributions {
        5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5
    };

    // Fill the plane
    for (unsigned int i(0); i<device::N_ROWS; ++i) {
	    unsigned short value = count_bits(column_values[i]);
	    if (value < min_contributions[n_wedge]) continue;
	    houghplanes.plane()[n_wedge][globalIndex + i] = value;
    }
}


template<typename bfield_t>
TRACCC_HOST_DEVICE
void fill_hough_plane_by_bin(const global_index_t globalIndex,
                const bfield_t& bfield,
                const edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
                const vecmem::data::vector_view<scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::view houghplanes_view) {
    if (globalIndex > device::N_COLS * device::N_ROWS) return;

    edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);
    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);

    const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
    const vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
    const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
    const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);

    const unsigned int x_bin = globalIndex / device::N_COLS;
    const unsigned int y_bin  = globalIndex % device::N_ROWS;

    // given the bin, take the edges in both dimensions
    const scalar min_x = houghplanes.edges_x()[x_bin];
    const scalar max_x = houghplanes.edges_x()[x_bin + 1];
    const scalar min_y = houghplanes.edges_y()[y_bin];
    const scalar max_y = houghplanes.edges_y()[y_bin + 1];

    std::array<unsigned long long, 25> binContents {};
    for (unsigned int i(0); i<n_spacepoints; ++i) {
      if (spacepoint_wedge_association[i] == 0ULL) continue;

      const auto covfie_field_at_sp = bfield.at(spacepoints.x()[i], spacepoints.y()[i], spacepoints.z()[i]);
      scalar min_value_y = transform(min_x, spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);
      scalar max_value_y = transform(max_x, spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);

      if (min_value_y > max_y) continue;
      if (max_value_y < min_y) continue;

      for (unsigned int wedge_idx(0); wedge_idx<25; ++wedge_idx) {
          if ((spacepoint_wedge_association[i] & (1ULL << wedge_idx)) == 0ULL) continue;
          binContents[wedge_idx] |= spacepoint_detector_id[i];
      }
    }

    static const std::array<unsigned short, 25> min_contributions {
        5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5
    };
    
    for (unsigned int wedge_idx(0); wedge_idx<25; ++wedge_idx) {
        if (binContents[wedge_idx] == 0ULL) continue;
	unsigned short value = count_bits(binContents[wedge_idx]);
	// if (value < min_contributions[wedge_idx]) continue;
	houghplanes.plane()[wedge_idx][globalIndex] = value;
    }
}


template <typename bfield_t>
TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
                const bfield_t& bfield,
                const edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
                const vecmem::data::vector_view<scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                edm::houghplane_collection<unsigned short, device::N_COLS, device::N_ROWS>::view houghplanes_view) {
        // fill_hough_plane_by_column<bfield_t>(globalIndex,
	fill_hough_plane_by_bin<bfield_t>(globalIndex,
                        bfield,
                        n_spacepoints,
                        spacepoints_view,
                        spacepoint_wedge_association_view,
                        spacepoint_detector_id_view,
                        spacepoint_r_view,
                        spacepoint_phi_view,
                        houghplanes_view);
}

} // namespace
