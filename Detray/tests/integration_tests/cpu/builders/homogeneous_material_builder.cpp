// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/builders/homogeneous_material_builder.hpp"

#include "detray/builders/cuboid_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>
#include <memory>
#include <vector>

using namespace detray;

namespace {

using scalar = test::scalar;
using point2 = test::point2;
using point3 = test::point3;

using metadata_t = test::default_metadata;
using detector_builder_t = detector_builder<metadata_t>;
using detector_t = typename detector_builder_t::detector_type;
using transform3 = dtransform3D<typename detector_t::algebra_type>;
using mat_id = typename detector_t::material::id;
using bin_index_t = axis::multi_bin<2u>;

using rectangle_factory = surface_factory<detector_t, rectangle2D>;
using hom_mat_factory_t = homogeneous_material_factory<detector_t>;
using map_factory_t = material_map_factory<detector_t, bin_index_t>;

constexpr scalar tol{std::numeric_limits<scalar>::epsilon()};

constexpr dindex n_volumes{3u};

/// Volume local surface indices
/// @{
// Homogeneous material: identical in all volumes
constexpr dindex shared_slab_sf{0u};
// Homogeneous material: different in every volume
constexpr dindex unique_slab_sf{1u};
/// @}

/// Build a detector with @c n_volumes cuboid volumes, in which some surfaces
/// carry identical material
detector_t build_detector(vecmem::memory_resource &mr,
                          const volume_builder_options &builder_opts) {
  detector_builder_t det_builder{};
  det_builder.set_name("test_detector");

  for (dindex v_idx = 0u; v_idx < n_volumes; ++v_idx) {
    auto vbuilder = det_builder.new_volume(volume_id::e_cuboid);
    vbuilder->add_volume_placement(point3{
        0.f, 0.f, static_cast<scalar>(v_idx) * 100.f * unit<scalar>::mm});

    auto mat_builder =
        det_builder.template decorate<homogeneous_material_builder<detector_t>>(
            v_idx);

    // Surfaces with homogeneous material
    auto hom_factory = std::make_shared<hom_mat_factory_t>(
        std::make_unique<rectangle_factory>());

    const auto vol_link{
        static_cast<typename detector_t::surface_type::navigation_link>(v_idx)};

    hom_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -20.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    hom_factory->add_material(mat_id::e_material_slab,
                              {1.f * unit<scalar>::mm, silicon<scalar>()});

    hom_factory->push_back({surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -10.f}), vol_link,
                            std::vector<scalar>{10.f, 10.f}});
    hom_factory->add_material(
        mat_id::e_material_slab,
        {(1.f + static_cast<scalar>(v_idx)) * unit<scalar>::mm,
         tungsten<scalar>()});

    // Add a portal box around the volume
    auto portal_generator =
        std::make_shared<cuboid_portal_generator<detector_t>>(0.1f *
                                                              unit<scalar>::mm);

    mat_builder->add_surfaces(hom_factory);
    mat_builder->add_surfaces(portal_generator);
  }

  return det_builder.build(mr, builder_opts);
}

/// @returns the material slab of the surface material at a local point
struct get_material_slab {
  template <typename mat_coll_t, concepts::index index_t>
  auto operator()(const mat_coll_t &mat_coll, const index_t idx,
                  const point2 &loc_p) const -> material_slab<scalar> {
    using material_t = typename mat_coll_t::value_type;

    if constexpr (std::same_as<material_t, material_slab<scalar>>) {
      return detail::material_accessor::get(mat_coll, idx, loc_p);
    } else {
      return {};
    }
  }
};

}  // anonymous namespace

/// Integration test: material builder as volume builder decorator
GTEST_TEST(detray_builders, decorator_homogeneous_material_builder) {
  using transform3 = typename detector_t::transform3_type;
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;

  using pt_cylinder_t = concentric_cylinder2D;
  using pt_cylinder_factory_t = surface_factory<detector_t, pt_cylinder_t>;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;
  using cylinder_factory = surface_factory<detector_t, cylinder2D>;

  using mat_factory_t = homogeneous_material_factory<detector_t>;

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  volume_builder_options builder_opts{};
  builder_opts.deduplicate(false);

  auto vbuilder =
      std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);
  auto mat_builder =
      homogeneous_material_builder<detector_t>{std::move(vbuilder)};

  EXPECT_TRUE(d.volumes().empty());

  // Add some portals first
  auto pt_cyl_factory = std::make_unique<pt_cylinder_factory_t>();

  pt_cyl_factory->push_back({surface_id::e_portal,
                             transform3(point3{0.f, 0.f, 0.f}), 1u,
                             std::vector<scalar>{10.f, -1500.f, 1500.f}});
  pt_cyl_factory->push_back({surface_id::e_portal,
                             transform3(point3{0.f, 0.f, 0.f}), 2u,
                             std::vector<scalar>{20.f, -1500.f, 1500.f}});

  // Then some passive and sensitive surfaces
  auto rect_factory = std::make_unique<rectangle_factory>();

  typename rectangle_factory::sf_data_collection rect_sf_data;
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -10.f}), 0u,
                            std::vector<scalar>{10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -20.f}), 0u,
                            std::vector<scalar>{10.f, 8.f});
  rect_sf_data.emplace_back(surface_id::e_sensitive,
                            transform3(point3{0.f, 0.f, -30.f}), 0u,
                            std::vector<scalar>{10.f, 8.f});
  rect_factory->push_back(std::move(rect_sf_data));

  auto trpz_factory = std::make_unique<trapezoid_factory>();

  trpz_factory->push_back({surface_id::e_sensitive,
                           transform3(point3{0.f, 0.f, 1000.f}), 0u,
                           std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}});

  auto cyl_factory = std::make_unique<cylinder_factory>();

  cyl_factory->push_back({surface_id::e_passive,
                          transform3(point3{0.f, 0.f, 0.f}), 0u,
                          std::vector<scalar>{5.f, -1300.f, 1300.f}});

  // Now add the material for each surface
  auto mat_pt_cyl_factory =
      std::make_shared<mat_factory_t>(std::move(pt_cyl_factory));
  mat_pt_cyl_factory->add_material(material_id::e_material_slab,
                                   {1.f * unit<scalar>::mm, silicon<scalar>()});
  mat_pt_cyl_factory->add_material(
      material_id::e_material_slab,
      {1.5f * unit<scalar>::mm, silicon<scalar>()});

  auto mat_rect_factory =
      std::make_shared<mat_factory_t>(std::move(rect_factory));
  mat_rect_factory->add_material(material_id::e_material_slab,
                                 {1.f * unit<scalar>::mm, silicon<scalar>()});
  mat_rect_factory->add_material(material_id::e_material_slab,
                                 {2.f * unit<scalar>::mm, silicon<scalar>()});
  mat_rect_factory->add_material(material_id::e_material_slab,
                                 {3.f * unit<scalar>::mm, silicon<scalar>()});

  auto mat_trpz_factory =
      std::make_shared<mat_factory_t>(std::move(trpz_factory));
  mat_trpz_factory->add_material(material_id::e_material_slab,
                                 {1.f * unit<scalar>::mm, silicon<scalar>()});

  auto mat_cyl_factory =
      std::make_shared<mat_factory_t>(std::move(cyl_factory));
  mat_cyl_factory->add_material(material_id::e_material_slab,
                                {1.5f * unit<scalar>::mm, tungsten<scalar>()});

  // Add surfaces and material to detector
  mat_builder.add_surfaces(mat_pt_cyl_factory, geo_ctx);
  mat_builder.add_surfaces(mat_rect_factory, geo_ctx);
  mat_builder.add_surfaces(mat_trpz_factory, geo_ctx);
  mat_builder.add_surfaces(mat_cyl_factory, geo_ctx);

  // Add the volume to the detector
  mat_builder.build(d, builder_opts);

  //
  // check results
  //
  const auto &vol = d.volumes().back();
  EXPECT_TRUE(d.volumes().size() == 1u);
  EXPECT_EQ(vol.index(), 0u);
  EXPECT_EQ(vol.id(), volume_id::e_cylinder);

  EXPECT_EQ(d.surfaces().size(), 7u);
  EXPECT_EQ(d.transform_store().size(), 8u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_concentric_cylinder2D>(),
            2u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2D>(), 0u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2D>(), 1u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 1u);

  EXPECT_EQ(d.material_store().template size<material_id::e_material_slab>(),
            7u);
  EXPECT_EQ(d.material_store().template size<material_id::e_material_rod>(),
            0u);

  for (auto [idx, sf_desc] : detray::views::enumerate(d.surfaces())) {
    const auto &mat_link = sf_desc.material();
    EXPECT_EQ(mat_link.id(), material_id::e_material_slab);
    EXPECT_EQ(mat_link.index(), idx);
  }

  for (const auto &mat_slab :
       d.material_store().template get<material_id::e_material_slab>()) {
    EXPECT_TRUE(mat_slab.get_material() == silicon<scalar>() ||
                mat_slab.get_material() == tungsten<scalar>());
  }
}

/// Integration test: homogeneous material on a sparse subset of surfaces
GTEST_TEST(detray_builders, homogeneous_material_on_sparse_surfaces) {
  using transform3 = typename detector_t::transform3_type;
  using material_id = typename detector_t::material::id;

  using rectangle_factory = surface_factory<detector_t, rectangle2D>;
  using mat_factory_t = homogeneous_material_factory<detector_t>;

  vecmem::host_memory_resource host_mr;
  detector_t d(host_mr);
  auto geo_ctx = typename detector_t::geometry_context{};
  const volume_builder_options builder_opts{};

  // Build a dummy volume first, so that the volume under test starts neither
  // at surface nor at material index zero of the detector containers: the
  // factory indices are volume local and have to be translated accordingly
  constexpr std::size_t n_dummy_surfaces{3u};
  constexpr std::size_t n_dummy_slabs{1u};
  {
    auto dummy_vbuilder =
        std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);
    auto dummy_mat_builder =
        homogeneous_material_builder<detector_t>{std::move(dummy_vbuilder)};

    auto dummy_factory = std::make_shared<rectangle_factory>();
    typename rectangle_factory::sf_data_collection dummy_sf_data;
    for (std::size_t i = 0u; i < n_dummy_surfaces; ++i) {
      dummy_sf_data.emplace_back(
          surface_id::e_sensitive,
          transform3(point3{0.f, 0.f, 10.f * static_cast<scalar>(i)}), 0u,
          std::vector<scalar>{10.f, 8.f});
    }
    dummy_factory->push_back(std::move(dummy_sf_data));
    dummy_mat_builder.add_surfaces(dummy_factory, geo_ctx);

    // Material on the last surface of the dummy volume only
    auto dummy_mat_factory = std::make_shared<mat_factory_t>();
    dummy_mat_factory->add_material(
        material_id::e_material_slab,
        {3.f * unit<scalar>::mm, tungsten<scalar>(), n_dummy_surfaces - 1u});
    dummy_mat_builder.add_surfaces(dummy_mat_factory, geo_ctx);

    dummy_mat_builder.build(d, builder_opts);
  }

  auto vbuilder =
      std::make_unique<volume_builder<detector_t>>(volume_id::e_cylinder);
  auto mat_builder =
      homogeneous_material_builder<detector_t>{std::move(vbuilder)};

  // Five sensitive surfaces, none of which carry material yet
  constexpr std::size_t n_surfaces{5u};

  auto rect_factory = std::make_shared<rectangle_factory>();
  typename rectangle_factory::sf_data_collection rect_sf_data;
  for (std::size_t i = 0u; i < n_surfaces; ++i) {
    rect_sf_data.emplace_back(
        surface_id::e_sensitive,
        transform3(point3{0.f, 0.f, -10.f * static_cast<scalar>(i)}), 0u,
        std::vector<scalar>{10.f, 8.f});
  }
  rect_factory->push_back(std::move(rect_sf_data));
  mat_builder.add_surfaces(rect_factory, geo_ctx);

  // Attach material to surfaces 1 and 4 only: there is a gap in between, the
  // block does not start at the first surface and the material is passed in
  // with explicit surface indices
  auto mat_factory = std::make_shared<mat_factory_t>();
  mat_factory->add_material(material_id::e_material_slab,
                            {1.f * unit<scalar>::mm, silicon<scalar>(), 1u});
  mat_factory->add_material(material_id::e_material_slab,
                            {2.f * unit<scalar>::mm, tungsten<scalar>(), 4u});
  mat_builder.add_surfaces(mat_factory, geo_ctx);

  mat_builder.build(d, builder_opts);

  // One slab per material entry: the gaps must not be padded with filler
  EXPECT_EQ(d.volumes().size(), 2u);
  EXPECT_EQ(d.surfaces().size(), n_dummy_surfaces + n_surfaces);
  EXPECT_EQ(d.material_store().template size<material_id::e_material_slab>(),
            n_dummy_slabs + 2u);

  const auto &slabs =
      d.material_store().template get<material_id::e_material_slab>();

  // The dummy volume is untouched by the second volume's material
  const auto &dummy_sf = d.surface(static_cast<dindex>(n_dummy_surfaces - 1u));
  ASSERT_TRUE(dummy_sf.has_material());
  EXPECT_EQ(slabs.at(dummy_sf.material().index()).get_material(),
            tungsten<scalar>());
  EXPECT_NEAR(slabs.at(dummy_sf.material().index()).thickness(),
              3.f * unit<scalar>::mm, tol);

  for (const auto [idx, sf_desc] : detray::views::enumerate(d.surfaces())) {
    if (idx < n_dummy_surfaces) {
      continue;
    }

    // The material was configured with surface indices local to the volume
    const std::size_t sf_idx{idx - n_dummy_surfaces};

    if (sf_idx == 1u || sf_idx == 4u) {
      ASSERT_TRUE(sf_desc.has_material());
      ASSERT_EQ(sf_desc.material().id(), material_id::e_material_slab);

      const auto &slab = slabs.at(sf_desc.material().index());
      if (sf_idx == 1u) {
        EXPECT_EQ(slab.get_material(), silicon<scalar>());
        EXPECT_NEAR(slab.thickness(), 1.f * unit<scalar>::mm, tol);
      } else {
        EXPECT_EQ(slab.get_material(), tungsten<scalar>());
        EXPECT_NEAR(slab.thickness(), 2.f * unit<scalar>::mm, tol);
      }
    } else {
      // Surfaces that were not given material must not have any
      EXPECT_FALSE(sf_desc.has_material());
    }
  }
}

/// Integration test to build an empty cuboid volume with material
GTEST_TEST(detray_builders, detector_builder_with_material) {
  using namespace detray;

  using transform3 = typename detector_t::transform3_type;
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;

  // Surface factories
  using trapezoid_factory = surface_factory<detector_t, trapezoid2D>;

  // detector builder
  detector_builder<typename detector_t::metadata> det_builder{};
  auto geo_ctx = typename detector_t::geometry_context{};
  const volume_builder_options builder_opts{};

  // Vanilla volume builder
  auto vbuilder = det_builder.new_volume(volume_id::e_cuboid);
  const auto vol_idx{
      static_cast<typename detector_t::surface_type::navigation_link>(
          vbuilder->vol_index())};

  // Add material
  auto mv_builder =
      det_builder.template decorate<homogeneous_material_builder<detector_t>>(
          vbuilder->vol_index());

  assert(mv_builder != nullptr);

  typename detector_t::point3_type t{0.f, 0.f, 20.f};
  mv_builder->add_volume_placement(t);

  // Add a sensitive surface
  auto trpz_factory = std::make_unique<trapezoid_factory>();
  // Add material to the surface
  auto mat_sf_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::move(trpz_factory));

  mat_sf_factory->push_back({surface_id::e_sensitive,
                             transform3(point3{0.f, 0.f, 1000.f}), vol_idx,
                             std::vector<scalar>{1.f, 3.f, 2.f, 0.25f}});
  mat_sf_factory->add_material(material_id::e_material_slab,
                               {1.f * unit<scalar>::mm, silicon<scalar>()});

  // Add a portal box around the cuboid volume with a min distance of 'env'
  constexpr auto env{0.1f * unit<scalar>::mm};
  auto portal_generator =
      std::make_unique<cuboid_portal_generator<detector_t>>(env);

  // Add homogeneous material to every portal
  auto mat_portal_factory =
      std::make_shared<homogeneous_material_factory<detector_t>>(
          std::move(portal_generator));
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {2.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {3.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {4.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {5.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {6.f * unit<scalar>::mm, silicon<scalar>()});
  mat_portal_factory->add_material(material_id::e_material_slab,
                                   {7.f * unit<scalar>::mm, silicon<scalar>()});

  mv_builder->add_surfaces(mat_sf_factory, geo_ctx);
  mv_builder->add_surfaces(mat_portal_factory);

  //
  // build the detector
  //
  vecmem::host_memory_resource host_mr;
  const detector_t d = det_builder.build(host_mr, builder_opts);
  const auto vol = tracking_volume{d, 0u};

  // check the results
  EXPECT_EQ(d.volumes().size(), 1u);
  EXPECT_EQ(vol.id(), volume_id::e_cuboid);
  EXPECT_EQ(vol.index(), 0u);

  // Check the volume placement
  typename detector_t::transform3_type trf{t};
  EXPECT_TRUE(vol.transform() == trf);
  EXPECT_TRUE(d.transform_store().at(0u) == trf);

  EXPECT_EQ(d.surfaces().size(), 7u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2D>(), 3u);
  EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2D>(), 1u);
  EXPECT_EQ(d.material_store().template size<material_id::e_material_slab>(),
            7u);

  // Check the material links
  for (const auto [idx, sf_desc] : detray::views::enumerate(d.surfaces())) {
    EXPECT_EQ(sf_desc.material().id(), material_id::e_material_slab);
    EXPECT_EQ(sf_desc.material().index(), idx);
  }

  // Check the material
  scalar thickness{1.f * unit<scalar>::mm};
  for (const auto &slab :
       d.material_store().template get<material_id::e_material_slab>()) {
    EXPECT_EQ(slab.get_material(), silicon<scalar>());
    EXPECT_NEAR(slab.thickness(), thickness, tol);
    thickness += 1.f * unit<scalar>::mm;
  }
}

/// Build a detector with identical material on surfaces in different volumes
/// with and without material deduplication
GTEST_TEST(detray_builders, homogeneous_material_deduplication) {
  vecmem::host_memory_resource host_mr;
  volume_builder_options builder_opts{};

  builder_opts.deduplicate(false);
  const detector_t ref_det = build_detector(host_mr, builder_opts);
  builder_opts.deduplicate(true);
  const detector_t det = build_detector(host_mr, builder_opts);

  EXPECT_TRUE(detail::check_consistency(ref_det));
  EXPECT_TRUE(detail::check_consistency(det));

  ASSERT_EQ(ref_det.volumes().size(), n_volumes);
  ASSERT_EQ(det.volumes().size(), n_volumes);
  ASSERT_EQ(ref_det.surfaces().size(), det.surfaces().size());

  // Without deduplication: One material entry per surface
  const auto &ref_store = ref_det.material_store();
  EXPECT_EQ(ref_store.template size<mat_id::e_material_slab>(), 2u * n_volumes);

  // With deduplication: The shared material is only present once
  const auto &store = det.material_store();
  EXPECT_EQ(store.template size<mat_id::e_material_slab>(), 1u + n_volumes);

  // Check the material links
  std::set<dindex> ref_slab_links;
  std::set<dindex> unique_slab_links;
  for (dindex v = 0u; v < n_volumes; ++v) {
    const auto vol_desc = det.volume(v);
    ASSERT_EQ(vol_desc, ref_det.volume(v));

    // Every surface has its own material entry in the ref detector
    for (dindex sf_idx : {shared_slab_sf, unique_slab_sf}) {
      const auto sf_desc = ref_det.surface(vol_desc.to_global_sf_index(sf_idx));
      const auto link = sf_desc.material();
      EXPECT_EQ(link.id(), mat_id::e_material_slab);
      // Insertion is successful if the index is uniuqe
      EXPECT_TRUE(ref_slab_links.insert(link.index()).second);
    }

    // The shared material is the first entry that was added
    auto sf_desc = det.surface(vol_desc.to_global_sf_index(shared_slab_sf));
    auto link = sf_desc.material();
    EXPECT_EQ(link.id(), mat_id::e_material_slab);
    EXPECT_EQ(link.index(), 0u);

    sf_desc = det.surface(vol_desc.to_global_sf_index(unique_slab_sf));
    link = sf_desc.material();
    EXPECT_EQ(link.id(), mat_id::e_material_slab);
    EXPECT_TRUE(unique_slab_links.insert(link.index()).second);
  }

  EXPECT_FALSE(unique_slab_links.contains(0u));

  // The material lookups are unchanged on every surface
  const std::array<point2, 5> loc_points{point2{0.f, 0.f}, point2{-9.f, -9.f},
                                         point2{3.5f, -2.f}, point2{9.f, 9.f},
                                         point2{-5.5f, 7.f}};

  std::size_t n_checked{0u};
  for (const auto ref_sf_desc : ref_det.surfaces()) {
    const geometry::surface ref_sf{ref_det, ref_sf_desc};
    const auto sf_desc = det.surface(ref_sf.index());
    const geometry::surface sf{det, sf_desc};

    // Check material link
    ASSERT_EQ(ref_sf.has_material(), sf.has_material());
    if (!sf.has_material()) {
      continue;
    }
    ASSERT_EQ(ref_sf_desc.material().id(), sf_desc.material().id());

    // Check lookup
    for (const auto &p : loc_points) {
      const auto ref_slab =
          ref_sf.template visit_material<get_material_slab>(p);
      const auto slab = sf.template visit_material<get_material_slab>(p);

      EXPECT_EQ(ref_slab, slab) << "surface " << ref_sf.index() << ", point ("
                                << p[0] << ", " << p[1] << ")";
      ++n_checked;
    }
  }
  // Have 2 material surfaces per volume
  EXPECT_EQ(n_checked, 2u * loc_points.size() * n_volumes);
}
