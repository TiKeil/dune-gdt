// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <algorithm>
#include <memory>
#include <tuple>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/refinement.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/spaces/fv/default.hh>


template <class GridViewType, size_t r>
struct FiniteVolumeSpace : public ::testing::Test
{
  using R = double;
  using SpaceType = Dune::GDT::FvSpace<GridViewType, R, r>;
  using D = typename SpaceType::D;
  static const constexpr size_t d = SpaceType::d;

  virtual std::shared_ptr<GridViewType> grid_view() = 0;

  std::shared_ptr<SpaceType> space;

  ~FiniteVolumeSpace() = default;

  void SetUp() override final
  {
    ASSERT_NE(grid_view(), nullptr);
    space = std::shared_ptr<SpaceType>(new SpaceType(*grid_view()));
  }

  void TearDown() override final
  {
    space.reset();
  }

  void basis_exists_on_each_element_with_correct_size()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view()))
      EXPECT_EQ(1, space->basis().localize(element)->size());
  }

  void basis_exists_on_each_element_with_correct_order()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view()))
      EXPECT_EQ(0, space->basis().localize(element)->order());
  }

  void mapper_reports_correct_num_DoFs_on_each_element()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    for (auto&& element : elements(*grid_view()))
      EXPECT_EQ(r, space->mapper().local_size(element));
  }

  void mapper_reports_correct_max_num_DoFs()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    size_t max_num_dofs = 0;
    for (auto&& element : elements(*grid_view()))
      max_num_dofs = std::max(max_num_dofs, space->mapper().local_size(element));
    EXPECT_LE(max_num_dofs, space->mapper().max_local_size());
  }

  void mapper_maps_correctly()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    // we want to check that the numbering is consecutive and that each global index exists only once
    std::set<size_t> global_indices;
    // we test both call variants
    std::set<size_t> map_to_global;
    for (auto&& element : Dune::elements(*grid_view())) {
      for (auto&& global_index : space->mapper().global_indices(element))
        global_indices.insert(global_index);
      for (size_t ii = 0; ii < space->mapper().local_size(element); ++ii)
        map_to_global.insert(space->mapper().global_index(element, ii));
    }
    // check for consecutive numbering
    EXPECT_EQ(0, *global_indices.begin());
    EXPECT_EQ(global_indices.size() - 1, *global_indices.rbegin());
    EXPECT_EQ(0, *map_to_global.begin());
    EXPECT_EQ(map_to_global.size() - 1, *map_to_global.rbegin());
    // check that the mapper is of the same opinion
    EXPECT_EQ(space->mapper().size(), global_indices.size());
    EXPECT_EQ(space->mapper().size(), map_to_global.size());
    // check that each global index is unique
    for (const auto& global_index : global_indices)
      EXPECT_EQ(1, global_indices.count(global_index));
    for (const auto& global_index : map_to_global)
      EXPECT_EQ(1, map_to_global.count(global_index));
  } // ... mapper_maps_correctly(...)

  void basis_is_finite_volume_basis()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    double tolerance = 1e-15;
    for (auto&& element : elements(*grid_view())) {
      const auto values = space->basis().localize(element)->evaluate(Dune::FieldVector<D, d>(0.));
      EXPECT_EQ(1, values.size());
      ASSERT_TRUE(Dune::XT::Common::FloatCmp::eq(values.at(0), Dune::FieldVector<R, r>(1), tolerance, tolerance));
    }
  } // ... basis_is_finite_volume_basis(...)

  void basis_jacobians_seem_to_be_correct()
  {
    ASSERT_NE(grid_view(), nullptr);
    ASSERT_NE(space, nullptr);
    double tolerance = 1e-15;
    for (auto&& element : elements(*grid_view())) {
      const auto grads = space->basis().localize(element)->jacobian(Dune::FieldVector<D, d>(0.));
      EXPECT_EQ(1, grads.size());
      ASSERT_TRUE(Dune::XT::Common::FloatCmp::eq(grads.at(0), decltype(grads[0])(0), tolerance, tolerance));
    }
  } // ... basis_jacobians_seem_to_be_correct(...)
}; // struct FiniteVolumeSpace


template <class G, size_t r>
struct FiniteVolumeSpaceOnSimplicialLeafView
    : public FiniteVolumeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, r>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  GridProviderType grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  FiniteVolumeSpaceOnSimplicialLeafView() //     (i) negative coordinates and not the same as the reference
      : grid_provider(Dune::XT::Grid::make_cube_grid<G>(-1.5, -1, 3).grid_ptr()) //                          element,
  { //                                                   (ii) at least 3 elements to have fully inner ones,
    grid_provider.global_refine(1); //                  (iii) refine at least once to obtain all kinds of orientations
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider.leaf_view());
  }

  ~FiniteVolumeSpaceOnSimplicialLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct FiniteVolumeSpaceOnSimplicialLeafView


using SimplicialGrids = ::testing::Types<ONED_1D,
                                         YASP_1D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                         ,
                                         ALU_2D_SIMPLEX_CONFORMING,
                                         ALU_2D_SIMPLEX_NONCONFORMING
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                         ,
                                         UG_2D
#endif
#if HAVE_DUNE_ALUGRID
                                         ,
                                         ALU_3D_SIMPLEX_CONFORMING,
                                         ALU_3D_SIMPLEX_NONCONFORMING
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                         ,
                                         UG_3D
#endif
                                         >;


template <class G>
using ScalarSimplicialFiniteVolumeSpace = FiniteVolumeSpaceOnSimplicialLeafView<G, 1>;
TYPED_TEST_CASE(ScalarSimplicialFiniteVolumeSpace, SimplicialGrids);
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, basis_is_finite_volume_basis)
{
  this->basis_is_finite_volume_basis();
}
TYPED_TEST(ScalarSimplicialFiniteVolumeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, size_t r>
struct FiniteVolumeSpaceOnCubicLeafView
    : public FiniteVolumeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, r>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  FiniteVolumeSpaceOnCubicLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    Dune::FieldVector<D, d> lower_left(-1.5); //  (i) negative coordinates and not the same as the reference element
    Dune::FieldVector<D, d> upper_right(-1.);
    std::array<unsigned int, d> num_elements; // (ii) at least 3 elements to have fully inner ones
    std::fill(num_elements.begin(), num_elements.end(), 3);
    grid_provider = std::make_shared<GridProviderType>(
        Dune::StructuredGridFactory<G>::createCubeGrid(lower_left, upper_right, num_elements));
    leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
  }

  ~FiniteVolumeSpaceOnCubicLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct FiniteVolumeSpaceOnCubicLeafView


using CubicGrids = ::testing::Types<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                    ,
                                    ALU_2D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                    ,
                                    UG_2D
#endif
                                    ,
                                    YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                    ,
                                    ALU_3D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                    ,
                                    UG_3D
#endif
                                    >;


template <class G>
using ScalarCubicFiniteVolumeSpace = FiniteVolumeSpaceOnCubicLeafView<G, 1>;
TYPED_TEST_CASE(ScalarCubicFiniteVolumeSpace, CubicGrids);
TYPED_TEST(ScalarCubicFiniteVolumeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(ScalarCubicFiniteVolumeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(ScalarCubicFiniteVolumeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(ScalarCubicFiniteVolumeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(ScalarCubicFiniteVolumeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(ScalarCubicFiniteVolumeSpace, basis_is_finite_volume_basis)
{
  this->basis_is_finite_volume_basis();
}
TYPED_TEST(ScalarCubicFiniteVolumeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, size_t r>
struct FiniteVolumeSpaceOnPrismLeafView
    : public FiniteVolumeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, r>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  FiniteVolumeSpaceOnPrismLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    if (d == 3) {
      Dune::GridFactory<G> factory;
      for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                            Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                            Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                            Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                            Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                            Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.})}) {
        factory.insertVertex(vertex);
      }
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism, 3), {0, 1, 2, 3, 4, 5});
      grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
      grid_provider->global_refine(1);
      leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
    } else {
      // cannot use ASSERT_... in a ctor
      EXPECT_TRUE(false) << "Does not make sense in " << d << "d!\n"
                         << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
      grid_provider = nullptr;
      leaf_view = nullptr;
    }
  } // FiniteVolumeSpaceOnPrismLeafView(...)

  ~FiniteVolumeSpaceOnPrismLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct FiniteVolumeSpaceOnPrismLeafView


using PrismGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_3D
#endif
    >;


template <class G>
using ScalarPrismFiniteVolumeSpace = FiniteVolumeSpaceOnPrismLeafView<G, 1>;
TYPED_TEST_CASE(ScalarPrismFiniteVolumeSpace, PrismGrids);
TYPED_TEST(ScalarPrismFiniteVolumeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(ScalarPrismFiniteVolumeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(ScalarPrismFiniteVolumeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(ScalarPrismFiniteVolumeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(ScalarPrismFiniteVolumeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(ScalarPrismFiniteVolumeSpace, basis_is_finite_volume_basis)
{
  this->basis_is_finite_volume_basis();
}
TYPED_TEST(ScalarPrismFiniteVolumeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, size_t r>
struct FiniteVolumeSpaceOnPyramidLeafView
    : public FiniteVolumeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, r>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  FiniteVolumeSpaceOnPyramidLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    if (d == 3) {
      Dune::GridFactory<G> factory;
      for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({0, 0, 0}),
                            Dune::XT::Common::FieldVector<D, d>({1, 0, 0}),
                            Dune::XT::Common::FieldVector<D, d>({0, 1, 0}),
                            Dune::XT::Common::FieldVector<D, d>({1, 1, 0}),
                            Dune::XT::Common::FieldVector<D, d>({0, 0, 1})}) {
        factory.insertVertex(vertex);
      }
      factory.insertElement(Dune::GeometryType(Dune::GeometryType::pyramid, 3), {0, 1, 2, 3, 4});
      grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
      grid_provider->global_refine(1);
      leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
    } else {
      // cannot use ASSERT_... in a ctor
      EXPECT_TRUE(false) << "Does not make sense in " << d << "d!\n"
                         << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
      grid_provider = nullptr;
      leaf_view = nullptr;
    }
  } // FiniteVolumeSpaceOnPyramidLeafView(...)

  ~FiniteVolumeSpaceOnPyramidLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct FiniteVolumeSpaceOnPyramidLeafView


using PyramidGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_3D
#endif
    >;


template <class G>
using ScalarPyramidFiniteVolumeSpace = FiniteVolumeSpaceOnPyramidLeafView<G, 1>;
TYPED_TEST_CASE(ScalarPyramidFiniteVolumeSpace, PyramidGrids);
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, basis_is_finite_volume_basis)
{
  this->basis_is_finite_volume_basis();
}
TYPED_TEST(ScalarPyramidFiniteVolumeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}


template <class G, size_t r>
struct FiniteVolumeSpaceOnMixedLeafView
    : public FiniteVolumeSpace<typename Dune::XT::Grid::GridProvider<G>::LeafGridViewType, r>
{
  using GridProviderType = Dune::XT::Grid::GridProvider<G>;
  using LeafGridViewType = typename GridProviderType::LeafGridViewType;

  std::shared_ptr<GridProviderType> grid_provider;
  std::shared_ptr<LeafGridViewType> leaf_view;

  FiniteVolumeSpaceOnMixedLeafView()
  {
    using D = typename G::ctype;
    static const constexpr size_t d = G::dimension;
    switch (d) {
      case 1: {
        // cannot use ASSERT_... in a ctor
        EXPECT_TRUE(false) << "Does not make sense in 1d (all cubes are simplices)!\n"
                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
        grid_provider = nullptr;
        leaf_view = nullptr;
        break;
      }
      case 2: {
        Dune::GridFactory<G> factory;
        for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.75, -1.25})}) {
          factory.insertVertex(vertex);
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), {3, 0, 4, 1});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), {4, 1, 5, 2});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), {4, 6, 3});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), {4, 5, 6});
        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
        grid_provider->global_refine(1);
        leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
        break;
      }
      case 3: {
        Dune::GridFactory<G> factory;
        for (auto&& vertex : {Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1., -1.}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.5, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1.25, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1., -1., -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.5, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1.25, -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.5, -1., -1.5}),
                              Dune::XT::Common::FieldVector<D, d>({-1.75, -1.25, -1.})}) {
          factory.insertVertex(vertex);
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), {3, 0, 4, 1, 9, 6, 10, 7});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), {4, 1, 5, 2, 10, 7, 11, 8});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), {4, 12, 3, 10});
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), {4, 5, 12, 10});
        grid_provider = std::make_shared<GridProviderType>(factory.createGrid());
        grid_provider->global_refine(1);
        leaf_view = std::make_shared<LeafGridViewType>(grid_provider->leaf_view());
        break;
      }
      default: {
        // cannot use ASSERT_... in a ctor
        EXPECT_TRUE(false) << "Not implemented yet for dimension " << d << "!\n"
                           << "=> ALL OTHER TESTS WILL FAIL FOR THIS GRID!!!";
        grid_provider = nullptr;
        leaf_view = nullptr;
      }
    }
  } // FiniteVolumeSpaceOnMixedLeafView(...)

  ~FiniteVolumeSpaceOnMixedLeafView() = default;

  std::shared_ptr<LeafGridViewType> grid_view() override final
  {
    return leaf_view;
  }
}; // struct FiniteVolumeSpaceOnMixedLeafView


using MixedGrids = ::testing::Types<
#if HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D,
    UG_3D
#endif
    >;


template <class G>
using ScalarMixedFiniteVolumeSpace = FiniteVolumeSpaceOnMixedLeafView<G, 1>;
TYPED_TEST_CASE(ScalarMixedFiniteVolumeSpace, MixedGrids);
TYPED_TEST(ScalarMixedFiniteVolumeSpace, basis_exists_on_each_element_with_correct_size)
{
  this->basis_exists_on_each_element_with_correct_size();
}
TYPED_TEST(ScalarMixedFiniteVolumeSpace, basis_exists_on_each_element_with_correct_order)
{
  this->basis_exists_on_each_element_with_correct_order();
}
TYPED_TEST(ScalarMixedFiniteVolumeSpace, mapper_reports_correct_num_DoFs_on_each_element)
{
  this->mapper_reports_correct_num_DoFs_on_each_element();
}
TYPED_TEST(ScalarMixedFiniteVolumeSpace, mapper_reports_correct_max_num_DoFs)
{
  this->mapper_reports_correct_max_num_DoFs();
}
TYPED_TEST(ScalarMixedFiniteVolumeSpace, mapper_maps_correctly)
{
  this->mapper_maps_correctly();
}
TYPED_TEST(ScalarMixedFiniteVolumeSpace, basis_is_finite_volume_basis)
{
  this->basis_is_finite_volume_basis();
}
TYPED_TEST(ScalarMixedFiniteVolumeSpace, basis_jacobians_seem_to_be_correct)
{
  this->basis_jacobians_seem_to_be_correct();
}
