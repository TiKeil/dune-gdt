// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016)

/**
  * This file is intended as a starting point for quick testing.
  */

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <functional>

// /dune-geometry/dune/geometry/referenceelements.hh

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/lambda/function.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/expression.hh>


#include <dune/gdt/assembler/global.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/operators/weighted-l2.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/constraints.hh>


template <class T>
struct YaspSeedInequalityComparator
{
  constexpr bool operator()(const T& lhs, const T& rhs) const
  {
    const auto& l = lhs.impl();
    const auto& r = rhs.impl();
    if (l.level() < r.level())
      return true;
    if (l.level() > r.level())
      return false;
    if (l.coord() < r.coord())
      return true;
    if (l.coord() > r.coord())
      return false;
    return l.offset() < r.offset();


    //    return (l.level() == r.level() && l.coord() == r.coord() && );
  }
};


using namespace Dune;
using namespace Dune::GDT;


/**
 * LOD WORK
 */

using namespace Dune::XT::Grid;

template <class DomainType>
struct DomainCompareType
{
  bool operator()(const DomainType& one, const DomainType& other) const
  {
    return XT::Common::FloatCmp::ne(one, other);
  }
};

template <class EntitySeed>
struct CompareType
{
  // this is for the l2_projection matrix
  bool operator()(const EntitySeed& one, const EntitySeed& other) const
  {
    return (one.isValid() == other.isValid());
    //    return XT::Common::FloatCmp::ne(one.geometry().center(), other.geometry().center());
  }
};

template <class MatrixType, class VectorType>
VectorType SchurComplementSolve(MatrixType& A, MatrixType& C, VectorType& rhs)
{
  const size_t rows = A.rows() + C.rows();
  const size_t cols = A.cols() + C.rows();

  VectorType full_rhs(rows);
  for (size_t ii = 0; ii < rhs.size(); ++ii)
    full_rhs.set_entry(ii, rhs.get_entry(ii));

  // solve straight forward

  // this can be done better (low right block is zero)
  Dune::XT::LA::SparsityPatternDefault full_pattern(rows);
  const auto As_pattern = A.pattern();
  for (size_t ii = 0; ii < As_pattern.size(); ++ii) {
    for (const size_t jj : As_pattern.inner(ii)) {
      full_pattern.insert(ii, jj);
    }
  }
  full_pattern.sort();

  const auto Cs_pattern = C.pattern();
  for (size_t ii = 0; ii < Cs_pattern.size(); ++ii) {
    for (const size_t jj : Cs_pattern.inner(ii)) {
      full_pattern.insert(A.rows() + ii, jj);
      full_pattern.insert(jj, A.cols() + ii);
    }
  }
  full_pattern.sort();
  MatrixType full_system(rows, cols, full_pattern);
  for (size_t ii = 0; ii < As_pattern.size(); ++ii) {
    for (const size_t jj : As_pattern.inner(ii)) {
      full_system.set_entry(ii, jj, A.get_entry(ii, jj));
    }
  }
  for (size_t ii = 0; ii < Cs_pattern.size(); ++ii) {
    for (const size_t jj : Cs_pattern.inner(ii)) {
      full_system.set_entry(A.rows() + ii, jj, C.get_entry(ii, jj));
      full_system.set_entry(jj, A.cols() + ii, C.get_entry(ii, jj));
    }
  }

  auto full_solution = XT::LA::solve(full_system, full_rhs);

  VectorType solution(A.cols());
  for (size_t ii = 0; ii < A.cols(); ++ii)
    solution.set_entry(ii, full_solution.get_entry(ii));

  //  std::cout << solution << std::endl;


  // pre processing Engwer et al

  Dune::XT::LA::SparsityPatternDefault new_dense_pattern(A.cols());
  for (size_t ii = 0; ii < A.cols(); ++ii) {
    for (size_t jj = 0; jj < C.rows(); ++jj)
      new_dense_pattern.inner(ii).push_back(jj);
  }

  MatrixType Y_pre(A.cols(), C.rows(), new_dense_pattern);
  for (size_t ii = 0; ii < C.rows(); ++ii) {
    VectorType c_col(A.cols());
    for (size_t jj = 0; jj < A.cols(); ++jj)
      c_col.set_entry(jj, C.get_entry(ii, jj));
    auto y_col = XT::LA::solve(A, c_col);
    for (size_t jj = 0; jj < A.cols(); ++jj)
      Y_pre.set_entry(jj, ii, y_col.get_entry(jj));
  }
  //  std::cout << Y_pre << std::endl;

  auto S = C * Y_pre;
  auto S_inv = XT::LA::invert_matrix(S);

  // post processing Engwer et al
  auto q = XT::LA::solve(A, rhs);

  auto lambda_pre = C * q;
  auto lambda = S_inv * lambda_pre;
  auto result = q - Y_pre * lambda;

  //  std::cout << result << std::endl;
  //  std::cout << result - solution << std::endl;
  return result;
}


struct LODTest : public ::testing::Test
{
  static const constexpr size_t d = 2;
  using MacroGridType = YaspGrid<d, EquidistantOffsetCoordinates<double, d>>;
  using LocalGridType = MacroGridType; // UGGrid<2>;

  using FunctionType = XT::Functions::IndicatorFunction<d>;

  using RangeReturnType = typename FunctionType::RangeReturnType;
  using DomainType = typename FunctionType::DomainType;

  template <class G, bool anything = true>
  struct get_local_layer
  {
    static const constexpr Layers type = Layers::level;
  };

  static const constexpr Layers local_layer = get_local_layer<LocalGridType>::type;

  void SetUp() override final
  {
    if (!macro_grid_)
      macro_grid_ = std::make_unique<GridProvider<MacroGridType>>(make_cube_grid<MacroGridType>(0., 1., 4, 0));
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    if (!dd_grid_)
      dd_grid_ = std::make_unique<DD::Glued<MacroGridType, LocalGridType, local_layer>>(
          *macro_grid_,
          2,
          /*prepare_glues=*/false,
          /*allow_for_broken_orientation_of_coupling_intersections=*/true);
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";
    for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
      EXPECT_EQ(dd_grid_->max_local_level(macro_entity), (local_layer == Layers::level) ? 2 : -1);
    }
  } // ... Setup()

  std::unique_ptr<GridProvider<MacroGridType>> macro_grid_;
  std::unique_ptr<DD::Glued<MacroGridType, LocalGridType, local_layer>> dd_grid_;
};


TEST_F(LODTest, DISABLED_standard_problem_and_saddle_point_problem)
{
  /**
   * This tests a standard FEM with Dirichlet constraints and it also tests the SchurComplement solver
   */

  //  auto micro_leaf_view = dd_grid_->global_grid_view(); // <- this is not working todo
  auto grid = XT::Grid::make_cube_grid<MacroGridType>(0, 1, 4);
  auto micro_leaf_view = grid.leaf_view();
  using GV = decltype(micro_leaf_view);
  using E = typename GV::template Codim<0>::Entity;
  ContinuousLagrangeSpace<GV, 1> space(micro_leaf_view);
  using SpaceType = ContinuousLagrangeSpace<GV, 1>;

  RangeReturnType value(1. - 0.05);
  std::vector<std::pair<XT::Common::FieldMatrix<double, d, 2>, RangeReturnType>> init;
  for (auto xx = 2. / 16.; xx < 1 - 1. / 16.; xx += 4. / 16.) {
    for (auto yy = 2. / 16.; yy < 1 - 1. / 16.; yy += 4. / 16.) {
      std::pair<XT::Common::FieldMatrix<double, d, 2>, RangeReturnType> part;
      part.second = value;
      part.first[0][0] = xx;
      part.first[0][1] = xx + 1. / 16.;
      part.first[1][0] = yy;
      part.first[1][1] = yy + 1. / 16.;
      init.emplace_back(part);
    }
  }

  XT::Common::FieldMatrix<double, d, d> eye(0.);
  for (auto ii = 0; ii < d; ++ii)
    eye[ii][ii] = 1;

  const XT::Functions::ConstantFunction<d> constant(0.05);
  const auto constant_grid_function = constant.as_grid_function<E>();
  const XT::Functions::ConstantFunction<d, d, d> eye_function(eye);
  const XT::Functions::IndicatorFunction<d> func(init);

  Dune::FieldVector<double, 1> new_value; // RangeType
  new_value[0] = 1 - 0.05;
  std::vector<std::pair<XT::Common::FieldMatrix<double, d, 2>, Dune::FieldVector<double, 1>>> new_init;
  for (auto xx = 2. / 16.; xx < 1 - 1. / 16.; xx += 4. / 16.) {
    for (auto yy = 2. / 16.; yy < 1 - 1. / 16.; yy += 4. / 16.) {
      std::pair<XT::Common::FieldMatrix<double, d, 2>, Dune::FieldVector<double, 1>> part;
      part.second = new_value;
      part.first[0][0] = xx;
      part.first[0][1] = xx + 1. / 16.;
      part.first[1][0] = yy;
      part.first[1][1] = yy + 1. / 16.;
      new_init.emplace_back(part);
    }
  }

  const XT::Functions::IndicatorGridFunction<E, 1> funci(new_init);
  funci.visualize(micro_leaf_view, "test_grid_indicator");

  auto coef = constant_grid_function + funci;
  coef.visualize(micro_leaf_view, "test_indicator");

  const XT::Functions::ConstantFunction<d> force(1.);

  const XT::LA::Backends la = XT::LA::Backends::istl_sparse;
  typedef typename XT::LA::Container<double, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<double, la>::VectorType VectorType;

  auto logger = XT::Common::TimedLogger().get("hi");
  logger.info() << "grid has " << space.grid_view().indexSet().size(0) << " elements" << std::endl;
  typedef typename SpaceType::GridViewType GridViewType;
  typedef XT::Grid::extract_intersection_t<GridViewType> IntersectionType;
  auto boundary_info = XT::Grid::BoundaryInfoFactory<IntersectionType>::create();
  logger.info() << "Assembling... " << std::endl;

  // this is a mass matrix (not required for the standard FEM)
  auto mass = make_matrix_operator<MatrixType>(micro_leaf_view, space);
  mass.append(LocalElementIntegralBilinearForm<E>(LocalElementProductIntegrand<E>()));

  auto op = make_matrix_operator<MatrixType>(micro_leaf_view, space);
  op.append(LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(coef, eye_function.as_grid_function<E>())));

  const LocalElementProductIntegrand<E> product_integrand;
  auto integrand = local_binary_to_unary_element_integrand(force.as_grid_function<E>(), product_integrand);
  const LocalElementIntegralFunctional<E> integral_functional(integrand);
  auto functional = make_vector_functional<VectorType>(space);
  functional.append(integral_functional);

  DirichletConstraints<IntersectionType, SpaceType> constraints(*boundary_info, space);

  auto assembler = make_global_assembler(space);
  assembler.append(functional);
  assembler.append(op);
  assembler.append(mass);

  //  // for now we have implemented the dirichlet constraints manually   DONE
  //  std::set<size_t> dirichlet_dofs;
  //  assembler.append([&](const auto& element) {
  //    std::set<size_t> local_DoFs;
  //    const auto& fe = space.finite_element(element.type());
  //    const auto& reference_element = ReferenceElements<double, 2>::general(element.geometry().type());
  //    const auto local_key_indices = fe.coefficients().local_key_indices();
  //    const auto intersection_it_end = space.grid_view().iend(element);
  //    for (auto intersection_it = space.grid_view().ibegin(element); intersection_it != intersection_it_end;
  //         ++intersection_it) {
  //      // only work on dirichlet ones
  //      const auto& intersection = *intersection_it;
  //      // actual dirichlet intersections + process boundaries for parallel runs
  //      if (boundary_info->type(intersection) == XT::Grid::DirichletBoundary()
  //          || (!intersection.neighbor() && !intersection.boundary())) {
  //        const auto intersection_index = intersection.indexInInside();
  //        for (const auto& local_DoF : local_key_indices[1][intersection_index])
  //          local_DoFs.insert(local_DoF);
  //        for (int ii = 0; ii < reference_element.size(intersection_index, 1, d); ++ii) {
  //          const auto element_vertex_id = reference_element.subEntity(intersection_index, 1, ii, d);
  //          for (const auto& local_DoF : local_key_indices[d][element_vertex_id])
  //            local_DoFs.insert(local_DoF);
  //        }
  //      }
  //    }
  //    if (local_DoFs.size() > 0) {
  //      for (const auto& local_DoF : local_DoFs) {
  //        dirichlet_dofs.insert(space.mapper().global_index(element, local_DoF));
  //      }
  //    }
  //  });

  assembler.append(constraints);
  assembler.assemble();

  logger.info() << "...Done " << std::endl;

  auto& system_matrix = op.matrix();
  auto& rhs_vector = functional.vector();

  auto& mass_matrix = mass.matrix();

  //  std::cout << dirichlet_dofs << std::endl;
  logger.info() << "mass matrix = \n" << mass_matrix << "\n\n" << std::endl;
  //    logger.info() << "system matrix = \n" << system_matrix << "\n\n" << std::endl;
  //    logger.info() << "rhs vector = \n" << rhs_vector << "\n\n" << std::endl;
  //  // logger.info() << "inverse = \n" << XT::LA::invert_matrix(system_matrix) << "\n\n" << std::endl;

  constraints.apply(op.matrix(), functional.vector());

  //  for (const auto& dof : dirichlet_dofs) {
  //    system_matrix.unit_row(dof);
  //    system_matrix.unit_col(dof);
  //    rhs_vector[dof] = 0.;
  //  }

  //    logger.info() << "NEW system matrix = \n" << system_matrix << "\n\n" << std::endl;
  //    logger.info() << "NEW rhs vector = \n" << rhs_vector << "\n\n" << std::endl;
  //    logger.info() << "NEW inverse = \n" << XT::LA::invert_matrix(system_matrix) << "\n\n" << std::endl;


  auto local_solution = XT::LA::solve(system_matrix, rhs_vector);
  std::cout << local_solution << std::endl;
  make_const_discrete_function(space, local_solution, "local_solution").visualize("local_solution");

  // visualize if you want

  // SCHUR COMPLEMENT TEST
  Dune::XT::LA::SparsityPatternDefault diagonal_pattern(9);
  for (size_t ii = 0; ii < diagonal_pattern.size(); ++ii)
    for (size_t jj = 0; jj < system_matrix.cols(); ++jj)
      diagonal_pattern.insert(ii, jj);
  diagonal_pattern.sort();

  auto C = XT::LA::eye_matrix<MatrixType>(diagonal_pattern.size(), system_matrix.cols(), diagonal_pattern);
  SchurComplementSolve<MatrixType, VectorType>(system_matrix, C, rhs_vector);
}

TEST_F(LODTest, DISABLED_local_interpolation)
{
  /**
   * This is for testing to evaluate a coarse function on a fine mesh
   */
  auto coarse_space = make_continuous_lagrange_space<1>(dd_grid_->macro_grid_view());
  auto coarse_basis = coarse_space.basis().localize();
  for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
    coarse_basis->bind(macro_entity);
    auto oversampling_view = dd_grid_->local_grid(macro_entity).leaf_view();
    XT::Functions::LambdaFunction<d> global_shape_function(
        coarse_basis->order(), [&](const auto& x_in_global_coordinates, const auto& param) {
          const auto x_in_macro_reference_coordinates = macro_entity.geometry().local(x_in_global_coordinates);
          return coarse_basis->evaluate_set(x_in_macro_reference_coordinates, param)[0];
        });
    global_shape_function.visualize(oversampling_view, "shape_function_0");
    DUNE_THROW(InvalidStateException, "this is enough");
  }
}

TEST_F(LODTest, l2_projection)
{
  /**
   * This is for the projection matrix that belongs to the saddle point problem
   */

  const XT::LA::Backends la = XT::LA::Backends::istl_sparse;
  typedef typename XT::LA::Container<double, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<double, la>::VectorType VectorType;

  const auto macro_grid_view = dd_grid_->macro_grid_view();
  const auto coarse_space = make_continuous_lagrange_space<1>(macro_grid_view);
  dd_grid_->setup_oversampling_grids(0, 1);

  using GV = decltype(macro_grid_view);
  using MacroElementType = typename GV::template Codim<0>::Entity;

  //  std::set<
  //      typename MacroElementType::EntitySeed,
  //      std::function<bool(const typename MacroElementType::EntitySeed&, const typename
  //      MacroElementType::EntitySeed)>>
  //  foo([&](const auto& lhs, const auto& rhs) {
  //    const auto l = macro_grid_view.grid().entity(lhs);
  //    const auto r = macro_grid_view.grid().entity(rhs);
  //    return macro_grid_view.indexSet().index(l) < macro_grid_view.indexSet().index(r);
  //  });


  // first we want to know which macro elements are part of the oversampling subdomain
  //  std::vector<DomainType> centers;
  std::vector<std::set<
      typename MacroElementType::EntitySeed,
      std::function<bool(const typename MacroElementType::EntitySeed&, const typename MacroElementType::EntitySeed)>>>
      macro_elements_which_cover_subdomain;
  for (size_t ii = 0; ii < macro_grid_view.size(0); ++ii)
    macro_elements_which_cover_subdomain.emplace_back([&](const auto& lhs, const auto& rhs) {
      const auto l = macro_grid_view.grid().entity(lhs);
      const auto r = macro_grid_view.grid().entity(rhs);
      return macro_grid_view.indexSet().index(l) < macro_grid_view.indexSet().index(r);
    });
  auto macro_search = make_entity_in_level_search(macro_grid_view);
  std::vector<FieldVector<double, d>> local_entity_center{FieldVector<double, d>(0.0)};
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    //    std::cout << "macro_element: " << subdomain_id << std::endl;
    for (auto&& fine_element : elements(dd_grid_->local_oversampling_grid(subdomain_id).leaf_view())) {
      local_entity_center[0] = fine_element.geometry().center();
      auto macro_search_result = macro_search(local_entity_center);
      DUNE_THROW_IF(macro_search_result.size() != 1,
                    InvalidStateException,
                    "macro_search_result.size() = " << macro_search_result.size());
      const auto& macro_element_ptr = macro_search_result[0];
      const auto& other_macro_element = *macro_element_ptr;
      macro_elements_which_cover_subdomain[subdomain_id].insert(other_macro_element.seed());
    }
  }

  // validation
  //  for (auto entity = macro_elements_which_cover_subdomain[0].begin();
  //       entity != macro_elements_which_cover_subdomain[0].end();
  //       ++entity) {
  //    std::cout << entity->geometry().center() << std::endl;
  //  }

  // now we want to knwo which dofs are in the oversampling subdomain
  std::map<size_t, std::set<size_t>> subdomain_to_DoF_ids;
  DynamicVector<size_t> global_indices;
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    for (const auto& entity_seed : macro_elements_which_cover_subdomain[subdomain_id]) {
      auto entity = macro_grid_view.grid().entity(entity_seed);
      const auto other_subdomain_id = dd_grid_->subdomain(entity);
      coarse_space.mapper().global_indices(macro_element, global_indices);
      for (size_t ii = 0; ii < coarse_space.mapper().local_size(macro_element); ++ii)
        subdomain_to_DoF_ids[other_subdomain_id].insert(global_indices[ii]);
    }
  }

  // validation
  //  for (auto dof = subdomain_to_DoF_ids[0].begin(); dof != subdomain_to_DoF_ids[0].end(); ++dof) {
  //    std::cout << *dof << std::endl;
  //  }

  // now we build P_l
  auto local_coarse_basis = coarse_space.basis().localize();
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    auto oversampling_grid = dd_grid_->local_oversampling_grid(subdomain_id);
    auto oversampling_view = oversampling_grid.leaf_view();
    using OversamplingGV = decltype(oversampling_view);
    using EntityType = typename OversamplingGV::template Codim<0>::Entity;
    const auto oversampled_fine_space = make_continuous_lagrange_space<1>(oversampling_view);
    auto walker = XT::Grid::make_walker(oversampling_view);

    // M_l
    auto l2_op = make_matrix_operator<MatrixType>(oversampled_fine_space);
    l2_op.append(LocalElementIntegralBilinearForm<EntityType>(LocalElementProductIntegrand<EntityType>()));
    walker.append(l2_op);

    // P_l
    XT::LA::SparsityPatternDefault Ps_pattern(
        coarse_space.mapper().size()); // this results to a lot of zeros in the matrix
    for (auto dof = subdomain_to_DoF_ids[subdomain_id].begin(); dof != subdomain_to_DoF_ids[subdomain_id].end();
         ++dof) {
      for (size_t jj = 0; jj < oversampled_fine_space.mapper().size(); ++jj)
        Ps_pattern.insert(*dof, jj);
    }
    Ps_pattern.sort();
    //    std::cout << Ps_pattern.size() << std::endl;
    MatrixType P_l(coarse_space.mapper().size(), oversampled_fine_space.mapper().size(), Ps_pattern);
    DynamicVector<size_t> DoF_ids;
    DynamicVector<size_t> inner_DoF_ids;
    walker.append([&](const auto& element) {
      const auto& lps = oversampled_fine_space.finite_element(element.geometry().type()).lagrange_points();
      oversampled_fine_space.mapper().global_indices(element, DoF_ids);
      for (size_t jj = 0; jj < lps.size(); ++jj) {
        const auto& lp_in_local_fine_reference_element_coordinates = lps[jj];
        const auto lp_in_global_coordinates = element.geometry().global(lp_in_local_fine_reference_element_coordinates);
        for (const auto& entity_seed : macro_elements_which_cover_subdomain[subdomain_id]) {
          auto entity = macro_grid_view.grid().entity(entity_seed);
          const auto lp_in_local_coarse_reference_element = entity.geometry().local(lp_in_global_coordinates);
          local_coarse_basis->bind(entity);
          if (ReferenceElements<double, d>::general(entity.geometry().type())
                  .checkInside(lp_in_local_coarse_reference_element)) {
            auto coarse_basis_values = local_coarse_basis->evaluate_set(lp_in_local_coarse_reference_element);
            //            std::cout << coarse_basis_values << std::endl;
            coarse_space.mapper().global_indices(entity, inner_DoF_ids);
            //            std::cout << "inner dof" << inner_DoF_ids << std::endl;
            for (size_t ii = 0; ii < coarse_basis_values.size(); ++ii) {
              // find index for P_l corresponding to ii
              P_l.set_entry(inner_DoF_ids[ii], DoF_ids[jj], coarse_basis_values[ii]);
            }
          }
        }
      }
    });
    walker.walk();
    //    std::cout << "P_l " <<  P_l << std::endl;

    XT::LA::SparsityPatternDefault dense_pattern(subdomain_to_DoF_ids[subdomain_id].size());
    for (size_t ii = 0; ii < subdomain_to_DoF_ids[subdomain_id].size(); ++ii) {
      for (size_t jj = 0; jj < oversampled_fine_space.mapper().size(); ++jj)
        dense_pattern.insert(ii, jj);
    }
    dense_pattern.sort();

    MatrixType P_l_restricted(
        subdomain_to_DoF_ids[subdomain_id].size(), oversampled_fine_space.mapper().size(), dense_pattern);
    auto iterator = subdomain_to_DoF_ids[subdomain_id].begin();
    for (size_t ii = 0; ii < subdomain_to_DoF_ids[subdomain_id].size(); ++ii) {
      for (size_t jj = 0; jj < oversampled_fine_space.mapper().size(); ++jj) {
        P_l_restricted.set_entry(ii, jj, P_l.get_entry(*(iterator), jj));
      }
      ++iterator;
    }
    //    std::cout << "P_l_restricted " << P_l_restricted << std::endl;

    auto P_l_transposed = XT::Common::transposed<MatrixType>(P_l_restricted);
    auto matrix_result = P_l_transposed * P_l_restricted;
    //    std::cout << "P.T * P : " << matrix_result << std::endl;

    const auto& M_l = l2_op.matrix();
    //    std::cout << "Mass:\n\n" << M_l << std::endl;
    //    std::cout << "M * P.T * P : " << M_l * matrix_result << std::endl;
    auto C_l = P_l_restricted * M_l;


    /**
     * just some stuff
     */

    auto C_l_transposed = XT::Common::transposed<MatrixType>(C_l);

    //    std::cout << "C_l_ " << C_l << std::endl;

    // validation
    VectorType test_vector(oversampled_fine_space.mapper().size(), 0.);
    test_vector[0] = 1.;
    auto project = C_l * test_vector;
    std::cout << "C_l * test " << project << std::endl;
    std::cout << "C_l_transposed times C_l" << C_l_transposed * C_l << std::endl;


    // better validation
    VectorType test_vector_2(oversampled_fine_space.mapper().size(), 0.);
    XT::Functions::LambdaFunction<d> global_shape_function(
        local_coarse_basis->order(), [&](const auto& x_in_global_coordinates, const auto& param = {}) {
          const auto x_in_macro_reference_coordinates = macro_element.geometry().local(x_in_global_coordinates);
          return local_coarse_basis->evaluate_set(x_in_macro_reference_coordinates, param)[0];
        });
    // find global DoFs
    //    DynamicVector<size_t> DoF_ids;
    std::vector<DomainType> global_lps(oversampled_fine_space.mapper().size());
    for (auto&& micro_element : elements(oversampling_view)) {
      const auto& lps = oversampled_fine_space.finite_element(micro_element.geometry().type()).lagrange_points();
      oversampled_fine_space.mapper().global_indices(micro_element, DoF_ids);
      for (size_t jj = 0; jj < lps.size(); ++jj) {
        const auto& lp_in_local_fine_reference_element_coordinates = lps[jj];
        const auto lp_in_global_coordinates =
            micro_element.geometry().global(lp_in_local_fine_reference_element_coordinates);
        if (ReferenceElements<double, d>::general(macro_element.geometry().type())
                .checkInside(macro_element.geometry().local(lp_in_global_coordinates))) {
          global_lps[DoF_ids[jj]] = lp_in_global_coordinates;
          test_vector_2[DoF_ids[jj]] = global_shape_function.evaluate(lp_in_global_coordinates);
        }
      }
    }

    //    std::cout << global_lps << std::endl;

    std::cout << "test_vector_2" << test_vector_2 << std::endl;
    //    auto proj_result = C_l * test_vector_2;
    //    std::cout << "C_l * test_vector_2 " << proj_result << std::endl;

    VectorType back_vector(P_l_transposed.cols(), 0);
    back_vector[0] = 1.;
    auto back_result = P_l_transposed * back_vector;
    std::cout << "P_l_transposed * back_vector " << back_result << std::endl;
    std::cout << "C_l * P_l_transposed" << C_l * P_l_transposed << std::endl;

    VectorType one_vector(P_l_restricted.cols(), 1.);
    auto one_result = P_l_restricted * one_vector;
    std::cout << "P_l_transposed * one_vector " << one_result << std::endl;
    auto proj_ = C_l * one_vector;
    std::cout << "C_l  * one_vector " << proj_ << std::endl;
    std::cout << "back to fine: " << P_l_transposed * proj_ << std::endl;

    DUNE_THROW(InvalidStateException, "this is enough");
  }
}

TEST_F(LODTest, DISABLED_saddle_rhs)
{
  /**
   * This is for the rhs of the saddle point problem
   */

  const XT::LA::Backends la = XT::LA::Backends::istl_sparse;
  typedef typename XT::LA::Container<double, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<double, la>::VectorType VectorType;

  const auto macro_grid_view = dd_grid_->macro_grid_view();
  const auto coarse_space = make_continuous_lagrange_space<1>(macro_grid_view);

  using GVT = decltype(macro_grid_view);
  using SpaceType = ContinuousLagrangeSpace<GridView<GVT>, 1, double>;

  auto coarse_basis = coarse_space.basis().localize();

  // bilinear form kappa
  XT::Functions::ExpressionFunction<d> kappa("x", {"-x[0]"}, 1);
  //  kappa.visualize(macro_grid_view, "hi");

  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    const auto local_grid = dd_grid_->local_grid(macro_element);
    const auto local_leaf_view = local_grid.leaf_view();

    using GV = decltype(local_leaf_view);
    using E = typename GV::template Codim<0>::Entity;

    const auto local_fine_space = make_continuous_lagrange_space<1>(local_leaf_view);

    // rhs
    coarse_basis->bind(macro_element);
    XT::Functions::LambdaFunction<d> global_shape_function(
        coarse_basis->order(),
        [&](const auto& x_in_global_coordinates, const auto& param) {
          const auto x_in_macro_reference_coordinates = macro_element.geometry().local(x_in_global_coordinates);
          return coarse_basis->evaluate_set(x_in_macro_reference_coordinates, param)[0];
        },
        "smooth_lambda_function",
        {},
        [&](const auto& x_in_global_coordinates, const auto& param) {
          const auto x_in_macro_reference_coordinates = macro_element.geometry().local(x_in_global_coordinates);
          using DerivativeRangeType = typename SpaceType::GlobalBasisType::LocalizedBasisType::DerivativeRangeType;
          std::vector<DerivativeRangeType> result;
          coarse_basis->jacobians(x_in_macro_reference_coordinates, result, param);
          return result[0];
        });

    //    global_shape_function.visualize(local_leaf_view, "shape_function_test");

    // ConstantFunction constant_function(1.);
    //    const auto constant_grid_function = constant_function.as_grid_function<E>();

    //    // this is for the bilinear form

    //    Dune::FieldVector<double, 1> new_value; // RangeType
    //    new_value[0] = 1 - 0.05;
    //    std::vector<std::pair<XT::Common::FieldMatrix<double, d, 2>, Dune::FieldVector<double, 1>>> new_init;
    //    for (auto xx = 2. / 16.; xx < 1 - 1. / 16.; xx += 4. / 16.) {
    //      for (auto yy = 2. / 16.; yy < 1 - 1. / 16.; yy += 4. / 16.) {
    //        std::pair<XT::Common::FieldMatrix<double, d, 2>, Dune::FieldVector<double, 1>> part;
    //        part.second = new_value;
    //        part.first[0][0] = xx;
    //        part.first[0][1] = xx + 1. / 16.;
    //        part.first[1][0] = yy;
    //        part.first[1][1] = yy + 1. / 16.;
    //        new_init.emplace_back(part);
    //      }
    //    }

    //    const XT::Functions::IndicatorGridFunction<E, 1> funci(new_init);
    //    auto coef = constant_grid_function + funci;
    const auto coef = kappa.as_grid_function<E>();

    XT::Common::FieldMatrix<double, d, d> eye(0.); // can I do this with XT::LA::eye_matrix?
    for (auto ii = 0; ii < d; ++ii)
      eye[ii][ii] = 1;
    const XT::Functions::ConstantFunction<d, d, d> eye_function(eye);

    const LocalEllipticIntegrand<E> elliptic_integrand(coef, eye_function.as_grid_function<E>());
    auto integrand =
        local_binary_to_unary_element_integrand(global_shape_function.as_grid_function<E>(), elliptic_integrand);
    const LocalElementIntegralFunctional<E> integral_functional(integrand);
    auto functional = make_vector_functional<VectorType>(local_fine_space);
    functional.append(integral_functional);
    auto walker = XT::Grid::make_walker(local_leaf_view);
    walker.append(functional);
    walker.walk();
    std::cout << functional.vector() << std::endl;
  }
}

TEST_F(LODTest, Corrector_problem)
{
  // PREPARATIONS
  const XT::LA::Backends la = XT::LA::Backends::istl_sparse;
  typedef typename XT::LA::Container<double, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<double, la>::VectorType VectorType;

  const auto macro_grid_view = dd_grid_->macro_grid_view();
  const auto coarse_space = make_continuous_lagrange_space<1>(macro_grid_view);
  dd_grid_->setup_oversampling_grids(0, 0);

  using Macro_GV = decltype(macro_grid_view);
  using MacroElementType = typename Macro_GV::template Codim<0>::Entity;


  /**
    * PREPARE our grid (THIS IS NOT SUPPOSED TO BE HERE. PUT THIS INTO GRID_GLUED)
    */
  // first we want to know which macro elements are part of the oversampling subdomain
  std::vector<DomainType> centers;
  std::vector<std::set<typename MacroElementType::EntitySeed, CompareType<typename MacroElementType::EntitySeed>>>
      macro_elements_which_cover_subdomain(macro_grid_view.size(0)); // set oder vector?
  auto macro_search = make_entity_in_level_search(macro_grid_view);
  std::vector<FieldVector<double, d>> local_entity_center{FieldVector<double, d>(0.0)};
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    //    std::cout << "macro_element: " << subdomain_id << std::endl;
    for (auto&& fine_element : elements(dd_grid_->local_oversampling_grid(subdomain_id).leaf_view())) {
      local_entity_center[0] = fine_element.geometry().center();
      auto macro_search_result = macro_search(local_entity_center);
      DUNE_THROW_IF(macro_search_result.size() != 1,
                    InvalidStateException,
                    "macro_search_result.size() = " << macro_search_result.size());
      const auto& macro_element_ptr = macro_search_result[0];
      const auto& other_macro_element = *macro_element_ptr;
      auto element_center = other_macro_element.geometry().center();
      if (std::find(centers.begin(), centers.end(), element_center) == centers.end()) {
        // this is a hack, write a better comparison method! We compare centers in order to detect multiple elements
        centers.push_back(element_center);
        macro_elements_which_cover_subdomain[subdomain_id].insert(other_macro_element.seed());
      }
    }
    centers.clear();
  }

  // now we want to knwo which dofs are in the oversampling subdomain
  std::map<size_t, std::set<size_t>> subdomain_to_DoF_ids;
  DynamicVector<size_t> global_indices;
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    for (const auto& entity_seed : macro_elements_which_cover_subdomain[subdomain_id]) {
      auto entity = macro_grid_view.grid().entity(entity_seed);
      const auto other_subdomain_id = dd_grid_->subdomain(entity);
      coarse_space.mapper().global_indices(macro_element, global_indices);
      for (size_t ii = 0; ii < coarse_space.mapper().local_size(macro_element); ++ii)
        subdomain_to_DoF_ids[other_subdomain_id].insert(global_indices[ii]);
    }
  }

  /**
   * main algorithm ! ******************************************************************************
   */

  // take one particular macro element and compute a corrector!
  auto local_coarse_basis = coarse_space.basis().localize();
  const size_t subdomain_id = 6;
  auto oversampling_grid = dd_grid_->local_oversampling_grid(subdomain_id);
  auto oversampling_view = oversampling_grid.leaf_view();
  using GV = decltype(oversampling_view);
  using E = typename GV::template Codim<0>::Entity;
  using SpaceType = ContinuousLagrangeSpace<GV, 1>;

  const auto oversampled_fine_space = make_continuous_lagrange_space<1>(oversampling_view);
  const auto local_grid = dd_grid_->local_grid(subdomain_id);

  // we need to derive macro_element from the subindex
  local_entity_center[0] = local_grid.leaf_view().begin<0>()->geometry().center();
  auto macro_search_result = macro_search(local_entity_center);
  const auto& macro_element_ptr = macro_search_result[0];
  const auto& macro_element = *macro_element_ptr;

  //  oversampling_grid.visualize("oversampling_grid_id_6");
  //  local_grid.visualize("local_grid_id_6");

  auto oversampling_walker = XT::Grid::make_walker(oversampling_view);

  // bilinear form kappa
  XT::Functions::ConstantFunction<d> kappa(1.);
  //  XT::Functions::ExpressionFunction<d> kappa("x", {"-x[0]"}, 1);
  const auto coef = kappa.as_grid_function<E>();

  /**
   * C_l MATRIX
   */

  // M_l
  auto l2_op = make_matrix_operator<MatrixType>(oversampled_fine_space);
  l2_op.append(LocalElementIntegralBilinearForm<E>(LocalElementProductIntegrand<E>(1.)));
  oversampling_walker.append(l2_op);

  // P_l
  XT::LA::SparsityPatternDefault Ps_pattern(
      coarse_space.mapper().size()); // this results to a lot of zeros in the matrix
  for (auto dof = subdomain_to_DoF_ids[subdomain_id].begin(); dof != subdomain_to_DoF_ids[subdomain_id].end(); ++dof) {
    for (size_t jj = 0; jj < oversampled_fine_space.mapper().size(); ++jj)
      Ps_pattern.insert(*dof, jj);
  }
  Ps_pattern.sort();
  //    std::cout << Ps_pattern.size() << std::endl;
  MatrixType P_l(coarse_space.mapper().size(), oversampled_fine_space.mapper().size(), Ps_pattern);
  DynamicVector<size_t> DoF_ids;
  DynamicVector<size_t> inner_DoF_ids;
  oversampling_walker.append([&](const auto& element) {
    const auto& lps = oversampled_fine_space.finite_element(element.geometry().type()).lagrange_points();
    oversampled_fine_space.mapper().global_indices(element, DoF_ids);
    for (size_t jj = 0; jj < lps.size(); ++jj) {
      const auto& lp_in_local_fine_reference_element_coordinates = lps[jj];
      const auto lp_in_global_coordinates = element.geometry().global(lp_in_local_fine_reference_element_coordinates);
      for (const auto& entity_seed : macro_elements_which_cover_subdomain[subdomain_id]) {
        auto entity = macro_grid_view.grid().entity(entity_seed);
        const auto lp_in_local_coarse_reference_element = entity.geometry().local(lp_in_global_coordinates);
        local_coarse_basis->bind(entity);
        if (ReferenceElements<double, 2>::general(entity.geometry().type())
                .checkInside(lp_in_local_coarse_reference_element)) {
          auto coarse_basis_values = local_coarse_basis->evaluate_set(lp_in_local_coarse_reference_element);
          //            std::cout << coarse_basis_values << std::endl;
          coarse_space.mapper().global_indices(entity, inner_DoF_ids);
          //            std::cout << "inner dof" << inner_DoF_ids << std::endl;
          for (size_t ii = 0; ii < coarse_basis_values.size(); ++ii) {
            // find index for P_l corresponding to ii
            P_l.set_entry(inner_DoF_ids[ii], DoF_ids[jj], coarse_basis_values[ii]);
          }
        }
      }
    }
  });

  /**
   * RIGHT HAND SIDE
   */

  // rhs
  auto coarse_basis = coarse_space.basis().localize();
  coarse_basis->bind(macro_element);
  XT::Functions::LambdaFunction<d> global_shape_function(
      coarse_basis->order(),
      [&](const auto& x_in_global_coordinates, const auto& param) {
        const auto x_in_macro_reference_coordinates = macro_element.geometry().local(x_in_global_coordinates);
        using RangeType = typename SpaceType::GlobalBasisType::LocalizedBasisType::RangeType;
        if (ReferenceElements<double, d>::general(macro_element.geometry().type())
                .checkInside(x_in_macro_reference_coordinates)) {
          return coarse_basis->evaluate_set(x_in_macro_reference_coordinates, param)[0];
        } else {
          return RangeType();
        }

      },
      "smooth_lambda_function",
      {},
      [&](const auto& x_in_global_coordinates, const auto& param) {
        const auto x_in_macro_reference_coordinates = macro_element.geometry().local(x_in_global_coordinates);
        using DerivativeRangeType = typename SpaceType::GlobalBasisType::LocalizedBasisType::DerivativeRangeType;
        std::vector<DerivativeRangeType> result;
        if (ReferenceElements<double, d>::general(macro_element.geometry().type())
                .checkInside(x_in_macro_reference_coordinates)) {
          coarse_basis->jacobians(x_in_macro_reference_coordinates, result, param);
          return result[0];
        } else {
          return DerivativeRangeType();
        }
      });

  XT::Common::FieldMatrix<double, d, d> eye(0.); // can I do this with XT::LA::eye_matrix?
  for (auto ii = 0; ii < d; ++ii)
    eye[ii][ii] = 1;
  const XT::Functions::ConstantFunction<d, d, d> eye_function(eye);

  const LocalEllipticIntegrand<E> elliptic_integrand(coef, eye_function.as_grid_function<E>());
  auto integrand =
      local_binary_to_unary_element_integrand(global_shape_function.as_grid_function<E>(), elliptic_integrand);
  const LocalElementIntegralFunctional<E> integral_functional(integrand);
  auto functional = make_vector_functional<VectorType>(oversampled_fine_space);
  functional.append(integral_functional);
  oversampling_walker.append(functional);

  /**
   * A_l MATRIX
   */

  typedef typename SpaceType::GridViewType GridViewType;
  typedef XT::Grid::extract_intersection_t<GridViewType> IntersectionType;
  auto boundary_info = XT::Grid::BoundaryInfoFactory<IntersectionType>::create();
  auto op = make_matrix_operator<MatrixType>(oversampled_fine_space);
  op.append(LocalElementIntegralBilinearForm<E>(LocalEllipticIntegrand<E>(coef, eye_function.as_grid_function<E>())));
  oversampling_walker.append(op);

  DirichletConstraints<IntersectionType, SpaceType> constraints(*boundary_info, oversampled_fine_space);
  oversampling_walker.append(constraints);

  /**
   * Putting things together
   */

  // waaaaaalk
  oversampling_walker.walk();

  // extract C_l
  XT::LA::SparsityPatternDefault dense_pattern(subdomain_to_DoF_ids[subdomain_id].size());
  for (size_t ii = 0; ii < subdomain_to_DoF_ids[subdomain_id].size(); ++ii) {
    for (size_t jj = 0; jj < oversampled_fine_space.mapper().size(); ++jj)
      dense_pattern.insert(ii, jj);
  }
  dense_pattern.sort();

  MatrixType P_l_restricted(
      subdomain_to_DoF_ids[subdomain_id].size(), oversampled_fine_space.mapper().size(), dense_pattern);
  auto iterator = subdomain_to_DoF_ids[subdomain_id].begin();
  for (size_t ii = 0; ii < subdomain_to_DoF_ids[subdomain_id].size(); ++ii) {
    for (size_t jj = 0; jj < oversampled_fine_space.mapper().size(); ++jj) {
      P_l_restricted.set_entry(ii, jj, P_l.get_entry(*(iterator), jj));
    }
    ++iterator;
  }

  const auto M_l = l2_op.matrix();
  //  std::cout << "\n\nM_l = " << M_l << "\n\n" << std::endl;
  auto C_l = P_l_restricted * M_l;

  std::cout << "foo = "
            << (XT::Common::transposed<MatrixType>(C_l) * C_l - XT::LA::eye_matrix<MatrixType>(C_l.cols(), C_l.cols()))
                   .sup_norm()
            << std::endl;

  auto P_l_transposed = XT::Common::transposed<MatrixType>(P_l_restricted);
  VectorType back_vector(P_l_transposed.cols(), 0);
  back_vector[0] = 1.;
  auto back_result = P_l_transposed * back_vector;
  //  std::cout << "P_l_transposed * back_vector " << back_result << std::endl;
  //  std::cout << "C_l * P_l_transposed" << C_l * P_l_transposed << std::endl;

  // extract A_l restricted
  auto& system_matrix = op.matrix();
  auto& rhs_vector = functional.vector();

  constraints.apply(op.matrix(), functional.vector());

  // apply SCHURCOMPLEMENT solver
  auto corrector_result = SchurComplementSolve<MatrixType, VectorType>(system_matrix, C_l, rhs_vector);
  std::cout << corrector_result.max() << std::endl;
  make_const_discrete_function(oversampled_fine_space, corrector_result, "local_corrector_problem")
      .visualize("local_corrector_problem");
}

TEST_F(LODTest, Global_corrected_matrix)
{
  /**
   * This will add the correctors to the corresponding values in the coarse matrix
   */
}
