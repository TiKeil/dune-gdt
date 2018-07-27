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

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/lambda/function.hh>
#include <dune/xt/functions/constant.hh>

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

using namespace Dune;
using namespace Dune::GDT;


/**
 * LOD WORK
 */

using namespace Dune::XT::Grid;

template <class MatrixType, class VectorType, size_t Nh, size_t NH>
void SchurComplementSolve(MatrixType& A, MatrixType& C, VectorType& rhs)
{
  VectorType full_rhs(Nh + NH);
  for (size_t ii = 0; ii < Nh; ++ii)
    full_rhs.set_entry(ii, rhs.get_entry(ii));

  // solve straight forward

  // this can be done better (low right block is zero)
  Dune::XT::LA::SparsityPatternDefault dense_pattern(Nh + NH);
  for (size_t ii = 0; ii < Nh + NH; ++ii) {
    for (size_t jj = 0; jj < Nh + NH; ++jj)
      dense_pattern.inner(ii).push_back(jj);
  }

  MatrixType full_system(Nh + NH, Nh + NH, dense_pattern);
  for (size_t ii = 0; ii < Nh; ++ii) {
    for (size_t jj = 0; jj < Nh; ++jj) {
      full_system.set_entry(ii, jj, A.get_entry(ii, jj));
    }
    for (size_t kk = 0; kk < NH; ++kk) {
      full_system.set_entry(Nh + kk, ii, C.get_entry(kk, ii));
      full_system.set_entry(ii, Nh + kk, C.get_entry(kk, ii));
    }
  }

  auto full_solution = XT::LA::solve(full_system, full_rhs);
  std::cout << full_solution << std::endl;

  VectorType solution(Nh);
  for (size_t ii = 0; ii < Nh; ++ii)
    solution.set_entry(ii, full_solution.get_entry(ii));

  std::cout << solution << std::endl;


  // pre processing Engwer et al

  Dune::XT::LA::SparsityPatternDefault new_dense_pattern(Nh);
  for (size_t ii = 0; ii < Nh; ++ii) {
    for (size_t jj = 0; jj < NH; ++jj)
      new_dense_pattern.inner(ii).push_back(jj);
  }

  MatrixType Y_pre(Nh, NH, new_dense_pattern);
  for (size_t ii = 0; ii < NH; ++ii) {
    VectorType c_col(Nh);
    for (size_t jj = 0; jj < Nh; ++jj)
      c_col.set_entry(jj, C.get_entry(ii, jj));
    auto y_col = XT::LA::solve(A, c_col);
    for (size_t jj = 0; jj < Nh; ++jj)
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

  std::cout << result << std::endl;
}


struct LODTest : public ::testing::Test
{
  using MacroGridType = YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>;
  using LocalGridType = MacroGridType; // UGGrid<2>;

  using FunctionType = XT::Functions::IndicatorFunction<2>;

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
  static const constexpr size_t d = 2;
};


TEST_F(LODTest, standard_problem)
{
  //  auto micro_leaf_view = dd_grid_->global_grid_view(); // <- this is not working
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
  std::set<size_t> dirichlet_dofs;
  assembler.append([&](const auto& element) {
    std::set<size_t> local_DoFs;
    const auto& fe = space.finite_element(element.type());
    const auto& reference_element = ReferenceElements<double, 2>::general(element.geometry().type());
    const auto local_key_indices = fe.coefficients().local_key_indices();
    const auto intersection_it_end = space.grid_view().iend(element);
    for (auto intersection_it = space.grid_view().ibegin(element); intersection_it != intersection_it_end;
         ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      // actual dirichlet intersections + process boundaries for parallel runs
      if (boundary_info->type(intersection) == XT::Grid::DirichletBoundary()
          || (!intersection.neighbor() && !intersection.boundary())) {
        const auto intersection_index = intersection.indexInInside();
        for (const auto& local_DoF : local_key_indices[1][intersection_index])
          local_DoFs.insert(local_DoF);
        for (int ii = 0; ii < reference_element.size(intersection_index, 1, d); ++ii) {
          const auto element_vertex_id = reference_element.subEntity(intersection_index, 1, ii, d);
          for (const auto& local_DoF : local_key_indices[d][element_vertex_id])
            local_DoFs.insert(local_DoF);
        }
      }
    }
    if (local_DoFs.size() > 0) {
      for (const auto& local_DoF : local_DoFs) {
        dirichlet_dofs.insert(space.mapper().global_index(element, local_DoF));
      }
    }
  });
  assembler.assemble();

  logger.info() << "...Done " << std::endl;

  auto& system_matrix = op.matrix();
  auto& rhs_vector = functional.vector();

  //  std::cout << dirichlet_dofs << std::endl;
  //    logger.info() << "system matrix = \n" << system_matrix << "\n\n" << std::endl;
  //    logger.info() << "rhs vector = \n" << rhs_vector << "\n\n" << std::endl;
  //  // logger.info() << "inverse = \n" << XT::LA::invert_matrix(system_matrix) << "\n\n" << std::endl;

  //  constraints.apply(op.matrix(), functional.vector());
  for (const auto& dof : dirichlet_dofs) {
    system_matrix.unit_row(dof);
    system_matrix.unit_col(dof);
    rhs_vector[dof] = 0.;
  }

  //    logger.info() << "NEW system matrix = \n" << system_matrix << "\n\n" << std::endl;
  //    logger.info() << "NEW rhs vector = \n" << rhs_vector << "\n\n" << std::endl;
  //    logger.info() << "NEW inverse = \n" << XT::LA::invert_matrix(system_matrix) << "\n\n" << std::endl;


  auto local_solution = XT::LA::solve(system_matrix, rhs_vector);
  std::cout << local_solution << std::endl;
  make_const_discrete_function(space, local_solution, "local_solution").visualize("local_solution");

  auto eye_matrix = XT::LA::eye_matrix<MatrixType>(4, 16);

  Dune::XT::LA::SparsityPatternDefault dense_pattern(4);
  for (size_t ii = 0; ii < 4; ++ii) {
    for (size_t jj = 0; jj < 16; ++jj)
      dense_pattern.inner(ii).push_back(jj);
  }

  MatrixType test_matrix(4, 25, dense_pattern);
  for (size_t ii = 0; ii < 4; ++ii)
    test_matrix.set_entry(ii, ii, 1);

  SchurComplementSolve<MatrixType, VectorType, 25, 4>(op.matrix(), test_matrix, rhs_vector);
}

TEST_F(LODTest, DISABLED_local_interpolation)
{
  auto coarse_space = make_continuous_lagrange_space<1>(dd_grid_->macro_grid_view());
  auto coarse_basis = coarse_space.basis().localize();
  for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
    coarse_basis->bind(macro_entity);
    auto local_view = dd_grid_->local_grid(macro_entity).leaf_view();
    XT::Functions::LambdaFunction<d> global_shape_function(
        coarse_basis->order(), [&](const auto& x_in_global_coordinates, const auto& param) {
          const auto x_in_macro_reference_coordinates = macro_entity.geometry().local(x_in_global_coordinates);
          return coarse_basis->evaluate_set(x_in_macro_reference_coordinates, param)[0];
        });
    global_shape_function.visualize(local_view, "shape_function_0");
    DUNE_THROW(InvalidStateException, "this is enough");
  }
}

// TEST_F(LODTest, l2_projection){
//  auto op = make_matrix_operator<MatrixType>(micro_leaf_view, oversampled_fine_space, wrapped_coarse_space);
//  op.append(LocalElementIntegralBilinearForm<E>(LocalElementProductIntegrand<E>(1.)));

//}

// TEST_F(LODTest, SchurComplement)
//{

//}
