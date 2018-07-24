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
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/functions/indicator.hh>
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

using namespace Dune;
using namespace Dune::GDT;


GTEST_TEST(empty, main)
{
  using G = YASP_1D_EQUIDISTANT_OFFSET;
  static const constexpr size_t d = G::dimension;
  auto grid = XT::Grid::make_cube_grid<G>(0, 1, 2);
  auto grid_view = grid.leaf_view();
  using GV = decltype(grid_view);
  using E = typename GV::template Codim<0>::Entity;
  ContinuousLagrangeSpace<GV, 1> space(grid_view);

  const XT::Functions::ConstantFunction<d> func(1);
  const LocalElementProductIntegrand<E> product_integrand(func.as_grid_function<E>());
  const LocalElementIntegralBilinearForm<E> local_op(product_integrand);

  auto op = make_matrix_operator<XT::LA::CommonDenseMatrix<double>>(grid_view, space);
  op.append(local_op);
  op.assemble();

  auto functional =
      make_l2_volume_vector_functional<XT::LA::CommonDenseVector<double>>(space, func.as_grid_function<E>());
  functional.assemble();

  auto assembler = make_global_assembler(space);
  assembler.append(functional);
  assembler.assemble();

  std::cout << "vector = " << functional.vector() << std::endl;
  std::cout << "matrix = \n" << op.matrix() << std::endl;

  XT::LA::CommonDenseVector<double> vector({1, 2, 3});
  auto df = make_discrete_function(space, vector);
  df.visualize("df");
  auto& dofs = df.dofs();
  std::cout << "dofs.vector() = " << dofs.vector() << std::endl;

  auto local_dofs = dofs.localize();
  for (auto&& element : elements(grid_view))
    local_dofs.bind(element);

  std::cout << "functional.apply(df) = " << functional.apply(df) << std::endl;
  std::cout << "op.apply(df).dofs().vector() = " << op.apply(df).dofs().vector() << std::endl;
  std::cout << "op.apply_inverse(df).dofs().vector() = " << op.apply_inverse(df).dofs().vector() << std::endl;
}

/**
 * LOD WORK
 */

using namespace Dune::XT::Grid;

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

  void setup_grids()
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
  } // ... setup()

  std::unique_ptr<GridProvider<MacroGridType>> macro_grid_;
  std::unique_ptr<DD::Glued<MacroGridType, LocalGridType, local_layer>> dd_grid_;
  static const constexpr size_t d = 2;
};

TEST_F(LODTest, standard_problem)
{
  this->setup_grids();
  //  auto micro_leaf_view = dd_grid_->global_grid_view();    // <- this is not working
  auto grid = XT::Grid::make_cube_grid<MacroGridType>(0, 1, 16);
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
  const XT::Functions::ConstantFunction<d, d, d> eye_function(eye);
  const XT::Functions::IndicatorFunction<d> func(init);

  Dune::FieldVector<double, 1> new_value; // RangeType
  new_value[0] = 1;
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

  auto coef = constant + func;
  coef.visualize(micro_leaf_view, "test_indicator");

  const XT::Functions::ConstantFunction<d> force(1);

  const XT::LA::Backends la = XT::LA::Backends::common_dense;
  typedef typename XT::LA::Container<double, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<double, la>::VectorType VectorType;

  auto logger = XT::Common::TimedLogger().get("hi");
  logger.info() << "grid has " << space.grid_view().indexSet().size(0) << " elements" << std::endl;
  typedef typename SpaceType::GridViewType GridViewType;
  typedef XT::Grid::extract_intersection_t<GridViewType> IntersectionType;
  //  auto boundary_info = XT::Grid::BoundaryInfoFactory<IntersectionType>::create(problem.boundary_info_cfg());
  logger.info() << "Assembling... " << std::endl;

  auto op = make_matrix_operator<MatrixType>(micro_leaf_view, space);
  op.append(LocalElementIntegralBilinearForm<E>(
      LocalEllipticIntegrand<E>(force.as_grid_function<E>(), eye_function.as_grid_function<E>())));

  const LocalElementProductIntegrand<E> product_integrand;
  auto integrand = local_binary_to_unary_element_integrand(force.as_grid_function<E>(), product_integrand);
  const LocalElementIntegralFunctional<E> integral_functional(integrand);
  auto functional = make_vector_functional<VectorType>(space);
  functional.append(integral_functional);

  auto assembler = make_global_assembler(space);
  assembler.append(functional);
  assembler.append(op);
  assembler.assemble();

  logger.info() << "...Done " << std::endl;

  auto& system_matrix = op.matrix();
  auto& rhs_vector = functional.vector();

  logger.info() << "system matrix = \n" << system_matrix << "\n\n" << std::endl;
  logger.info() << "rhs vector = \n" << rhs_vector << "\n\n" << std::endl;
  logger.info() << "inverse = \n" << XT::LA::invert_matrix(system_matrix) << "\n\n" << std::endl;

  auto local_solution = XT::LA::solve(system_matrix, rhs_vector);
  make_const_discrete_function(space, local_solution, "local_solution").visualize("local_solution");
}
