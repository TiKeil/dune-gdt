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

// /dune-geometry/dune/geometry/referenceelements.hh

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
#include <dune/xt/grid/intersection.hh>
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

template <class MacroElementType>
struct CompareType
{
  // this is for the l2_projection matrix
  bool operator()(const MacroElementType& one, const MacroElementType& other) const
  {
    //    return !(one.equals(other));
    return one.geometry().center() != other.geometry().center();
  }
};


template <class MatrixType, class VectorType>
void SchurComplementSolve(MatrixType& A, MatrixType& C, VectorType& rhs)
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

  std::cout << solution << std::endl;


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

  std::cout << result << std::endl;
  std::cout << result - solution << std::endl;
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

  // this is left to do
  DirichletConstraints<IntersectionType, SpaceType> constraints(*boundary_info, space);

  auto assembler = make_global_assembler(space);
  assembler.append(functional);
  assembler.append(op);
  assembler.append(mass);

  // for now we have implemented the dirichlet constraints manually
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

  auto& mass_matrix = mass.matrix();

  //  std::cout << dirichlet_dofs << std::endl;
  //  logger.info() << "mass matrix = \n" << mass_matrix << "\n\n" << std::endl;
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

  // first we want to know which macro elements are part of the oversampling subdomain
  std::vector<DomainType> centers;
  std::vector<std::set<MacroElementType, CompareType<MacroElementType>>> macro_elements_which_cover_subdomain(
      macro_grid_view.size(0)); // set oder vector?
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    //    std::cout << "macro_element: " << subdomain_id << std::endl;
    for (auto&& fine_element : elements(dd_grid_->local_oversampling_grid(subdomain_id).leaf_view())) {
      for (auto corner_id = 0; corner_id < fine_element.geometry().corners(); ++corner_id) {
        auto corner = fine_element.geometry().corner(corner_id);
        //        std::cout << "corner: " << corner << std::endl;

        for (auto&& other_macro_element : elements(macro_grid_view)) {
          int is_inside = false;
          const auto intersection_it_end = macro_grid_view.iend(other_macro_element);
          for (auto intersection_it = macro_grid_view.ibegin(other_macro_element);
               intersection_it != intersection_it_end;
               ++intersection_it) {
            const auto& intersection = *intersection_it;
            if (XT::Grid::contains(intersection, corner)) {
              is_inside = true;
            }
          }
          if (is_inside) {
            // this is a hack, write a better comparison method! We compare centers in order to detect multiple macro
            // elements
            auto element_center = other_macro_element.geometry().center();
            if (std::find(centers.begin(), centers.end(), element_center) == centers.end()) {
              centers.push_back(element_center);
              macro_elements_which_cover_subdomain[subdomain_id].insert(other_macro_element);
            }
          }
        }
      }
    }
    centers.clear();
    //    std::cout << macro_elements_which_cover_subdomain[subdomain_id].size() << std::endl;
  }
  // validation
  //  for (auto entity = macro_elements_which_cover_subdomain[0].begin();
  //       entity != macro_elements_which_cover_subdomain[0].end();
  //       ++entity) {
  //    std::cout << entity->geometry().center() << std::endl;
  //  }

  // now we want to knwo which dof are in the oversampling subdomain
  std::map<size_t, std::set<size_t>> subdomain_to_DoF_ids;
  DynamicVector<size_t> global_indices;
  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    for (auto entity = macro_elements_which_cover_subdomain[subdomain_id].begin();
         entity != macro_elements_which_cover_subdomain[subdomain_id].end();
         ++entity) {
      const auto other_subdomain_id = dd_grid_->subdomain(*entity);
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
    l2_op.append(LocalElementIntegralBilinearForm<EntityType>(LocalElementProductIntegrand<EntityType>(1.)));
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
        //        std::cout << "lp_in_global_coordinates" << lp_in_global_coordinates << std::endl;
        for (auto entity = macro_elements_which_cover_subdomain[subdomain_id].begin();
             entity != macro_elements_which_cover_subdomain[subdomain_id].end();
             ++entity) {
          //          if (entity->leaf_view().contains(lp_in_global_coordinates))
          //            std::cout << entity->geometry().center() << std::endl;
          const auto lp_in_local_coarse_reference_element = entity->geometry().local(lp_in_global_coordinates);
          local_coarse_basis->bind(*entity);
          if (ReferenceElements<double, 2>::general(entity->geometry().type())
                  .checkInside(lp_in_local_coarse_reference_element)) {
            auto coarse_basis_values = local_coarse_basis->evaluate_set(lp_in_local_coarse_reference_element);
            //            std::cout << coarse_basis_values << std::endl;
            coarse_space.mapper().global_indices(*entity, inner_DoF_ids);
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

    const auto M_l = l2_op.matrix();
    auto C_l = P_l_restricted * M_l;
    //    std::cout << "C_l_ " << C_l << std::endl;
    //    VectorType test(oversampled_fine_space.mapper().size());
    //    test.set_entry(0, 1.);
    //    std::cout << "C_l * test " << C_l * test << std::endl;
    //    DUNE_THROW(InvalidStateException, "this is enough");
  }

  //  using GV = decltype(LocalGrid.leaf_view());
  //  using E = typename GV::template Codim<0>::Entity;

  //  LocalGrid.visualize("LocalGrid");

  //  auto grid = XT::Grid::make_cube_grid<MacroGridType>(0, 1, 4);
  //  grid.visualize("GlobalGrid");

  //  FieldVector<double, 2> lower_left;
  //  FieldVector<double, 2> upper_right;
  //  lower_left[0] = 0.25;
  //  lower_left[1] = 0;
  //  upper_right[0] = 1;
  //  upper_right[1] = 0.75;
  //  std::array<unsigned int, 2> num;
  //  num[0] = 3;
  //  num[1] = 3;
  //  auto wrapped_grid = XT::Grid::make_cube_grid<MacroGridType>(lower_left, upper_right, num);
  //  wrapped_grid.visualize("WrappedGrid");

  //  auto oversampled_fine_space = make_continuous_lagrange_space<1>(LocalGrid.leaf_view());
  //  auto wrapped_coarse_space = make_continuous_lagrange_space<1>(wrapped_grid.leaf_view());
  //  auto op = make_matrix_operator<MatrixType>(LocalGrid.leaf_view(), oversampled_fine_space,
  //  wrapped_coarse_space);

  //  auto grid = XT::Grid::make_cube_grid<MacroGridType>(0, 1, 16);
  //  auto micro_leaf_view = grid.leaf_view();
  //  using GVT = decltype(micro_leaf_view);
  //  using ET = typename GVT::template Codim<0>::Entity;
  //  ContinuousLagrangeSpace<GVT, 1> space(micro_leaf_view);
  //  auto l2_op = make_matrix_operator<MatrixType>(micro_leaf_view, space);
  //  l2_op.append(LocalElementIntegralBilinearForm<ET>(LocalElementProductIntegrand<ET>(1.)));
  //  l2_op.assemble();

  //  const int Nh = 16*16;
  //  const int NH = 16;
  //  Dune::XT::LA::SparsityPatternDefault new_dense_pattern(Nh);
  //  for (size_t ii = 0; ii < Nh; ++ii) {
  //    for (size_t jj = 0; jj < NH; ++jj)
  //      new_dense_pattern.inner(ii).push_back(jj);
  //  }

  //  MatrixType P_h(Nh, NH, new_dense_pattern);
  //  for (size_t ii = 0; ii < Nh; ++ii) {
  //    for (size_t jj = 0; jj < NH; ++jj){
  //      int t = 0;
  //    }
  //  }
}

Test_F(LODTest, saddle_rhs)
{
  /**
   * This is for the rhs of the saddle point problem
   */

  const XT::LA::Backends la = XT::LA::Backends::istl_sparse;
  typedef typename XT::LA::Container<double, la>::MatrixType MatrixType;
  typedef typename XT::LA::Container<double, la>::VectorType VectorType;

  const auto macro_grid_view = dd_grid_->macro_grid_view();
  const auto coarse_space = make_continuous_lagrange_space<1>(macro_grid_view);

  for (auto&& macro_element : elements(macro_grid_view)) {
    const auto subdomain_id = dd_grid_->subdomain(macro_element);
    const auto local_grid = dd_grid_->local_grid(macro_element);
    const auto local_leaf_view = local_grid.leaf_view();

    using GV = decltype(local_leaf_view);
    using E = typename GV::template Codim<0>::Entity;

    const auto local_fine_space = make_continuous_lagrange_space<1>(local_leaf_view);

    // Functional
    //    const LocalElementProductIntegrand<E> product_integrand;
    //    auto integrand = local_binary_to_unary_element_integrand(force.as_grid_function<E>(), product_integrand);
    //    const LocalElementIntegralFunctional<E> integral_functional(integrand);
    //    auto functional = make_vector_functional<VectorType>(local_fine_space);
    //    functional.append(integral_functional);

    //    auto coarse_space = make_continuous_lagrange_space<1>(dd_grid_->macro_grid_view());
    //    auto coarse_basis = coarse_space.basis().localize();
    //    for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
    //      coarse_basis->bind(macro_entity);
    //      auto oversampling_view = dd_grid_->local_grid(macro_entity).leaf_view();
    //      XT::Functions::LambdaFunction<d> global_shape_function(
    //          coarse_basis->order(), [&](const auto& x_in_global_coordinates, const auto& param) {
    //            const auto x_in_macro_reference_coordinates = macro_entity.geometry().local(x_in_global_coordinates);
    //            return coarse_basis->evaluate_set(x_in_macro_reference_coordinates, param)[0];
    //          });
    //      global_shape_function.visualize(oversampling_view, "shape_function_0");
    //    rhs_func.append(LocalElementIntegralFunctional<E>(?????<E>(1.)));
  }
}

// TEST_F(LODTest, SchurComplement)
//{

//}
