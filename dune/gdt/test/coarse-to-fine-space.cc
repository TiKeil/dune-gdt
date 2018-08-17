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

//#include <dune/xt/la/container.hh>
//#include <dune/xt/la/container/eye-matrix.hh>
//#include <dune/xt/la/matrix-inverter.hh>
//#include <dune/xt/la/solver.hh>
//#include <dune/xt/la/container/pattern.hh>
//#include <dune/xt/common/matrix.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>
//#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/dd/glued.hh>
//#include <dune/xt/grid/intersection.hh>
#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/lambda/grid-function.hh>
//#include <dune/xt/functions/constant.hh>
//#include <dune/xt/functions/expression.hh>


//#include <dune/gdt/assembler/global.hh>
//#include <dune/gdt/discretefunction/default.hh>
//#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
//#include <dune/gdt/functionals/l2.hh>
#include <dune/gdt/local/integrands/product.hh>
//#include <dune/gdt/local/integrands/elliptic.hh>
//#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/operators/matrix-based.hh>
//#include <dune/gdt/operators/weighted-l2.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
//#include <dune/gdt/spaces/constraints.hh>

using namespace Dune;
using namespace Dune::GDT;


template <class CoarseGridView, class FineGridView, size_t r = 1, size_t rC = 1, class R = double>
class CoarseSpaceAsFineSpaceWrapper : public SpaceInterface<FineGridView, r, rC, R>
{
  using BaseType = SpaceInterface<FineGridView, r, rC, R>;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::MapperType;
  using typename BaseType::GlobalBasisType;
  using typename BaseType::FiniteElementType;
  using CoarseSpaceType = SpaceInterface<CoarseGridView, r, rC, R>;
  using BaseType::d;

private:
  class CoarseSpaceAsFineSpaceWrapperBasis : public GlobalBasisInterface<FineGridView, r, rC, R>
  {
    using BaseType = GlobalBasisInterface<FineGridView, r, rC, R>;
    using CoarseGlobalBasisType = GlobalBasisInterface<CoarseGridView, r, rC, R>;

  public:
    using typename BaseType::E;
    using typename BaseType::D;
    using BaseType::d;
    using typename BaseType::GridViewType;
    using typename BaseType::ShapeFunctionsType;
    using typename BaseType::LocalizedBasisType;

    CoarseSpaceAsFineSpaceWrapperBasis(const CoarseGlobalBasisType& coarse_basis)
      : coarse_basis_(coarse_basis)
    {
    }

    const GridViewType& grid_view() const
    {
      return coarse_basis_.grid_view();
    }

    const ShapeFunctionsType& shape_functions(const GeometryType& geometry_type) const
    {
      DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
      return coarse_basis_.shape_functions(geometry_type);
    }

    std::unique_ptr<LocalizedBasisType> localize() const
    {
      return std::make_unique<LocalizedCoarseSpaceAsFineSpaceWrapperBasis>(*this);
    }

  private:
    class LocalizedCoarseSpaceAsFineSpaceWrapperBasis : public XT::Functions::ElementFunctionSetInterface<E, r, rC, R>
    {
      using ThisType = LocalizedCoarseSpaceAsFineSpaceWrapperBasis;
      using BaseType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

    public:
      using typename BaseType::ElementType;
      using typename BaseType::DomainType;
      using typename BaseType::RangeType;
      using typename BaseType::DerivativeRangeType;

      LocalizedCoarseSpaceAsFineSpaceWrapperBasis(const CoarseSpaceAsFineSpaceWrapperBasis& self)
        : BaseType()
        , self_(self)
        , coarse_search_(self_.grid_view())
        , current_fine_center_(1)
        , localized_coarse_basis_(self.coarse_basis_.localize())
      {
      }

      LocalizedCoarseSpaceAsFineSpaceWrapperBasis(const ThisType&) = default;
      LocalizedCoarseSpaceAsFineSpaceWrapperBasis(ThisType&&) = default;

      ThisType& operator=(const ThisType&) = delete;
      ThisType& operator=(ThisType&&) = delete;

      size_t max_size(const XT::Common::Parameter& param = {}) const override final
      {
        return localized_coarse_basis_->max_size(param);
      }

    protected:
      void post_bind(const ElementType& elemnt) override final
      {
        current_fine_center_[0] = elemnt.geometry().center();
        auto coarse_search_result = coarse_search_(current_fine_center_);
        DUNE_THROW_IF(coarse_search_result.size() != 1,
                      Exceptions::basis_error,
                      "coarse_search_result.size() = " << coarse_search_result.size());
        const auto& coarse_element_ptr = coarse_search_result[0];
        const auto& coarse_element = *coarse_element_ptr;
        localized_coarse_basis_->bind(coarse_element);
      }

    public:
      size_t size(const XT::Common::Parameter& param = {}) const override final
      {
        return localized_coarse_basis_->size(param);
      }

      int order(const XT::Common::Parameter& param = {}) const override final
      {
        return localized_coarse_basis_->size(param);
      }

      void evaluate(const DomainType& point_in_fine_reference_element,
                    std::vector<RangeType>& result,
                    const XT::Common::Parameter& param = {}) const override final
      {
        DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
        auto point_in_global_coordinates = this->element().geometry().global(point_in_fine_reference_element);
        auto point_in_coarse_reference_element =
            localized_coarse_basis_->element().geometry().local(point_in_global_coordinates);
        localized_coarse_basis_->evaluate(point_in_coarse_reference_element, result, param);
      }

      void jacobians(const DomainType& point_in_fine_reference_element,
                     std::vector<DerivativeRangeType>& result,
                     const XT::Common::Parameter& param = {}) const override final
      {
        DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
        auto point_in_global_coordinates = this->element().geometry().global(point_in_fine_reference_element);
        auto point_in_coarse_reference_element =
            localized_coarse_basis_->element().geometry().local(point_in_global_coordinates);
        localized_coarse_basis_->jacobians(point_in_coarse_reference_element, result, param);
      } // ... jacobian(...)

    private:
      const CoarseSpaceAsFineSpaceWrapperBasis& self_;
      XT::Grid::EntityInlevelSearch<CoarseGridView> coarse_search_;
      std::vector<FieldVector<D, d>> current_fine_center_;
      std::unique_ptr<typename CoarseGlobalBasisType::LocalizedBasisType> localized_coarse_basis_;
    }; // class LocalizedCoarseSpaceAsFineSpaceWrapperBasis

    const CoarseGlobalBasisType& coarse_basis_;
  }; // class CoarseSpaceAsFineSpaceWrapperBasis

  class CoarseSpaceAsFineSpaceWrapperMapper : public MapperInterface<FineGridView>
  {
    using BaseType = MapperInterface<FineGridView>;
    using CoarseMapperType = MapperInterface<CoarseGridView>;
    using MacroElementType = typename CoarseGridView::template Codim<0>::Entity;

  public:
    using typename BaseType::GridViewType;
    using typename BaseType::ElementType;
    using typename BaseType::D;
    using BaseType::d;

    CoarseSpaceAsFineSpaceWrapperMapper(const CoarseMapperType& coarse_mapper, const FineGridView fine_grid_view)
      : wrapped_space_(coarse_mapper)
      , fine_grid_view_(fine_grid_view)
      , max_local_size_(0)
    {
      // This can be done in a single loop. we seperated it for clarity
      // find all coarse elements which are touched by a fine element
      GridViewType macro_grid_view = wrapped_space_.grid_view();
      std::set<typename MacroElementType::EntitySeed,
               std::function<bool(const typename MacroElementType::EntitySeed&,
                                  const typename MacroElementType::EntitySeed)>>
      macro_elements_which_cover_subdomain([&](const auto& lhs, const auto& rhs) {
        const auto ls = macro_grid_view.grid().entity(lhs);
        const auto rs = macro_grid_view.grid().entity(rhs);
        return macro_grid_view.indexSet().index(ls) < macro_grid_view.indexSet().index(rs);
      });
      auto macro_search = XT::Grid::make_entity_in_level_search(macro_grid_view);
      std::vector<FieldVector<double, d>> local_entity_center{FieldVector<double, d>(0.0)};
      for (auto&& fine_element : elements(fine_grid_view_)) {
        local_entity_center[0] = fine_element.geometry().center();
        auto macro_search_result = macro_search(local_entity_center);
        DUNE_THROW_IF(macro_search_result.size() != 1,
                      InvalidStateException,
                      "macro_search_result.size() = " << macro_search_result.size());
        const auto& macro_element_ptr = macro_search_result[0];
        const auto& other_macro_element = *macro_element_ptr;
        macro_elements_which_cover_subdomain.insert(other_macro_element.seed());
      }

      // record all global DoF indices of these elements
      DynamicVector<size_t> global_indices;
      for (const auto& entity_seed : macro_elements_which_cover_subdomain) {
        auto entity = macro_grid_view.grid().entity(entity_seed);
        wrapped_space_.global_indices(entity, global_indices);
        for (size_t ii = 0; ii < wrapped_space_.local_size(entity); ++ii)
          all_global_dofs_.insert(global_indices[ii]);
      }

      // map global DoF to fine elements
      for (auto&& fine_element : elements(fine_grid_view_)) {
        auto fine_id = fine_grid_view_.indexSet().index(fine_element);
        local_entity_center[0] = fine_element.geometry().center();
        auto macro_search_result = macro_search(local_entity_center);
        DUNE_THROW_IF(macro_search_result.size() != 1,
                      InvalidStateException,
                      "macro_search_result.size() = " << macro_search_result.size());
        const auto& macro_element_ptr = macro_search_result[0];
        const auto& macro_element = *macro_element_ptr;
        wrapped_space_.global_indices(macro_element, global_indices);
        DoF_indices_[fine_id] = global_indices;
        max_local_size_ = std::max(max_local_size_, global_indices.size());
      }
    }

    const GridViewType& grid_view() const
    {
      return wrapped_space_.grid_view();
    }

    const LocalFiniteElementCoefficientsInterface<D, d>& local_coefficients(const GeometryType& geometry_type) const
    {
      DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
      return *local_coefficients_;
    }

    size_t size() const
    {
      return all_global_dofs_.size();
    }

    size_t max_local_size() const
    {
      return max_local_size_;
    }

    size_t local_size(const ElementType& element) const
    {
      return DoF_indices_.at(fine_grid_view_.indexSet().index(element)).size();
    }

    size_t global_index(const ElementType& element, const size_t local_index) const
    {
      return DoF_indices_.at(fine_grid_view_.indexSet().index(element))[local_index];
    }

    void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const
    {
      indices = DoF_indices_.at(fine_grid_view_.indexSet().index(element));
    }

  private:
    const CoarseMapperType& wrapped_space_;
    const FineGridView& fine_grid_view_;
    size_t max_local_size_;
    std::map<size_t, DynamicVector<size_t>> DoF_indices_;
    std::set<size_t> all_global_dofs_;

    //    DynamicVector<size_t> all_global_dofs_;
    std::unique_ptr<LocalFiniteElementCoefficientsInterface<D, d>> local_coefficients_; // todo

  }; // class CoarseSpaceAsFineSpaceWrapperMapper

public:
  CoarseSpaceAsFineSpaceWrapper(const CoarseSpaceType& coarse_space, FineGridView fine_grid_view)
    : coarse_space_(coarse_space)
    , fine_grid_view_(fine_grid_view)
    , mapper_(new CoarseSpaceAsFineSpaceWrapperMapper(coarse_space_.mapper(), fine_grid_view_))
    , basis_(new CoarseSpaceAsFineSpaceWrapperBasis(coarse_space_.basis()))
    , fe_(nullptr)
  {
    this->create_communicator();
  }

  const GridViewType& grid_view() const
  {
    return fine_grid_view_;
  }

  const MapperType& mapper() const
  {
    //    DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
    return *mapper_;
  }

  const GlobalBasisType& basis() const
  {
    return *basis_;
  }

  const FiniteElementType& finite_element(const GeometryType& /*geometry_type*/) const
  {
    DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
    return *fe_;
  }

  SpaceType type() const
  {
    return coarse_space_.type();
  }

  int min_polorder() const
  {
    return coarse_space_.min_polorder();
  }

  int max_polorder() const
  {
    return coarse_space_.max_polorder();
  }

  bool continuous(const int diff_order) const
  {
    DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
    return coarse_space_.continuous(diff_order);
  }

  bool continuous_normal_components() const
  {
    DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
    return coarse_space_.continuous_normal_components();
  }

  bool is_lagrangian() const
  {
    DUNE_THROW(Exceptions::space_error, "This does not make sense yet!");
    return coarse_space_.is_lagrangian();
  }

private:
  const CoarseSpaceType& coarse_space_;
  FineGridView fine_grid_view_;
  std::unique_ptr<MapperType> mapper_;
  std::unique_ptr<CoarseSpaceAsFineSpaceWrapperBasis> basis_;
  std::unique_ptr<FiniteElementType> fe_;
}; // class CoarseSpaceAsFineSpaceWrapper


struct LODTest : public ::testing::Test
{
  static const constexpr size_t d = 2;
  using MacroGridType = YASP_2D_EQUIDISTANT_OFFSET;
  using LocalGridType = MacroGridType; // UGGrid<2>;

  using FunctionType = XT::Functions::IndicatorFunction<d>;

  using RangeReturnType = typename FunctionType::RangeReturnType;
  using DomainType = typename FunctionType::DomainType;

  template <class G, bool anything = true>
  struct get_local_layer
  {
    static const constexpr XT::Grid::Layers type = XT::Grid::Layers::level;
  };

  static const constexpr XT::Grid::Layers local_layer = get_local_layer<LocalGridType>::type;

  void SetUp() override final
  {
    if (!macro_grid_)
      macro_grid_ = std::make_unique<XT::Grid::GridProvider<MacroGridType>>(
          XT::Grid::make_cube_grid<MacroGridType>(0., 1., 4, 0));
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    if (!dd_grid_)
      dd_grid_ = std::make_unique<XT::Grid::DD::Glued<MacroGridType, LocalGridType, local_layer>>(
          *macro_grid_,
          2,
          /*prepare_glues=*/false,
          /*allow_for_broken_orientation_of_coupling_intersections=*/true);
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";
    for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
      EXPECT_EQ(dd_grid_->max_local_level(macro_entity), (local_layer == XT::Grid::Layers::level) ? 2 : -1);
    }
  } // ... Setup()

  std::unique_ptr<XT::Grid::GridProvider<MacroGridType>> macro_grid_;
  std::unique_ptr<XT::Grid::DD::Glued<MacroGridType, LocalGridType, local_layer>> dd_grid_;
};


TEST_F(LODTest, coarse_as_fine_space)
{
  /**
   * This is for testing to evaluate a coarse function on a fine mesh
   */

  dd_grid_->setup_oversampling_grids(0, 1);

  auto coarse_view = dd_grid_->macro_grid_view();
  using CGV = decltype(coarse_view);
  using E = XT::Grid::extract_entity_t<CGV>;
  auto coarse_space = make_continuous_lagrange_space<1>(coarse_view);
  auto coarse_local_basis = coarse_space.basis().localize();
  XT::Functions::LambdaGridFunction<E> first_coarse_shape_function(
      /*order_lambda=*/[&](const auto& param) { return coarse_local_basis->order(param); },
      /*post_bind_lambda=*/[&](const auto& e) { coarse_local_basis->bind(e); },
      /*evaluate_lambda=*/[&](const auto& x,
                              const auto& param) { return coarse_local_basis->evaluate_set(x, param)[0]; });
  first_coarse_shape_function.visualize(coarse_view, "coarse_shape_function_0");
  for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
    auto oversampling_view = dd_grid_->local_oversampling_grid(macro_entity).leaf_view();
    using FGV = decltype(oversampling_view);
    using FE = XT::Grid::extract_entity_t<FGV>;
    CoarseSpaceAsFineSpaceWrapper<CGV, FGV, 1> wrapper(coarse_space, oversampling_view);
    auto fine_local_basis = wrapper.basis().localize();
    XT::Functions::LambdaGridFunction<FE> first_fine_shape_function(
        /*order_lambda=*/[&](const auto& param) { return fine_local_basis->order(param); },
        /*post_bind_lambda=*/[&](const auto& e) { fine_local_basis->bind(e); },
        /*evaluate_lambda=*/[&](const auto& x,
                                const auto& param) { return fine_local_basis->evaluate_set(x, param)[0]; });
    first_fine_shape_function.visualize(oversampling_view, "fine_shape_function_0");

    auto fine_space = make_continuous_lagrange_space<1>(oversampling_view);
    auto subdomain_l2_op =
        make_matrix_operator<XT::LA::CommonDenseMatrix<double>>(oversampling_view, wrapper, fine_space);
    subdomain_l2_op.append(LocalElementIntegralBilinearForm<FE>(LocalElementProductIntegrand<FE>(1.)));
    subdomain_l2_op.assemble();


    DUNE_THROW(InvalidStateException, "This is enough");
  }
}
