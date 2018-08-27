// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_OPERATORS_ADVECTION_DG_HH

#include <cmath>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker/filters.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/identity.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/advection-dg.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>

#include "interfaces.hh"
#include "localizable-operator.hh"

namespace Dune {
namespace GDT {


/**
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 */
template <class AssemblyGridView,
          class SV,
          size_t m = 1,
          class SF = double,
          class SGV = AssemblyGridView,
          class RF = SF,
          class RGV = AssemblyGridView,
          class RV = SV>
class AdvectionDgOperator : public OperatorInterface<SV,
                                                     SGV,
                                                     m,
                                                     1,
                                                     SF,
                                                     typename XT::Common::multiplication_promotion<SF, RF>::type,
                                                     m,
                                                     1,
                                                     RF,
                                                     RGV,
                                                     RV>
{
  // No need to check the rest, is done in OperatorInterface.
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert((AssemblyGridView::dimension == SGV::dimension) && (AssemblyGridView::dimension == RGV::dimension), "");

  using ThisType = AdvectionDgOperator<AssemblyGridView, SV, m, SF, SGV, RF, RGV, RV>;
  using BaseType = OperatorInterface<SV,
                                     SGV,
                                     m,
                                     1,
                                     SF,
                                     typename XT::Common::multiplication_promotion<SF, RF>::type,
                                     m,
                                     1,
                                     RF,
                                     RGV,
                                     RV>;

protected:
  using AGV = AssemblyGridView;
  static const constexpr size_t d = AssemblyGridView::dimension;
  using D = typename AssemblyGridView::ctype;

public:
  using AssemblyGridViewType = AssemblyGridView;
  using NumericalFluxType = NumericalFluxInterface<d, m, RF>;

  using I = XT::Grid::extract_intersection_t<AGV>;
  using BoundaryTreatmentByCustomNumericalFluxOperatorType =
      LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator<I, SV, SGV, m, SF, RF, RGV, RV>;
  using BoundaryTreatmentByCustomExtrapolationOperatorType =
      LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator<I, SV, SGV, m, SF, RF, RGV, RV>;

  using typename BaseType::SourceSpaceType;
  using typename BaseType::SourceVectorType;

  using typename BaseType::RangeSpaceType;
  using typename BaseType::RangeVectorType;
  using typename BaseType::RangeFunctionType;

  AdvectionDgOperator(
      const AssemblyGridViewType& assembly_grid_view,
      const NumericalFluxType& numerical_flux,
      const SourceSpaceType& source_space,
      const RangeSpaceType& range_space,
      const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>(),
      const size_t jump_indicator_component = 0)
    : BaseType(numerical_flux.parameter_type())
    , assembly_grid_view_(assembly_grid_view)
    , numerical_flux_(numerical_flux.copy())
    , source_space_(source_space)
    , range_space_(range_space)
    , periodicity_exception_(periodicity_exception.copy())
  {
  }

  AdvectionDgOperator(ThisType&& source) = default;

  bool linear() const override final
  {
    return numerical_flux_->linear();
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  /// \name Non-periodic boundary treatment
  /// \{

  ThisType&
  append(typename BoundaryTreatmentByCustomNumericalFluxOperatorType::LambdaType numerical_boundary_treatment_flux,
         const int numerical_boundary_treatment_flux_order,
         const XT::Common::ParameterType& boundary_treatment_parameter_type = {},
         const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<AGV>())
  {
    boundary_treatments_by_custom_numerical_flux_.emplace_back(
        new BoundaryTreatmentByCustomNumericalFluxOperatorType(numerical_boundary_treatment_flux,
                                                               numerical_boundary_treatment_flux_order,
                                                               boundary_treatment_parameter_type),
        filter.copy());
    return *this;
  } // ... append(...)

  ThisType& append(typename BoundaryTreatmentByCustomExtrapolationOperatorType::LambdaType extrapolation,
                   const XT::Common::ParameterType& extrapolation_parameter_type = {},
                   const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<AGV>())
  {
    boundary_treatments_by_custom_extrapolation_.emplace_back(
        new BoundaryTreatmentByCustomExtrapolationOperatorType(
            *numerical_flux_, extrapolation, extrapolation_parameter_type),
        filter.copy());
    return *this;
  }

  /// \}

protected:
  void append_standard_contributions(LocalizableOperatorBase<AGV, SV, m, 1, SF, SGV, m, 1, RF, RGV, RV>& localizable_op,
                                     const XT::Common::Parameter& param) const
  {
    // element contributions
    localizable_op.append(LocalAdvectionDgVolumeOperator<SV, SGV, m, SF, RF, RGV, RV>(numerical_flux_->flux()), param);
    // contributions from inner intersections
    localizable_op.append(LocalAdvectionDgCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>(*numerical_flux_),
                          param,
                          XT::Grid::ApplyOn::InnerIntersectionsOnce<AGV>());
    // contributions from periodic boundaries
    localizable_op.append(LocalAdvectionDgCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>(*numerical_flux_),
                          param,
                          *(XT::Grid::ApplyOn::PeriodicBoundaryIntersectionsOnce<AGV>() && !(*periodicity_exception_)));
    // contributions from other boundaries by custom numerical flux
    for (const auto& boundary_treatment : boundary_treatments_by_custom_numerical_flux_) {
      const auto& boundary_op = *boundary_treatment.first;
      const auto& filter = *boundary_treatment.second;
      localizable_op.append(boundary_op, param, filter);
    }
    // contributions from other boundaries by custom extrapolation
    for (const auto& boundary_treatment : boundary_treatments_by_custom_extrapolation_) {
      const auto& boundary_op = *boundary_treatment.first;
      const auto& filter = *boundary_treatment.second;
      localizable_op.append(boundary_op, param, filter);
    }
  } // ... append_standard_contributions(...)

  void append_local_mass_matrix_inversion(
      LocalizableOperatorBase<AGV, SV, m, 1, SF, SGV, m, 1, RF, RGV, RV>& localizable_op,
      RangeFunctionType& range_function) const
  {
    localizable_op.append([&](const auto& element) {
      // (creating these objects before the grid walk and reusing them would be more efficient, but not thread safe)
      auto local_range = range_function.local_discrete_function(element);
      const auto& range_basis = local_range->basis();
      using E = XT::Grid::extract_entity_t<RGV>;
      const LocalElementIntegralBilinearForm<E, m, 1, RF, RF> local_l2_bilinear_form(
          LocalElementProductIntegrand<E, m, 1, RF, RF>(1.));
      local_range->dofs().assign_from(XT::LA::solve(
          XT::LA::convert_to<XT::LA::CommonDenseMatrix<RF>>(local_l2_bilinear_form.apply2(range_basis, range_basis)),
          XT::LA::convert_to<XT::LA::CommonDenseVector<RF>>(local_range->dofs()),
          {{"type", XT::LA::SolverOptions<XT::LA::CommonDenseMatrix<RF>>::types().at(0)},
           {"post_check_solves_system", "-1"}}));
    });
  } // ... append_local_mass_matrix_inversion(...)

public:
  void
  apply(const SourceVectorType& source, RangeVectorType& range, const XT::Common::Parameter& param = {}) const override
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    range.set_all(0);
    auto source_function = make_discrete_function(this->source_space_, source);
    auto range_function = make_discrete_function(this->range_space_, range);
    // set up the actual operator
    auto localizable_op = make_localizable_operator(this->assembly_grid_view_, source_function, range_function);
    this->append_standard_contributions(localizable_op, param);
    // and apply it in a first grid walk
    localizable_op.assemble(/*use_tbb=*/true);
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
    // prepare inversion of the mass matrix
    this->append_local_mass_matrix_inversion(localizable_op, range_function);
    // and actually do it in a second grid walk
    localizable_op.walk(/*use_tbb=*/true); // Do not call assemble() more than once, will not do anything!
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
  } // ... apply(...)

protected:
  const AssemblyGridViewType& assembly_grid_view_;
  const std::unique_ptr<const NumericalFluxType> numerical_flux_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::unique_ptr<XT::Grid::IntersectionFilter<AssemblyGridViewType>> periodicity_exception_;
  std::list<std::pair<std::unique_ptr<BoundaryTreatmentByCustomNumericalFluxOperatorType>,
                      std::unique_ptr<XT::Grid::IntersectionFilter<AGV>>>>
      boundary_treatments_by_custom_numerical_flux_;
  std::list<std::pair<std::unique_ptr<BoundaryTreatmentByCustomExtrapolationOperatorType>,
                      std::unique_ptr<XT::Grid::IntersectionFilter<AGV>>>>
      boundary_treatments_by_custom_extrapolation_;
}; // class AdvectionDgOperator


template <class VectorType, class AGV, size_t m, class RF, class SGV, class SF, class RGV>
std::enable_if_t<XT::LA::is_vector<VectorType>::value,
                 AdvectionDgOperator<AGV, VectorType, m, SF, SGV, RF, RGV, VectorType>>
make_advection_dg_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<AGV::dimension, m, RF>& numerical_flux,
    const SpaceInterface<SGV, m, 1, SF>& source_space,
    const SpaceInterface<RGV, m, 1, RF>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>())
{
  return AdvectionDgOperator<AGV, VectorType, m, SF, SGV, RF, RGV, VectorType>(
      assembly_grid_view, numerical_flux, source_space, range_space, periodicity_exception);
}


/**
 * \note See AdvectionDgOperator for a description of the template arguments.
 *
 * \sa AdvectionDgOperator
 * \sa [FK2007, Sec. 6.2] for nu_1
 */
template <class AGV,
          class SV,
          size_t m = 1,
          class SF = double,
          class SGV = AGV,
          class RF = SF,
          class RGV = AGV,
          class RV = SV>
class AdvectionDgArtificialViscosityOperator : public AdvectionDgOperator<AGV, SV, m, SF, SGV, RF, RGV, RV>
{
  using ThisType = AdvectionDgArtificialViscosityOperator<AGV, SV, m, SF, SGV, RF, RGV, RV>;
  using BaseType = AdvectionDgOperator<AGV, SV, m, SF, SGV, RF, RGV, RV>;
  using BaseType::d;
  using typename BaseType::D;

public:
  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::NumericalFluxType;

  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;

  using typename BaseType::SourceVectorType;
  using typename BaseType::RangeVectorType;

  AdvectionDgArtificialViscosityOperator(
      const AssemblyGridViewType& assembly_grid_view,
      const NumericalFluxType& numerical_flux,
      const SourceSpaceType& source_space,
      const RangeSpaceType& range_space,
      const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>(),
      const double nu_1 = 0.2,
      const double alpha_1 = 1.,
      const size_t jump_indicator_component = 0)
    : BaseType(assembly_grid_view, numerical_flux, source_space, range_space, periodicity_exception)
    , jump_indicator_component_(jump_indicator_component)
    , nu_1_(nu_1)
    , alpha_1_(alpha_1)
  {
  }

  AdvectionDgArtificialViscosityOperator(ThisType&& source) = default;

  using BaseType::apply;

  void apply(const SourceVectorType& source,
             RangeVectorType& range,
             const XT::Common::Parameter& param = {}) const override final
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    range.set_all(0);
#if 0
    auto src = source.copy();
    auto source_function = make_discrete_function(source_space_, src);
    // apply local P0 projection
    auto walker = XT::Grid::make_walker(assembly_grid_view_);
    walker.append([&](const auto& element) {
      auto local_source = source_function.local_discrete_function(element);
      using E = XT::Grid::extract_entity_t<AGV>;
      const auto average = LocalElementIntegralFunctional<E>(LocalElementIdentityIntegrand<E>()).apply(*local_source)[0]
                           / element.geometry().volume();
      local_source->dofs().assign_from(source_space_.finite_element(element.geometry().type())
                                           .interpolation()
                                           .interpolate([&](const auto& /*xx*/) { return average; }, 0));
    });
    walker.walk(/*use_tbb=*/true);
#endif
    auto source_function = make_discrete_function(this->source_space_, source);
    auto range_function = make_discrete_function(this->range_space_, range);
    // set up the actual operator
    auto localizable_op = make_localizable_operator(this->assembly_grid_view_, source_function, range_function);
    this->append_standard_contributions(localizable_op, param);
    // compute jump indicators for artificial viscosity at detected shocks, see DF2015, p. 449, (8.180)
    // we need a thread safe vector with one entry per grid element and just use a scalar FV discrete function
    const auto fv_space = make_finite_volume_space<1>(this->assembly_grid_view_);
    auto jump_indicators = make_discrete_function<XT::LA::CommonDenseVector<RF>>(fv_space);
    localizable_op.append([&](const auto& element) {
      // compute jump indicator (8.176)
      double element_jump_indicator = 0;
      double element_boundary_without_domain_boundary = (d == 1) ? 1. : 0.;
      const auto local_source_element = source_function.local_discrete_function(element);
      for (auto&& intersection : intersections(this->assembly_grid_view_, element)) {
        if (intersection.neighbor() && !intersection.boundary()) {
          if (d > 1)
            element_boundary_without_domain_boundary += XT::Grid::diameter(intersection);
          const auto neighbor = intersection.outside();
          const auto local_source_neighbor = source_function.local_discrete_function(neighbor);
          const auto integration_order =
              std::pow(std::max(local_source_element->order(), local_source_neighbor->order()), 2);
          for (auto&& quadrature_point :
               QuadratureRules<D, d - 1>::rule(intersection.geometry().type(), integration_order)) {
            const auto point_in_reference_intersection = quadrature_point.position();
            const auto point_in_reference_element =
                intersection.geometryInInside().global(point_in_reference_intersection);
            const auto point_in_reference_neighbor =
                intersection.geometryInOutside().global(point_in_reference_intersection);
            const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
            const auto quadrature_weight = quadrature_point.weight();
            const auto value_on_element =
                local_source_element->evaluate(point_in_reference_element)[jump_indicator_component_];
            const auto value_on_neighbor =
                local_source_neighbor->evaluate(point_in_reference_neighbor)[jump_indicator_component_];
            element_jump_indicator +=
                integration_factor * quadrature_weight * std::pow(value_on_element - value_on_neighbor, 2);
          }
        }
        element_jump_indicator /= element_boundary_without_domain_boundary * element.geometry().volume();
      }
      // compute smoothed discrete jump indicator (8.180)
      double smoothed_discrete_jump_indicator = 0;
      const double xi_min = 0.5;
      const double xi_max = 1.5;
      if (element_jump_indicator < xi_min)
        smoothed_discrete_jump_indicator = 0;
      else if (!(element_jump_indicator < xi_max))
        smoothed_discrete_jump_indicator = 1;
      else
        smoothed_discrete_jump_indicator =
            0.5 * std::sin(M_PI * (element_jump_indicator - (xi_max - xi_min)) / (2 * (xi_max - xi_min))) + 0.5;
      jump_indicators.local_discrete_function(element)->dofs()[0] = smoothed_discrete_jump_indicator;
    });
    // do the actual (first) grid walk: the above operators will be applied and afterwards cleared
    localizable_op.assemble(/*use_tbb=*/true);
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
    // compute artificial viscosity shock capturing [see DF2015, p. 450, (8.183, 8.1.84)] by appending a grid element
    // functor, use the localizable_op as a XT::Grid::Walker
    localizable_op.append([&](const auto& element) {
      // evaluate artificial viscosity form (8.183)
      const auto local_source = source_function.local_discrete_function(element);
      auto local_range = range_function.local_discrete_function(element);
      const auto& local_basis = local_range->basis();
      const auto h = element.geometry().volume();
      const auto jump_indicator_element = jump_indicators.local_discrete_function(element)->dofs()[0];
      for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.type(), 2 * local_basis.order())) {
        const auto point_in_reference_element = quadrature_point.position();
        const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
        const auto quadrature_weight = quadrature_point.weight();
        const auto source_jacobian = local_source->jacobian(point_in_reference_element);
        const auto basis_jacobians = local_basis.jacobians_of_set(point_in_reference_element);
        // compute beta_h
        for (size_t ii = 0; ii < local_basis.size(); ++ii)
          local_range->dofs()[ii] += integration_factor * quadrature_weight * nu_1_ * std::pow(h, alpha_1_)
                                     * jump_indicator_element * (source_jacobian[0] * basis_jacobians[ii][0]);
      }
    });
    // do the actual (second) grid walk
    localizable_op.walk(/*use_tbb=*/true); // Do not call assemble() more than once, will not do anything!
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
    // apply the inverse mass matrix by appending a grid element functor, use the localizable_op as a XT::Grid::Walker
    this->append_local_mass_matrix_inversion(localizable_op, range_function);
    // do the actual (third) grid walk
    localizable_op.walk(/*use_tbb=*/true); // Do not call assemble() more than once, will not do anything!
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
  } // ... apply(...)

private:
  const size_t jump_indicator_component_;
  const double nu_1_;
  const double alpha_1_;
}; // class AdvectionDgArtificialViscosityOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_DG_HH
