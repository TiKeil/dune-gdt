// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
#define DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/algorithms/triangular_solves.hh>
#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/interfaces.hh>

#include "../quadrature.hh"
#include "reconstructed_function.hh"
#include "slopelimiters.hh"

namespace Dune {
namespace GDT {


template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
class LocalReconstructionFvOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t polOrder = 1;
  static constexpr size_t stencil_size = 2 * polOrder + 1;
  static constexpr std::array<size_t, 3> stencil = {
      {2 * polOrder + 1, dimDomain > 1 ? 2 * polOrder + 1 : 1, dimDomain > 2 ? 2 * polOrder + 1 : 1}};
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;
  typedef typename XT::LA::EigenSolver<MatrixType> EigenSolverType;
  typedef typename XT::LA::EigenSolverOptions<MatrixType> EigenSolverOptionsType;
  typedef typename XT::LA::MatrixInverterOptions<MatrixType> MatrixInverterOptionsType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> Quadrature1dType;
  typedef typename GridLayerType::Intersection IntersectionType;
  typedef FieldVector<IntersectionType, 2 * dimDomain> IntersectionVectorType;
  typedef FieldVector<FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>, stencil[0]> ValuesType;
  typedef typename IntersectionType::Geometry::LocalCoordinate IntersectionLocalCoordType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  typedef typename Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain> JacobianRangeType;
  typedef typename XT::Grid::BoundaryInfo<IntersectionType> BoundaryInfoType;
  using ReconstructedFunctionType =
      ReconstructedLocalizableFunction<GridLayerType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>;

public:
  explicit LocalReconstructionFvOperator(const std::vector<RangeType>& source_values,
                                         const AnalyticalFluxType& analytical_flux,
                                         const BoundaryValueType& boundary_values,
                                         const GridLayerType& grid_layer,
                                         const XT::Common::Parameter& param,
                                         const Quadrature1dType& quadrature,
                                         ReconstructedFunctionType& reconstructed_function)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(grid_layer)
    , param_(param)
    , quadrature_(quadrature)
    , reconstructed_function_(reconstructed_function)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    param_.set("boundary", {0.});
    is_instantiated_ = true;
  }

  ~LocalReconstructionFvOperator()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    const Stencil stencil(source_values_, stencil_sizes_, grid_layer_, boundary_values_, entity);
    // In a MPI parallel run, if entity is in overlap, we do not have to reconstruct
    if (!stencil.valid())
      return;
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity))
      intersections[intersection.indexInInside()] = intersection;
    // get jacobians
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_function_.values()[entity_index];
    static thread_local FieldVector<RangeType, dimDomain> tau_;
    if (local_initialization_count_ != initialization_count_) {
      if (!jacobian())
        jacobian() = XT::Common::make_unique<JacobianRangeType>();
      if (!eigenvectors()) {
        eigenvectors() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
        QR() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
      }
      const auto& u_entity = source_values_[entity_index];
      helper<dimDomain>::get_jacobian(analytical_flux_, u_entity, *(jacobian()), param_, entity);
      helper<dimDomain>::get_eigenvectors(*(jacobian()), *(eigenvectors()), *(QR()), tau_, permutations());
      if (analytical_flux_.is_affine()) {
        jacobian() = nullptr;
        ++local_initialization_count_;
      }
    }

    for (size_t dd = 0; dd < dimDomain; ++dd)
      helper<dimDomain>::reconstruct(dd,
                                     values,
                                     *(eigenvectors()),
                                     *(QR()),
                                     tau_,
                                     permutations(),
                                     quadrature_,
                                     reconstructed_values_map,
                                     intersections);
  } // void apply_local(...)

  static void reset()
  {
    ++initialization_count_;
  }

private:
  // quadrature rule containing left and right interface points
  static Quadrature1dType left_right_quadrature()
  {
    Quadrature1dType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper;

  template <class anything>
  struct helper<1, anything>
  {
    static void reconstruct(size_t /*dd*/,
                            const ValuesType& values,
                            FieldVector<MatrixType, dimDomain>& eigenvectors,
                            FieldVector<MatrixType, dimDomain>& QR,
                            FieldVector<RangeType, dimDomain>& tau,
                            FieldVector<FieldVector<int, dimRange>, dimDomain>& permutations,
                            const Quadrature1dType& /*quadrature*/,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<RangeType, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        apply_inverse_eigenvectors(QR[0], tau[0], permutations[0], values[ii][0][0], char_values[ii]);

      // reconstruction in x direction
      FieldVector<RangeType, 2> reconstructed_values;
      // quadrature rule containing left and right interface points
      const auto left_and_right_boundary_point = left_right_quadrature();
      slope_reconstruction(char_values, reconstructed_values, left_and_right_boundary_point);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        // convert back to non-characteristic variables
        eigenvectors[0].mv(reconstructed_values[ii], value);
        auto quadrature_point = FieldVector<DomainFieldType, dimDomain - 1>();
        reconstructed_values_map.insert(
            std::make_pair(intersections[ii].geometryInInside().global(quadrature_point), value));
      } // ii
    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             JacobianRangeType& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, ret[0], param);
    }

    static void get_eigenvectors(const JacobianRangeType& jacobian,
                                 FieldVector<MatrixType, dimDomain>& eigenvectors,
                                 FieldVector<MatrixType, dimDomain>& QR,
                                 FieldVector<RangeType, dimDomain>& tau,
                                 FieldVector<FieldVector<int, dimRange>, dimDomain>& permutations)
    {
      static XT::Common::Configuration eigensolver_options = create_eigensolver_options();
      const auto eigensolver = EigenSolverType(jacobian, &eigensolver_options);
      eigenvectors[0] = eigensolver.real_eigenvectors();
      QR[0] = eigenvectors[0];
      std::fill(permutations[0].begin(), permutations[0].end(), 0.);
      XT::LA::qr(QR[0], tau[0], permutations[0]);
    } // ... get_eigenvectors(...)
  }; // struct helper<1,...

  template <class anything>
  struct helper<3, anything>
  {
    static void reconstruct(size_t dd,
                            const ValuesType& values,
                            FieldVector<MatrixType, dimDomain>& eigenvectors,
                            FieldVector<MatrixType, dimDomain>& QR,
                            FieldVector<RangeType, dimDomain>& tau,
                            FieldVector<FieldVector<int, dimRange>, dimDomain>& permutations,
                            const Quadrature1dType& quadrature,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y- and z-direction accordingly. For that purpose, define new
      // coordinates x', y', z'. First reorder x, y, z to x', y', z' and convert to x'-characteristic variables.
      // Reordering is done such that the indices in char_values are in the order z', y', x'
      size_t curr_dir = dd;
      FieldVector<FieldVector<FieldVector<RangeType, stencil_size>, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            if (dd == 0)
              apply_inverse_eigenvectors(
                  QR[curr_dir], tau[curr_dir], permutations[curr_dir], values[ii][jj][kk], char_values[kk][jj][ii]);
            else if (dd == 1)
              apply_inverse_eigenvectors(
                  QR[curr_dir], tau[curr_dir], permutations[curr_dir], values[ii][jj][kk], char_values[ii][kk][jj]);
            else if (dd == 2)
              apply_inverse_eigenvectors(
                  QR[curr_dir], tau[curr_dir], permutations[curr_dir], values[ii][jj][kk], char_values[jj][ii][kk]);
          }
        }
      }

      // reconstruction in x' direction
      // first index: z' direction
      // second index: dimension 2 for left/right interface
      // third index: y' direction
      FieldVector<FieldVector<FieldVector<RangeType, stencil_size>, 2>, stencil_size> x_reconstructed_values;
      std::vector<RangeType> result(2);
      const auto left_and_right_boundary_point = left_right_quadrature();
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          slope_reconstruction(char_values[kk][jj], result, left_and_right_boundary_point);
          x_reconstructed_values[kk][0][jj] = result[0];
          x_reconstructed_values[kk][1][jj] = result[1];
        }
      }

      // convert from x'-characteristic variables to y'-characteristic variables
      size_t next_dir = (dd + 1) % 3;
      RangeType tmp_vec;
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t ii = 0; ii < 2; ++ii) {
          for (size_t jj = 0; jj < stencil_size; ++jj) {
            tmp_vec = x_reconstructed_values[kk][ii][jj];
            eigenvectors[curr_dir].mv(tmp_vec, x_reconstructed_values[kk][ii][jj]);
            tmp_vec = x_reconstructed_values[kk][ii][jj];
            apply_inverse_eigenvectors(
                QR[next_dir], tau[next_dir], permutations[next_dir], tmp_vec, x_reconstructed_values[kk][ii][jj]);
          }
        }
      }

      // reconstruction in y' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      // third index: z' direction
      const auto& num_quad_points = quadrature.size();
      FieldVector<std::vector<FieldVector<RangeType, stencil_size>>, 2> y_reconstructed_values(
          (std::vector<FieldVector<RangeType, stencil_size>>(num_quad_points)));
      result.resize(num_quad_points);
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t ii = 0; ii < 2; ++ii) {
          slope_reconstruction(x_reconstructed_values[kk][ii], result, quadrature);
          for (size_t jj = 0; jj < num_quad_points; ++jj)
            y_reconstructed_values[ii][jj][kk] = result[jj];
        } // ii
      } // kk

      // convert from y'-characteristic variables to z'-characteristic variables
      curr_dir = next_dir;
      next_dir = (dd + 2) % 3;
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            tmp_vec = y_reconstructed_values[ii][jj][kk];
            eigenvectors[curr_dir].mv(tmp_vec, y_reconstructed_values[ii][jj][kk]);
            tmp_vec = y_reconstructed_values[ii][jj][kk];
            apply_inverse_eigenvectors(
                QR[next_dir], tau[next_dir], permutations[next_dir], tmp_vec, y_reconstructed_values[ii][jj][kk]);
          }
        }
      }

      // reconstruction in z' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      // third index: quadrature_points in z' direction
      FieldVector<std::vector<std::vector<RangeType>>, 2> reconstructed_values(
          std::vector<std::vector<RangeType>>(num_quad_points, std::vector<RangeType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < num_quad_points; ++jj)
          slope_reconstruction(y_reconstructed_values[ii][jj], reconstructed_values[ii][jj], quadrature);

      // convert coordinates on face to local entity coordinates and store
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < num_quad_points; ++kk) {
            // convert back to non-characteristic variables
            eigenvectors[next_dir].mv(reconstructed_values[ii][jj][kk], tmp_vec);
            IntersectionLocalCoordType quadrature_point = {quadrature[jj].position(), quadrature[kk].position()};
            reconstructed_values_map.insert(
                std::make_pair(intersections[2 * dd + ii].geometryInInside().global(quadrature_point), tmp_vec));
          } // kk
        } // jj
      } // ii

    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             JacobianRangeType& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, ret, param);
    }

    static void get_eigenvectors(const JacobianRangeType& jacobian,
                                 FieldVector<MatrixType, dimDomain>& eigenvectors,
                                 FieldVector<MatrixType, dimDomain>& QR,
                                 FieldVector<RangeType, dimDomain>& tau,
                                 FieldVector<FieldVector<int, dimRange>, dimDomain>& permutations)
    {
      static XT::Common::Configuration eigensolver_options = create_eigensolver_options();
      for (size_t ii = 0; ii < dimDomain; ++ii) {
        const auto eigensolver = EigenSolverType(jacobian[ii], eigensolver_options);
        eigenvectors[ii] = eigensolver.real_eigenvectors();
        QR[ii] = eigenvectors[ii];
        std::fill(permutations[ii].begin(), permutations[ii].end(), 0.);
        XT::LA::qr(&(QR[ii][0][0]), dimRange, dimRange, &(permutations[ii][0]), &(tau[ii][0]));
      }
    } // ... get_eigenvectors(...)
  }; // struct helper<3,...>

  static XT::Common::Configuration create_eigensolver_options()
  {
    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[1]);
    // XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
    eigensolver_options["assert_eigendecomposition"] = "1e-6";
    eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
    eigensolver_options["disable_checks"] = "true";
    XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
    matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
    matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
    eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
    return eigensolver_options;
  } // ... create_eigensolver_options()

  // Calculate A^{-1} x, where we have a QR decomposition with column pivoting A = QRP^T.
  // A^{-1} x = (QRP^T)^{-1} x = P R^{-1} Q^T x
  static void apply_inverse_eigenvectors(const MatrixType& QR,
                                         const RangeType& tau,
                                         const FieldVector<int, dimRange>& permutations,
                                         const RangeType& x,
                                         RangeType& ret)
  {
    RangeType work;
    XT::LA::solve_qr_factorized(QR, tau, permutations, ret, x, &work);
  }

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   std::vector<RangeType>& result,
                                   const Quadrature1dType& quadrature)
  {
    assert(stencil_size == 3);
    std::vector<DomainFieldType> points(quadrature.size());
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      points[ii] = quadrature[ii].position();
    assert(result.size() == points.size());
    const auto& u_left = cell_values[0];
    const auto& u_entity = cell_values[1];
    const auto& u_right = cell_values[2];
    const auto slope_left = u_entity - u_left;
    const auto slope_right = u_right - u_entity;
    auto slope_center = u_right - u_left;
    slope_center *= 0.5;
    const auto slope = internal::ChooseLimiter<slope_limiter>::limit(slope_left, slope_right, slope_center);
    std::fill(result.begin(), result.end(), u_entity);
    for (size_t ii = 0; ii < points.size(); ++ii) {
      auto slope_scaled = slope;
      slope_scaled *= points[ii] - 0.5;
      result[ii] += slope_scaled;
    }
  }

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   FieldVector<RangeType, 2>& result,
                                   const Quadrature1dType& quadrature)
  {
    std::vector<RangeType> result_vec(2);
    slope_reconstruction(cell_values, result_vec, quadrature);
    std::copy(result_vec.begin(), result_vec.end(), result.begin());
  }

  const std::vector<RangeType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const Quadrature1dType quadrature_;
  ReconstructedFunctionType& reconstructed_function_;
  //  static thread_local std::unique_ptr<JacobianRangeType> jacobian_;
  //  static thread_local std::unique_ptr<FieldVector<SparseMatrixType, dimDomain>> eigenvectors_;
  //  static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> Q_;
  //  static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> R_;
  //  static thread_local FieldVector<FieldVector<size_t, dimRange>, dimDomain> permutations_;

  // work around gcc bug 66944
  static std::unique_ptr<JacobianRangeType>& jacobian()
  {
    static thread_local std::unique_ptr<JacobianRangeType> jacobian_;
    return jacobian_;
  }

  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& eigenvectors()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> eigenvectors_;
    return eigenvectors_;
  }

  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& QR()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> QR_;
    return QR_;
  }

  static FieldVector<FieldVector<int, dimRange>, dimDomain>& permutations()
  {
    static thread_local FieldVector<FieldVector<int, dimRange>, dimDomain> permutations_;
    return permutations_;
  }

  static std::atomic<size_t> initialization_count_;
  static thread_local size_t local_initialization_count_;
  static bool is_instantiated_;
}; // class LocalReconstructionFvOperator

template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
constexpr std::array<size_t, 3>
    LocalReconstructionFvOperator<AnalyticalFluxType, BoundaryValueType, GridLayerType, slope_limiter>::stencil;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<FieldVector<
//    typename LocalReconstructionFvOperator<GridLayerType,
//                                           AnalyticalFluxType,
//                                           BoundaryValueType,
//                                           polOrder,
//                                           slope_limiter>::SparseMatrixType,
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
//        dimDomain>>
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
//        eigenvectors_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<FieldVector<
//    typename LocalReconstructionFvOperator<GridLayerType,
//                                           AnalyticalFluxType,
//                                           BoundaryValueType,
//                                           polOrder,
//                                           slope_limiter>::CscSparseMatrixType,
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
//        dimDomain>>
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::Q_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<FieldVector<
//    typename LocalReconstructionFvOperator<GridLayerType,
//                                           AnalyticalFluxType,
//                                           BoundaryValueType,
//                                           polOrder,
//                                           slope_limiter>::CscSparseMatrixType,
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
//        dimDomain>>
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::R_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local FieldVector<FieldVector<size_t,
//                                     LocalReconstructionFvOperator<GridLayerType,
//                                                                   AnalyticalFluxType,
//                                                                   BoundaryValueType,
//                                                                   polOrder,
//                                                                   slope_limiter>::dimRange>,
//                         LocalReconstructionFvOperator<GridLayerType,
//                                                       AnalyticalFluxType,
//                                                       BoundaryValueType,
//                                                       polOrder,
//                                                       slope_limiter>::dimDomain>
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
//        permutations_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<typename LocalReconstructionFvOperator<GridLayerType,
//                                                                    AnalyticalFluxType,
//                                                                    BoundaryValueType,
//                                                                    polOrder,
//                                                                    slope_limiter>::JacobianRangeType>
//    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
//        jacobian_;

template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
std::atomic<size_t> LocalReconstructionFvOperator<AnalyticalFluxType, BoundaryValueType, GridLayerType, slope_limiter>::
    initialization_count_(1);

template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
thread_local size_t LocalReconstructionFvOperator<AnalyticalFluxType, BoundaryValueType, GridLayerType, slope_limiter>::
    local_initialization_count_(0);

template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
bool LocalReconstructionFvOperator<AnalyticalFluxType, BoundaryValueType, GridLayerType, slope_limiter>::
    is_instantiated_(false);


#if 0
template <class AnalyticalFluxType, class BoundaryValueType, class GridLayerType, SlopeLimiters slope_limiter>
class LocalReconstructionFvOperator : public XT::Grid::Functor::Codim0<GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  static constexpr size_t dimDomain = GridLayerType::ctype;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t polOrder = 1;
  static constexpr size_t stencil_size = 2 * polOrder + 1;
  static constexpr std::array<size_t, 3> stencil = {
      {2 * polOrder + 1, dimDomain > 1 ? 2 * polOrder + 1 : 1, dimDomain > 2 ? 2 * polOrder + 1 : 1}};
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;
  typedef typename XT::LA::CommonSparseOrDenseMatrixCsr<RangeFieldType> SparseMatrixType;
  typedef typename XT::LA::CommonSparseOrDenseMatrixCsc<RangeFieldType> CscSparseMatrixType;
  typedef typename XT::LA::EigenSolver<MatrixType> EigenSolverType;
  typedef typename XT::LA::EigenSolverOptions<MatrixType> EigenSolverOptionsType;
  typedef typename XT::LA::MatrixInverterOptions<MatrixType> MatrixInverterOptionsType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> Quadrature1dType;
  typedef typename GridLayerType::Intersection IntersectionType;
  typedef FieldVector<IntersectionType, 2 * dimDomain> IntersectionVectorType;
  typedef FieldVector<FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>, stencil[0]> ValuesType;
  typedef typename GridLayerType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  typedef typename Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimRange>, dimDomain> JacobianRangeType;
  typedef typename XT::Grid::BoundaryInfo<IntersectionType> BoundaryInfoType;

public:
  explicit LocalReconstructionFvOperator(
      const std::vector<RangeType>& source_values,
      const AnalyticalFluxType& analytical_flux,
      const BoundaryValueType& boundary_values,
      const GridLayerType& grid_layer,
      const XT::Common::Parameter& param,
      const Quadrature1dType& quadrature_1d,
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(grid_layer)
    , param_(param)
    , quadrature_1d_(quadrature_1d)
    , reconstructed_values_(reconstructed_values)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    param_.set("boundary", {0.});
    is_instantiated_ = true;
  }

  ~LocalReconstructionFvOperator()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    // get cell averages on stencil
    FieldVector<int, 3> offsets(0);
    const RangeType nan(std::numeric_limits<double>::quiet_NaN());
    ValuesType values(
        (FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>(FieldVector<RangeType, stencil[2]>(nan))));
    StencilIterator::apply(source_values_, boundary_values_, values, entity, grid_layer_, -1, offsets);
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity))
      intersections[intersection.indexInInside()] = intersection;
    // get jacobians
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_values_[entity_index];
    if (local_initialization_count_ != initialization_count_) {
      if (!jacobian())
        jacobian() = XT::Common::make_unique<JacobianRangeType>();
      if (!eigenvectors()) {
        eigenvectors() = XT::Common::make_unique<FieldVector<SparseMatrixType, dimDomain>>();
        Q() = XT::Common::make_unique<FieldVector<CscSparseMatrixType, dimDomain>>(
            CscSparseMatrixType(dimRange, dimRange));
        R() = XT::Common::make_unique<FieldVector<CscSparseMatrixType, dimDomain>>(
            CscSparseMatrixType(dimRange, dimRange));
      }
      const auto& u_entity = values[stencil[0] / 2][stencil[1] / 2][stencil[2] / 2];
      helper<dimDomain>::get_jacobian(analytical_flux_, u_entity, *(jacobian()), param_, entity);
      helper<dimDomain>::get_eigenvectors(*(jacobian()), *(eigenvectors()), *(Q()), *(R()), permutations());
      if (analytical_flux_.is_linear()) {
        jacobian() = nullptr;
        ++local_initialization_count_;
      }
    }

    for (size_t dd = 0; dd < dimDomain; ++dd)
      helper<dimDomain>::reconstruct(dd,
                                     values,
                                     *(eigenvectors()),
                                     *(Q()),
                                     *(R()),
                                     permutations(),
                                     quadrature_1d_,
                                     reconstructed_values_map,
                                     intersections);
  } // void apply_local(...)

  static void reset()
  {
    ++initialization_count_;
  }

private:
  // quadrature rule containing left and right interface points
  static Quadrature1dType left_right_quadrature()
  {
    Quadrature1dType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper;

  template <class anything>
  struct helper<1, anything>
  {
    static void reconstruct(size_t /*dd*/,
                            const ValuesType& values,
                            FieldVector<SparseMatrixType, dimDomain>& eigenvectors,
                            FieldVector<CscSparseMatrixType, dimDomain>& Q,
                            FieldVector<CscSparseMatrixType, dimDomain>& R,
                            FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations,
                            const Quadrature1dType& /*quadrature*/,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<RangeType, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        apply_inverse_eigenvectors(Q[0], R[0], permutations[0], values[ii][0][0], char_values[ii]);

      // reconstruction in x direction
      FieldVector<RangeType, 2> reconstructed_values;
      // quadrature rule containing left and right interface points
      const auto left_and_right_boundary_point = left_right_quadrature();
      slope_reconstruction(char_values, reconstructed_values, left_and_right_boundary_point);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        // convert back to non-characteristic variables
        eigenvectors[0].mv(reconstructed_values[ii], value);
        auto quadrature_point = FieldVector<DomainFieldType, dimDomain - 1>();
        reconstructed_values_map.insert(
            std::make_pair(intersections[ii].geometryInInside().global(quadrature_point), value));
      } // ii
    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             JacobianRangeType& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, ret[0], param);
    }

    static void get_eigenvectors(const JacobianRangeType& jacobian,
                                 FieldVector<SparseMatrixType, dimDomain>& eigenvectors,
                                 FieldVector<CscSparseMatrixType, dimDomain>& QR,
                                 FieldVector<RangeType, dimDomain>& tau,
                                 FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations)
    {
      static XT::Common::Configuration eigensolver_options = create_eigensolver_options();
      const auto eigensolver = EigenSolverType(jacobian, eigensolver_options);
      thread_local auto eigenvectors_dense = XT::Common::make_unique<MatrixType>(eigensolver.real_eigenvectors());
      *eigenvectors_dense = eigensolver.real_eigenvectors();
      eigenvectors[0] = SparseMatrixType(*eigenvectors_dense, true);
      XT::LA::qr(QR[0], tau[0], permutations[0]);
      R[0] = CscSparseMatrixType(*eigenvectors_dense, true, size_t(0));
      Q[0] = CscSparseMatrixType(*Qdense, true, size_t(0));
    } // ... get_eigenvectors(...)
  }; // struct helper<1,...

  static XT::Common::Configuration create_eigensolver_options()
  {
    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[0]);
    eigensolver_options["assert_eigendecomposition"] = "1e-6";
    eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
    XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
    matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
    matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
    eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
    return eigensolver_options;
  } // ... create_eigensolver_options()


  template <class anything>
  struct helper<2, anything>
  {

    static void reconstruct(size_t dd,
                            const ValuesType& values,
                            FieldVector<SparseMatrixType, dimDomain>& eigenvectors,
                            FieldVector<CscSparseMatrixType, dimDomain>& Q,
                            FieldVector<CscSparseMatrixType, dimDomain>& R,
                            FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations,
                            const Quadrature1dType& quadrature,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y-direction accordingly. For that purpose, define new
      // coordinates x', y'. First reorder x, y to x', y' and convert to x'-characteristic variables.
      // Reordering is done such that the indices in char_values are in the order y', x'.
      size_t curr_dir = dd;
      FieldVector<FieldVector<RangeType, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          if (dd == 0)
            apply_inverse_eigenvectors(
                Q[curr_dir], R[curr_dir], permutations[curr_dir], values[ii][jj][0], char_values[jj][ii]);
          else if (dd == 1)
            apply_inverse_eigenvectors(
                Q[curr_dir], R[curr_dir], permutations[curr_dir], values[ii][jj][0], char_values[ii][jj]);
        }
      }

      // reconstruction in x' direction
      // first index: dimension 2 for left/right interface
      // second index: y' direction
      FieldVector<FieldVector<RangeType, stencil_size>, 2> x_reconstructed_values;
      std::vector<RangeType> result(2);
      const auto left_and_right_boundary_point = left_right_quadrature();
      for (size_t jj = 0; jj < stencil_size; ++jj) {
        slope_reconstruction(char_values[jj], result, left_and_right_boundary_point);
        x_reconstructed_values[0][jj] = result[0];
        x_reconstructed_values[1][jj] = result[1];
      }

      // convert from x'-characteristic variables to y'-characteristic variables
      size_t next_dir = (dd + 1) % 2;
      RangeType tmp_vec;
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          tmp_vec = x_reconstructed_values[ii][jj];
          eigenvectors[curr_dir].mv(tmp_vec, x_reconstructed_values[ii][jj]);
          tmp_vec = x_reconstructed_values[ii][jj];
          apply_inverse_eigenvectors(
              Q[next_dir], R[next_dir], permutations[next_dir], tmp_vec, x_reconstructed_values[ii][jj]);
        }
      }
      curr_dir = next_dir;

      const auto& num_quad_points = quadrature.size();
      // reconstruction in y' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      FieldVector<std::vector<RangeType>, 2> reconstructed_values((std::vector<RangeType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        slope_reconstruction(x_reconstructed_values[ii], reconstructed_values[ii], quadrature);

      // convert coordinates on face to local entity coordinates and store
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          // convert back to non-characteristic variables
          eigenvectors[curr_dir].mv(reconstructed_values[ii][jj], tmp_vec);
          auto quadrature_point = quadrature[jj].position();
          reconstructed_values_map.insert(
              std::make_pair(intersections[2 * dd + ii].geometryInInside().global(quadrature_point), tmp_vec));
        } // jj
      } // ii
    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             JacobianRangeType& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      helper<3, anything>::get_jacobian(analytical_flux, u, ret, param, entity);
    }

    static void get_eigenvectors(const JacobianRangeType& jacobian,
                                 FieldVector<SparseMatrixType, dimDomain>& eigenvectors,
                                 FieldVector<CscSparseMatrixType, dimDomain>& Q,
                                 FieldVector<CscSparseMatrixType, dimDomain>& R,
                                 FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations)
    {
      helper<3, anything>::get_eigenvectors(jacobian, eigenvectors, Q, R, permutations);
    }
  }; // helper<2,...

  template <class anything>
  struct helper<3, anything>
  {

    static void reconstruct(size_t dd,
                            const ValuesType& values,
                            FieldVector<SparseMatrixType, dimDomain>& eigenvectors,
                            FieldVector<CscSparseMatrixType, dimDomain>& Q,
                            FieldVector<CscSparseMatrixType, dimDomain>& R,
                            FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations,
                            const Quadrature1dType& quadrature,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      // We always treat dd as the x-direction and set y- and z-direction accordingly. For that purpose, define new
      // coordinates x', y', z'. First reorder x, y, z to x', y', z' and convert to x'-characteristic variables.
      // Reordering is done such that the indices in char_values are in the order z', y', x'
      size_t curr_dir = dd;
      FieldVector<FieldVector<FieldVector<RangeType, stencil_size>, stencil_size>, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            if (dd == 0)
              apply_inverse_eigenvectors(
                  Q[curr_dir], R[curr_dir], permutations[curr_dir], values[ii][jj][kk], char_values[kk][jj][ii]);
            else if (dd == 1)
              apply_inverse_eigenvectors(
                  Q[curr_dir], R[curr_dir], permutations[curr_dir], values[ii][jj][kk], char_values[ii][kk][jj]);
            else if (dd == 2)
              apply_inverse_eigenvectors(
                  Q[curr_dir], R[curr_dir], permutations[curr_dir], values[ii][jj][kk], char_values[jj][ii][kk]);
          }
        }
      }

      // reconstruction in x' direction
      // first index: z' direction
      // second index: dimension 2 for left/right interface
      // third index: y' direction
      FieldVector<FieldVector<FieldVector<RangeType, stencil_size>, 2>, stencil_size> x_reconstructed_values;
      std::vector<RangeType> result(2);
      const auto left_and_right_boundary_point = left_right_quadrature();
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          slope_reconstruction(char_values[kk][jj], result, left_and_right_boundary_point);
          x_reconstructed_values[kk][0][jj] = result[0];
          x_reconstructed_values[kk][1][jj] = result[1];
        }
      }

      // convert from x'-characteristic variables to y'-characteristic variables
      size_t next_dir = (dd + 1) % 3;
      RangeType tmp_vec;
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t ii = 0; ii < 2; ++ii) {
          for (size_t jj = 0; jj < stencil_size; ++jj) {
            tmp_vec = x_reconstructed_values[kk][ii][jj];
            eigenvectors[curr_dir].mv(tmp_vec, x_reconstructed_values[kk][ii][jj]);
            tmp_vec = x_reconstructed_values[kk][ii][jj];
            apply_inverse_eigenvectors(
                Q[next_dir], R[next_dir], permutations[next_dir], tmp_vec, x_reconstructed_values[kk][ii][jj]);
          }
        }
      }

      // reconstruction in y' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      // third index: z' direction
      const auto& num_quad_points = quadrature.size();
      FieldVector<std::vector<FieldVector<RangeType, stencil_size>>, 2> y_reconstructed_values(
          (std::vector<FieldVector<RangeType, stencil_size>>(num_quad_points)));
      result.resize(num_quad_points);
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t ii = 0; ii < 2; ++ii) {
          slope_reconstruction(x_reconstructed_values[kk][ii], result, quadrature);
          for (size_t jj = 0; jj < num_quad_points; ++jj)
            y_reconstructed_values[ii][jj][kk] = result[jj];
        } // ii
      } // kk

      // convert from y'-characteristic variables to z'-characteristic variables
      curr_dir = next_dir;
      next_dir = (dd + 2) % 3;
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < stencil_size; ++kk) {
            tmp_vec = y_reconstructed_values[ii][jj][kk];
            eigenvectors[curr_dir].mv(tmp_vec, y_reconstructed_values[ii][jj][kk]);
            tmp_vec = y_reconstructed_values[ii][jj][kk];
            apply_inverse_eigenvectors(
                Q[next_dir], R[next_dir], permutations[next_dir], tmp_vec, y_reconstructed_values[ii][jj][kk]);
          }
        }
      }

      // reconstruction in z' direction
      // first index: left/right interface
      // second index: quadrature_points in y' direction
      // third index: quadrature_points in z' direction
      FieldVector<std::vector<std::vector<RangeType>>, 2> reconstructed_values(
          std::vector<std::vector<RangeType>>(num_quad_points, std::vector<RangeType>(num_quad_points)));
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < num_quad_points; ++jj)
          slope_reconstruction(y_reconstructed_values[ii][jj], reconstructed_values[ii][jj], quadrature);

      // convert coordinates on face to local entity coordinates and store
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < num_quad_points; ++kk) {
            // convert back to non-characteristic variables
            eigenvectors[next_dir].mv(reconstructed_values[ii][jj][kk], tmp_vec);
            IntersectionLocalCoordType quadrature_point = {quadrature[jj].position(), quadrature[kk].position()};
            reconstructed_values_map.insert(
                std::make_pair(intersections[2 * dd + ii].geometryInInside().global(quadrature_point), tmp_vec));
          } // kk
        } // jj
      } // ii
    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             JacobianRangeType& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, ret, param);
    }

    static void get_eigenvectors(const JacobianRangeType& jacobian,
                                 FieldVector<SparseMatrixType, dimDomain>& eigenvectors,
                                 FieldVector<CscSparseMatrixType, dimDomain>& Q,
                                 FieldVector<CscSparseMatrixType, dimDomain>& R,
                                 FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations)
    {
      static XT::Common::Configuration eigensolver_options = create_eigensolver_options();


      for (size_t ii = 0; ii < dimDomain; ++ii) {
        const auto eigensolver = EigenSolverType(jacobian[ii], eigensolver_options);
        thread_local auto eigenvectors_dense = XT::Common::make_unique<MatrixType>(eigensolver.real_eigenvectors());
        *eigenvectors_dense = eigensolver.real_eigenvectors();
        eigenvectors[ii] = SparseMatrixType(*eigenvectors_dense, true);
        thread_local std::unique_ptr<MatrixType> Qdense = XT::Common::make_unique<MatrixType>(0.);
        XT::LA::qr_decomposition(*eigenvectors_dense, permutations[ii], *Qdense);
        R[ii] = CscSparseMatrixType(*eigenvectors_dense, true, size_t(0));
        Q[ii] = CscSparseMatrixType(*Qdense, true, size_t(0));
      } // ii
    } // ... get_eigenvectors(...)

  }; // struct helper<3, ...

  // Calculate A^{-1} x, where we have a QR decomposition with column pivoting A = QRP^T.
  // A^{-1} x = (QRP^T)^{-1} x = P R^{-1} Q^T x
  static void apply_inverse_eigenvectors(const CscSparseMatrixType& Q,
                                         const CscSparseMatrixType& R,
                                         const FieldVector<size_t, dimRange>& permutations,
                                         const RangeType& x,
                                         RangeType& ret)
  {
    // calculate ret = Q^T x
    Q.mtv(x, ret);
    // calculate ret = R^{-1} ret
    auto y = ret;
    XT::LA::solve_upper_triangular(R, ret, y);
    // calculate ret = P ret
    y = ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[permutations[ii]] = y[ii];
  }

  class StencilIterator
  {
  public:
    static constexpr size_t stencil_x = stencil[0];
    static constexpr size_t stencil_y = stencil[1];
    static constexpr size_t stencil_z = stencil[2];

    static void apply(const std::vector<RangeType>& source_values,
                      const BoundaryValueType& boundary_values,
                      const BoundaryInfoType& boundary_info,
                      FieldVector<FieldVector<FieldVector<RangeType, stencil_z>, stencil_y>, stencil_x>& values,
                      const EntityType& entity,
                      const GridLayerType& grid_layer,
                      const int direction,
                      FieldVector<int, 3> offsets)
    {
      const auto& entity_index = grid_layer.indexSet().index(entity);
      values[stencil_x / 2 + offsets[0]][stencil_y / 2 + offsets[1]][stencil_z / 2 + offsets[2]] =
          source_values[entity_index];
      std::vector<int> boundary_dirs;
      for (const auto& intersection : Dune::intersections(grid_layer, entity)) {
        const auto& intersection_index = intersection.indexInInside();
        if (!end_of_stencil(intersection_index, offsets)) {
          auto new_offsets = offsets;
          if (intersection.boundary() && !intersection.neighbor()) {
            boundary_dirs.push_back(intersection_index);
            RangeType boundary_value;
            const auto& boundary_type = boundary_info.type(intersection);
            if (boundary_type == dirichlet_boundary_) {
              boundary_value = std::get<0>(boundary_values)
                                   ->local_function(entity)
                                   ->evaluate(entity.geometry().local(intersection.geometry().center()));
            } else if (boundary_type == reflecting_boundary_) {
              boundary_value = std::get<1>(boundary_values).evaluate(intersection, source_values[entity_index]);
            }
            while (!end_of_stencil(intersection_index, new_offsets)) {
              walk(intersection_index, new_offsets);
              values[stencil_x / 2 + new_offsets[0]][stencil_y / 2 + new_offsets[1]][stencil_z / 2 + new_offsets[2]] =
                  boundary_value;
            }
          } else if (intersection.neighbor() && direction_allowed(direction, intersection_index)) {
            const auto& outside = intersection.outside();
            walk(intersection_index, new_offsets);
            StencilIterator::apply(source_values,
                                   boundary_values,
                                   boundary_info,
                                   values,
                                   outside,
                                   grid_layer,
                                   intersection_index,
                                   new_offsets);
          }
        } // if (!end_of_stencil(...))
      } // intersections

      // TODO: improve multiple boundary handling, currently everything is filled with the boundary value in the first
      // direction, do not use NaNs
      assert(boundary_dirs.size() <= 3);
      if (boundary_dirs.size() > 1) {
        walk(boundary_dirs[0], offsets);
        const auto& boundary_value =
            values[stencil_x / 2 + offsets[0]][stencil_y / 2 + offsets[1]][stencil_z / 2 + offsets[2]];
        for (auto& values_yz : values)
          for (auto& values_z : values_yz)
            for (auto& value : values_z)
              if (std::isnan(value[0]))
                value = boundary_value;
      }
    } // void apply(...)

  private:
    static void walk(const int dir, FieldVector<int, 3>& offsets)
    {
      dir % 2 ? offsets[dir / 2]++ : offsets[dir / 2]--;
    }

    // Direction is allowed if end of stencil is not reached and direction is not visited by another iterator.
    // Iterators never change direction, they may only spawn new iterators in the directions that have a higher
    // index (i.e. iterators that walk in x direction will spawn iterators going in y and z direction,
    // iterators going in y direction will only spawn iterators in z-direction and z iterators only walk
    // without emitting new iterators).
    static bool direction_allowed(const int dir, const int new_dir)
    {
      return (polOrder > 0 && dir == -1) || new_dir == dir || new_dir / 2 > dir / 2;
    }

    static bool end_of_stencil(const int dir, const FieldVector<int, 3>& offsets)
    {
      return (polOrder == 0 || (dir != -1 && size_t(std::abs(offsets[dir / 2])) >= stencil[dir / 2] / 2));
    }
  }; // class StencilIterator<...

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   std::vector<RangeType>& result,
                                   const Quadrature1dType& quadrature)
  {
    assert(stencil_size == 3);
    std::vector<DomainFieldType> points(quadrature.size());
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      points[ii] = quadrature[ii].position();
    assert(result.size() == points.size());
    const auto& u_left = cell_values[0];
    const auto& u_entity = cell_values[1];
    const auto& u_right = cell_values[2];
    const auto slope_left = u_entity - u_left;
    const auto slope_right = u_right - u_entity;
    auto slope_center = u_right - u_left;
    slope_center *= 0.5;
    const auto slope = internal::ChooseLimiter<slope_limiter>::limit(slope_left, slope_right, slope_center);
    std::fill(result.begin(), result.end(), u_entity);
    for (size_t ii = 0; ii < points.size(); ++ii) {
      auto slope_scaled = slope;
      slope_scaled *= points[ii] - 0.5;
      result[ii] += slope_scaled;
    }
  }

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   FieldVector<RangeType, 2>& result,
                                   const Quadrature1dType& quadrature)
  {
    std::vector<RangeType> result_vec(2);
    slope_reconstruction(cell_values, result_vec, quadrature);
    std::copy(result_vec.begin(), result_vec.end(), result.begin());
  }

  const std::vector<RangeType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const Quadrature1dType quadrature_1d_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values_;

  // work around gcc bug 66944
  static std::unique_ptr<JacobianRangeType>& jacobian()
  {
    static thread_local std::unique_ptr<JacobianRangeType> jacobian_;
    return jacobian_;
  }

  static std::unique_ptr<FieldVector<SparseMatrixType, dimDomain>>& eigenvectors()
  {
    static thread_local std::unique_ptr<FieldVector<SparseMatrixType, dimDomain>> eigenvectors_;
    return eigenvectors_;
  }

  static std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>>& Q()
  {
    static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> Q_;
    return Q_;
  }

  static std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>>& R()
  {
    static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> R_;
    return R_;
  }

  static FieldVector<FieldVector<size_t, dimRange>, dimDomain>& permutations()
  {
    static thread_local FieldVector<FieldVector<size_t, dimRange>, dimDomain> permutations_;
    return permutations_;
  }

  static std::atomic<size_t> initialization_count_;
  static thread_local size_t local_initialization_count_;
  static bool is_instantiated_;
}; // class LocalReconstructionFvOperator

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
constexpr std::array<size_t, 3>
    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
        stencil;

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
std::atomic<size_t>
    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
        initialization_count_(1);

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
thread_local size_t
    LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
        local_initialization_count_(0);

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
bool LocalReconstructionFvOperator<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder, slope_limiter>::
    is_instantiated_(false);
#endif


template <class SourceType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
class LocalReconstructionFvOperatorBlocked
    : public XT::Grid::Functor::Codim0<typename SourceType::SpaceType::GridLayerType>
{
  // stencil is (i-r, i+r) in all dimensions, where r = polOrder + 1
  typedef typename SourceType::SpaceType SpaceType;
  typedef typename SpaceType::GridLayerType GridLayerType;
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  static constexpr size_t dimRange = BoundaryValueType::dimRange;
  static constexpr size_t block_size = (dimDomain == 1) ? 2 : 4;
  static constexpr size_t num_blocks = dimRange / block_size;
  static constexpr size_t stencil_size = 2 * polOrder + 1;
  static constexpr std::array<size_t, 3> stencil = {
      {2 * polOrder + 1, dimDomain > 1 ? 2 * polOrder + 1 : 1, dimDomain > 2 ? 2 * polOrder + 1 : 1}};
  typedef typename GridLayerType::template Codim<0>::Entity EntityType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename BoundaryValueType::DomainType DomainType;
  typedef typename BoundaryValueType::DomainFieldType DomainFieldType;
  typedef typename BoundaryValueType::RangeType RangeType;
  typedef typename BoundaryValueType::RangeFieldType RangeFieldType;
  typedef FieldVector<RangeFieldType, block_size> LocalRangeType;
  typedef FieldMatrix<RangeFieldType, block_size, block_size> LocalMatrixType;
  typedef FieldVector<LocalMatrixType, num_blocks> MatrixType;
  typedef typename XT::LA::EigenSolver<LocalMatrixType> EigenSolverType;
  typedef typename XT::LA::EigenSolverOptions<LocalMatrixType> EigenSolverOptionsType;
  typedef typename XT::LA::MatrixInverterOptions<LocalMatrixType> MatrixInverterOptionsType;
  typedef Dune::QuadratureRule<DomainFieldType, 1> QuadratureType;
  typedef typename GridLayerType::Intersection IntersectionType;
  typedef FieldVector<IntersectionType, 2 * dimDomain> IntersectionVectorType;
  typedef FieldVector<FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>, stencil[0]> ValuesType;
  typedef typename GridLayerType::Intersection::Geometry::LocalCoordinate IntersectionLocalCoordType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  typedef typename Dune::FieldVector<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>
      JacobianRangeType;
  typedef typename XT::Grid::BoundaryInfo<IntersectionType> BoundaryInfoType;

public:
  explicit LocalReconstructionFvOperatorBlocked(
      SourceType& source,
      const std::vector<RangeType>& source_values,
      const AnalyticalFluxType& analytical_flux,
      const BoundaryValueType& boundary_values,
      const XT::Common::Parameter& param,
      const bool is_linear,
      const QuadratureType& quadrature,
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values)
    : source_values_(source_values)
    , analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , grid_layer_(source.space().grid_layer())
    , param_(param)
    , is_linear_(is_linear)
    , quadrature_(quadrature)
    , reconstructed_values_(reconstructed_values)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    param_.set("boundary", {0.});
    is_instantiated_ = true;
  }

  ~LocalReconstructionFvOperatorBlocked()
  {
    is_instantiated_ = false;
  }

  void apply_local(const EntityType& entity)
  {
    // get cell averages on stencil
    FieldVector<int, 3> offsets(0);
    ValuesType values((FieldVector<FieldVector<RangeType, stencil[2]>, stencil[1]>(
        FieldVector<RangeType, stencil[2]>(RangeType(std::numeric_limits<double>::quiet_NaN())))));
    StencilIterator::apply(source_values_, boundary_values_, values, entity, grid_layer_, -1, offsets);
    // get intersections
    FieldVector<typename GridLayerType::Intersection, 2 * dimDomain> intersections;
    for (const auto& intersection : Dune::intersections(grid_layer_, entity))
      intersections[intersection.indexInInside()] = intersection;
    // get jacobians
    const auto& entity_index = grid_layer_.indexSet().index(entity);
    auto& reconstructed_values_map = reconstructed_values_[entity_index];
    static thread_local FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain> tau_;
    if (local_initialization_count_ != initialization_count_) {
      if (!jacobian())
        jacobian() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
      if (!eigenvectors()) {
        eigenvectors() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
        QR() = XT::Common::make_unique<FieldVector<MatrixType, dimDomain>>();
      }
      const auto& u_entity = values[stencil[0] / 2][stencil[1] / 2][stencil[2] / 2];
      helper<dimDomain>::get_jacobian(analytical_flux_, u_entity, *(jacobian()), param_, entity);
      get_eigenvectors(*(jacobian()), *(eigenvectors()), *(QR()), tau_, permutations());
      if (is_linear_) {
        jacobian() = nullptr;
        ++local_initialization_count_;
      }
    }

    for (size_t dd = 0; dd < dimDomain; ++dd)
      helper<dimDomain>::reconstruct(dd,
                                     values,
                                     *(eigenvectors()),
                                     *(QR()),
                                     tau_,
                                     permutations(),
                                     quadrature_,
                                     reconstructed_values_map,
                                     intersections);
  } // void apply_local(...)

  static void reset()
  {
    ++initialization_count_;
  }

private:
  // quadrature rule containing left and right interface points
  static QuadratureType left_right_quadrature()
  {
    QuadratureType ret;
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(0., 0.5));
    ret.push_back(Dune::QuadraturePoint<DomainFieldType, 1>(1., 0.5));
    return ret;
  }

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper;

  static void mv(const MatrixType& Q, const RangeType& x, RangeType& ret)
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < block_size; ++ll)
        for (size_t mm = 0; mm < block_size; ++mm)
          ret[offset + ll] += Q[jj][ll][mm] * x[offset + mm];
    } // jj
  }

  static void mtv(const MatrixType& Q, const RangeType& x, RangeType& ret)
  {
    std::fill(ret.begin(), ret.end(), 0.);
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto offset = block_size * jj;
      for (size_t ll = 0; ll < block_size; ++ll)
        for (size_t mm = 0; mm < block_size; ++mm)
          ret[offset + ll] += Q[jj][mm][ll] * x[offset + mm];
    } // jj
  }

  template <class anything>
  struct helper<1, anything>
  {
    static void reconstruct(size_t /*dd*/,
                            const ValuesType& values,
                            FieldVector<MatrixType, dimDomain>& eigenvectors,
                            FieldVector<MatrixType, dimDomain>& QR,
                            FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain>& tau,
                            FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations,
                            const QuadratureType& /*quadrature*/,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<RangeType, stencil_size> char_values;
      for (size_t ii = 0; ii < stencil_size; ++ii)
        apply_inverse_eigenvectors(QR[0], tau[0], permutations[0], values[ii][0][0], char_values[ii]);

      // reconstruction in x direction
      FieldVector<RangeType, 2> reconstructed_values;
      // quadrature rule containing left and right interface points
      const auto left_and_right_boundary_point = left_right_quadrature();
      slope_reconstruction(char_values, reconstructed_values, left_and_right_boundary_point);

      // convert coordinates on face to local entity coordinates and store
      RangeType value;
      for (size_t ii = 0; ii < 2; ++ii) {
        // convert back to non-characteristic variables
        mv(eigenvectors[0], reconstructed_values[ii], value);
        auto quadrature_point = FieldVector<DomainFieldType, dimDomain - 1>();
        reconstructed_values_map.insert(
            std::make_pair(intersections[ii].geometryInInside().global(quadrature_point), value));
      } // ii
    } // static void reconstruct(...)

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             FieldVector<MatrixType, dimDomain>& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      thread_local auto jac = std::make_unique<FieldMatrix<RangeFieldType, dimRange, dimRange>>();
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, *jac, param);
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        const auto offset = jj * block_size;
        for (size_t ll = 0; ll < block_size; ++ll)
          for (size_t mm = 0; mm < block_size; ++mm)
            ret[0][jj][ll][mm] = (*jac)[offset + ll][offset + mm];
      } // jj
    } // static void get_jacobian(...)
  }; // struct helper<1,...>

  template <class anything>
  struct helper<3, anything>
  {
    static void reconstruct(size_t dd,
                            const StencilType& stencil,
                            FieldVector<MatrixType, dimDomain>& eigenvectors,
                            FieldVector<MatrixType, dimDomain>& QR,
                            FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain>& tau,
                            FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations,
                            const QuadratureType& quadrature,
                            std::map<DomainType, RangeType, XT::Common::FieldVectorLess>& reconstructed_values_map,
                            const IntersectionVectorType& intersections)
    {
      FieldVector<size_t, dimDomain> sizes = stencil_sizes_;

      // Transform values on stencil to characteristic variables of the first coordinate direction dd
      std::vector<RangeType> reconstructed_values(stencil.num_values());
      for (size_t ii = 0; ii < stencil.num_values(); ++ii)
        apply_inverse_eigenvectors(dd, *(stencil.values()[ii]), reconstructed_values[ii]);

      size_t curr_dir, last_dir;
      RangeType tmp_value;
      for (size_t dir = 0; dir < dimDomain; ++dir) {
        curr_dir = (dd + dir) % dimDomain;
        // Transform to characteristic variables of the current reconstruction direction.
        if (dir > 0) {
          for (auto& value : reconstructed_values) {
            apply_eigenvectors(last_dir, value, tmp_value);
            apply_inverse_eigenvectors(curr_dir, tmp_value, value);
          }
        } // if (dir > 0)
        // perform the actual reconstruction
        const auto& curr_quadrature = dir > 0 ? quadrature_[curr_dir] : left_right_quadrature();
        slope_reconstruction(curr_dir, curr_quadrature, reconstructed_values);
        sizes[curr_dir] = curr_quadrature.size();
        last_dir = curr_dir;
      }
      // convert back to non-characteristic variables
      for (auto& value : reconstructed_values) {
        tmp_value = value;
        apply_eigenvectors(last_dir, tmp_value, value);
      }

      // Convert coordinates on face to local entity coordinates and store reconstructed values
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < num_quad_points; ++jj) {
          for (size_t kk = 0; kk < num_quad_points; ++kk) {
            // convert back to non-characteristic variables
            IntersectionLocalCoordType quadrature_point = {quadrature[jj].position(), quadrature[kk].position()};
            reconstructed_values_map.insert(
                std::make_pair(intersections[2 * dd + ii].geometryInInside().global(quadrature_point), tmp_vec));
          } // kk
        } // jj
      } // ii

      // reconstruction in x' direction
      // calculate values on left and right intersection in x' direction
      size_t remaining_size = stencil.num_values() / stencil_sizes[curr_dir];
      std::vector<RangeType> x_reconstructed_values(2 * remaining_size);
      std::vector<RangeType> result(2);
      const auto left_and_right_boundary_point = left_right_quadrature();
      for (size_t kk = 0; kk < stencil_size; ++kk) {
        for (size_t jj = 0; jj < stencil_size; ++jj) {
          slope_reconstruction(char_values[kk][jj], result, left_and_right_boundary_point);
          x_reconstructed_values[kk][0][jj] = result[0];
          x_reconstructed_values[kk][1][jj] = result[1];
        }
      }

    } // reconstruct

    static void get_jacobian(const AnalyticalFluxType& analytical_flux,
                             const StateRangeType& u,
                             FieldVector<MatrixType, dimDomain>& ret,
                             const XT::Common::Parameter& param,
                             const EntityType& entity)
    {
      thread_local auto jac =
          std::make_unique<FieldVector<FieldMatrix<RangeFieldType, dimRange, dimRange>, dimDomain>>();
      const DomainType x_in_inside_coords = entity.geometry().local(entity.geometry().center());
      const auto local_func = analytical_flux.local_function(entity);
      local_func->partial_u(x_in_inside_coords, u, *jac, param);
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        for (size_t jj = 0; jj < num_blocks; ++jj) {
          const auto offset = jj * block_size;
          for (size_t ll = 0; ll < block_size; ++ll)
            for (size_t mm = 0; mm < block_size; ++mm)
              ret[dd][jj][ll][mm] = (*jac)[dd][offset + ll][offset + mm];
        } // jj
      } // dd
    } // static void get_jacobian(...)
  }; // struct helper<3,...>

  static void
  get_eigenvectors(const FieldVector<MatrixType, dimDomain>& jacobian,
                   FieldVector<MatrixType, dimDomain>& eigenvectors,
                   FieldVector<MatrixType, dimDomain>& QR,
                   FieldVector<FieldVector<LocalRangeType, num_blocks>, dimDomain>& tau,
                   FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations)
  {
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t jj = 0; jj < num_blocks; ++jj) {
        static XT::Common::Configuration eigensolver_options = create_eigensolver_options();
        const auto eigensolver = EigenSolverType(jacobian[dd][jj], &eigensolver_options);
        eigenvectors[dd][jj] = eigensolver.real_eigenvectors();
        QR[dd][jj] = eigenvectors[dd][jj];
        XT::LA::qr(QR[dd][jj], tau[dd][jj], permutations[dd][jj]);
      } // jj
    } // dd
  } // ... get_eigenvectors(...)

  static XT::Common::Configuration create_eigensolver_options()
  {
    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[0]);
    // XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
    eigensolver_options["assert_eigendecomposition"] = "1e-6";
    eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
    eigensolver_options["disable_checks"] = "true";
    XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
    matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
    matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
    eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
    return eigensolver_options;
  } // ... create_eigensolver_options()


  // Calculate A^{-1} x, where we have a QR decomposition with column pivoting A = QRP^T.
  // A^{-1} x = (QRP^T)^{-1} x = P R^{-1} Q^T x
  static void apply_inverse_eigenvectors(const MatrixType& QR,
                                         const FieldVector<LocalRangeType, num_blocks>& tau,
                                         const FieldVector<FieldVector<int, block_size>, num_blocks>& permutations,
                                         const RangeType& x,
                                         RangeType& ret)
  {
    LocalRangeType work;
    LocalRangeType tmp_ret, tmp_x;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      for (size_t ll = 0; ll < block_size; ++ll)
        tmp_x[ll] = x[jj * block_size + ll];
      XT::LA::solve_qr_factorized(QR[jj], tau[jj], permutations[jj], tmp_ret, tmp_x, &work);
      for (size_t ll = 0; ll < block_size; ++ll)
        ret[jj * block_size + ll] = tmp_ret[ll];
    }
  }

  class StencilIterator
  {
  public:
    static constexpr size_t stencil_x = stencil[0];
    static constexpr size_t stencil_y = stencil[1];
    static constexpr size_t stencil_z = stencil[2];

    static void apply(const std::vector<RangeType>& source_values,
                      const BoundaryValueType& boundary_values,
                      FieldVector<FieldVector<FieldVector<RangeType, stencil_z>, stencil_y>, stencil_x>& values,
                      const EntityType& entity,
                      const GridLayerType& grid_layer,
                      const int direction,
                      FieldVector<int, 3> offsets)
    {
      const auto& entity_index = grid_layer.indexSet().index(entity);
      values[stencil_x / 2 + offsets[0]][stencil_y / 2 + offsets[1]][stencil_z / 2 + offsets[2]] =
          source_values[entity_index];
      std::vector<int> boundary_dirs;
      for (const auto& intersection : Dune::intersections(grid_layer, entity)) {
        const auto& intersection_index = intersection.indexInInside();
        if (!end_of_stencil(intersection_index, offsets)) {
          auto new_offsets = offsets;
          if (intersection.boundary() && !intersection.neighbor()) {
            boundary_dirs.push_back(intersection_index);
            RangeType boundary_value = boundary_values.local_function(entity)->evaluate(
                intersection, entity.geometry().local(intersection.geometry().center()), source_values[entity_index]);
            while (!end_of_stencil(intersection_index, new_offsets)) {
              walk(intersection_index, new_offsets);
              values[stencil_x / 2 + new_offsets[0]][stencil_y / 2 + new_offsets[1]][stencil_z / 2 + new_offsets[2]] =
                  boundary_value;
            }
          } else if (intersection.neighbor() && direction_allowed(direction, intersection_index)) {
            const auto& outside = intersection.outside();
            walk(intersection_index, new_offsets);
            StencilIterator::apply(
                source_values, boundary_values, values, outside, grid_layer, intersection_index, new_offsets);
          }
        } // if (!end_of_stencil(...))
      } // intersections

      // TODO: improve multiple boundary handling, currently everything is filled with the boundary value in the first
      // direction, do not use NaNs
      assert(boundary_dirs.size() <= 3);
      if (boundary_dirs.size() > 1) {
        walk(boundary_dirs[0], offsets);
        const auto& boundary_value =
            values[stencil_x / 2 + offsets[0]][stencil_y / 2 + offsets[1]][stencil_z / 2 + offsets[2]];
        for (auto& values_yz : values)
          for (auto& values_z : values_yz)
            for (auto& value : values_z)
              if (std::isnan(value[0]))
                value = boundary_value;
      }
    } // void apply(...)

  private:
    static void walk(const int dir, FieldVector<int, 3>& offsets)
    {
      dir % 2 ? offsets[dir / 2]++ : offsets[dir / 2]--;
    }

    // Direction is allowed if end of stencil is not reached and direction is not visited by another iterator.
    // Iterators never change direction, they may only spawn new iterators in the directions that have a higher
    // index (i.e. iterators that walk in x direction will spawn iterators going in y and z direction,
    // iterators going in y direction will only spawn iterators in z-direction and z iterators only walk
    // without emitting new iterators).
    static bool direction_allowed(const int dir, const int new_dir)
    {
      return (polOrder > 0 && dir == -1) || new_dir == dir || new_dir / 2 > dir / 2;
    }

    static bool end_of_stencil(const int dir, const FieldVector<int, 3>& offsets)
    {
      return (polOrder == 0 || (dir != -1 && size_t(std::abs(offsets[dir / 2])) >= stencil[dir / 2] / 2));
    }
  }; // class StencilIterator<...

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   std::vector<RangeType>& result,
                                   const QuadratureType& quadrature)
  {
    assert(stencil_size == 3);
    std::vector<DomainFieldType> points(quadrature.size());
    for (size_t ii = 0; ii < quadrature.size(); ++ii)
      points[ii] = quadrature[ii].position();
    assert(result.size() == points.size());
    const auto& u_left = cell_values[0];
    const auto& u_entity = cell_values[1];
    const auto& u_right = cell_values[2];
    const auto slope_left = u_entity - u_left;
    const auto slope_right = u_right - u_entity;
    auto slope_center = u_right - u_left;
    slope_center *= 0.5;
    const auto slope = internal::ChooseLimiter<slope_limiter>::limit(slope_left, slope_right, slope_center);
    std::fill(result.begin(), result.end(), u_entity);
    for (size_t ii = 0; ii < points.size(); ++ii) {
      auto slope_scaled = slope;
      slope_scaled *= points[ii] - 0.5;
      result[ii] += slope_scaled;
    }
  }

  static void slope_reconstruction(const FieldVector<RangeType, stencil_size>& cell_values,
                                   FieldVector<RangeType, 2>& result,
                                   const QuadratureType& quadrature)
  {
    std::vector<RangeType> result_vec(2);
    slope_reconstruction(cell_values, result_vec, quadrature);
    std::copy(result_vec.begin(), result_vec.end(), result.begin());
  }

  const std::vector<RangeType>& source_values_;
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const GridLayerType& grid_layer_;
  XT::Common::Parameter param_;
  const bool is_linear_;
  const QuadratureType quadrature_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>& reconstructed_values_;
  //  static thread_local std::unique_ptr<JacobianRangeType> jacobian_;
  //  static thread_local std::unique_ptr<FieldVector<SparseMatrixType, dimDomain>> eigenvectors_;
  //  static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> Q_;
  //  static thread_local std::unique_ptr<FieldVector<CscSparseMatrixType, dimDomain>> R_;
  //  static thread_local FieldVector<FieldVector<size_t, dimRange>, dimDomain> permutations_;

  // work around gcc bug 66944
  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& jacobian()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> jacobian_;
    return jacobian_;
  }

  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& eigenvectors()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> eigenvectors_;
    return eigenvectors_;
  }

  static std::unique_ptr<FieldVector<MatrixType, dimDomain>>& QR()
  {
    static thread_local std::unique_ptr<FieldVector<MatrixType, dimDomain>> QR_;
    return QR_;
  }

  static FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain>& permutations()
  {
    static thread_local FieldVector<FieldVector<FieldVector<int, block_size>, num_blocks>, dimDomain> permutations_;
    return permutations_;
  }

  static std::atomic<size_t> initialization_count_;
  static thread_local size_t local_initialization_count_;
  static bool is_instantiated_;
}; // class LocalReconstructionFvOperatorBlocked

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
constexpr std::array<size_t, 3> LocalReconstructionFvOperatorBlocked<GridLayerType,
                                                                     AnalyticalFluxType,
                                                                     BoundaryValueType,
                                                                     polOrder,
                                                                     slope_limiter>::stencil;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<FieldVector<
//    typename LocalReconstructionFvOperatorBlocked<GridLayerType,
//                                           AnalyticalFluxType,
//                                           BoundaryValueType,
//                                           polOrder,
//                                           slope_limiter>::SparseMatrixType,
//    LocalReconstructionFvOperatorBlocked<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder,
//    slope_limiter>::
//        dimDomain>>
//    LocalReconstructionFvOperatorBlocked<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder,
//    slope_limiter>::
//        eigenvectors_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<FieldVector<
//    typename LocalReconstructionFvOperatorBlocked<GridLayerType,
//                                           AnalyticalFluxType,
//                                           BoundaryValueType,
//                                           polOrder,
//                                           slope_limiter>::CscSparseMatrixType,
//    LocalReconstructionFvOperatorBlocked<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder,
//    slope_limiter>::
//        dimDomain>>
//    LocalReconstructionFvOperatorBlocked<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder,
//    slope_limiter>::QR_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local FieldVector<FieldVector<size_t,
//                                     LocalReconstructionFvOperatorBlocked<GridLayerType,
//                                                                   AnalyticalFluxType,
//                                                                   BoundaryValueType,
//                                                                   polOrder,
//                                                                   slope_limiter>::dimRange>,
//                         LocalReconstructionFvOperatorBlocked<GridLayerType,
//                                                       AnalyticalFluxType,
//                                                       BoundaryValueType,
//                                                       polOrder,
//                                                       slope_limiter>::dimDomain>
//    LocalReconstructionFvOperatorBlocked<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder,
//    slope_limiter>::
//        permutations_;

// template <class GridLayerType,
//          class AnalyticalFluxType,
//          class BoundaryValueType,
//          size_t polOrder,
//          SlopeLimiters slope_limiter>
// thread_local std::unique_ptr<typename LocalReconstructionFvOperatorBlocked<GridLayerType,
//                                                                    AnalyticalFluxType,
//                                                                    BoundaryValueType,
//                                                                    polOrder,
//                                                                    slope_limiter>::JacobianRangeType>
//    LocalReconstructionFvOperatorBlocked<GridLayerType, AnalyticalFluxType, BoundaryValueType, polOrder,
//    slope_limiter>::
//        jacobian_;

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
std::atomic<size_t> LocalReconstructionFvOperatorBlocked<GridLayerType,
                                                         AnalyticalFluxType,
                                                         BoundaryValueType,
                                                         polOrder,
                                                         slope_limiter>::initialization_count_(1);

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
thread_local size_t LocalReconstructionFvOperatorBlocked<GridLayerType,
                                                         AnalyticalFluxType,
                                                         BoundaryValueType,
                                                         polOrder,
                                                         slope_limiter>::local_initialization_count_(0);

template <class GridLayerType,
          class AnalyticalFluxType,
          class BoundaryValueType,
          size_t polOrder,
          SlopeLimiters slope_limiter>
bool LocalReconstructionFvOperatorBlocked<GridLayerType,
                                          AnalyticalFluxType,
                                          BoundaryValueType,
                                          polOrder,
                                          slope_limiter>::is_instantiated_(false);


template <class AnalyticalFluxImp, class BoundaryValueImp, SlopeLimiters slope_limiter, class Traits>
class LinearReconstructionOperator;


namespace internal {


template <class AnalyticalFluxImp, class BoundaryValueImp, SlopeLimiters slope_lim = SlopeLimiters::minmod>
struct LinearReconstructionOperatorTraits
{
  using AnalyticalFluxType = AnalyticalFluxImp;
  using BoundaryValueType = BoundaryValueImp;
  static const SlopeLimiters slope_limiter = slope_lim;
  using DomainFieldType = typename BoundaryValueType::DomainFieldType;
  using FieldType = DomainFieldType;
  using JacobianType = NoJacobian;
  static constexpr size_t dimDomain = BoundaryValueType::dimDomain;
  using ProductQuadratureType = QuadratureRule<DomainFieldType, dimDomain - 1>;
  using Quadrature1dType = Dune::QuadratureRule<DomainFieldType, 1>;
  using derived_type = LinearReconstructionOperator<AnalyticalFluxType,
                                                    BoundaryValueType,
                                                    slope_limiter,
                                                    LinearReconstructionOperatorTraits>;
};


} // namespace internal


template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          SlopeLimiters slope_limiter = SlopeLimiters::minmod,
          class Traits =
              internal::LinearReconstructionOperatorTraits<AnalyticalFluxImp, BoundaryValueImp, slope_limiter>>
class LinearReconstructionOperator : public OperatorInterface<Traits>
{
public:
  using AnalyticalFluxType = typename Traits::AnalyticalFluxType;
  using BoundaryValueType = typename Traits::BoundaryValueType;
  using Quadrature1dType = typename Traits::Quadrature1dType;
  using ProductQuadratureType = typename Traits::ProductQuadratureType;
  using DomainFieldType = typename Traits::DomainFieldType;
  static constexpr size_t dimDomain = Traits::dimDomain;

  LinearReconstructionOperator(const AnalyticalFluxType& analytical_flux,
                               const BoundaryValueType& boundary_values,
                               const Quadrature1dType& quadrature_1d = default_1d_quadrature<DomainFieldType>(1))
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , quadrature_1d_(quadrature_1d)
    , product_quadrature_(product_quadrature_on_intersection<DomainFieldType, dimDomain>(quadrature_1d))
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    static_assert(is_reconstructed_localizable_function<RangeType>::value,
                  "RangeType has to be derived from ReconstructedLocalizableFunction!");
    // evaluate cell averages
    const auto& grid_layer = source.space().grid_layer();
    const auto& index_set = grid_layer.indexSet();
    std::vector<typename SourceType::RangeType> source_values(index_set.size(0));
    for (const auto& entity : Dune::elements(grid_layer)) {
      const auto& entity_index = index_set.index(entity);
      const auto& local_source = source.local_function(entity);
      source_values[entity_index] = local_source->evaluate(entity.geometry().local(entity.geometry().center()));
    }

    // do reconstruction
    auto local_reconstruction_operator = LocalReconstructionFvOperator<AnalyticalFluxType,
                                                                       BoundaryValueType,
                                                                       typename SourceType::SpaceType::GridLayerType,
                                                                       slope_limiter>(
        source_values, analytical_flux_, boundary_values_, grid_layer, param, quadrature_1d_, range);
    auto walker = XT::Grid::Walker<typename SourceType::SpaceType::GridLayerType>(grid_layer);
    walker.append(local_reconstruction_operator);
    walker.walk(true);
  } // void apply(...)

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueType& boundary_values_;
  const Quadrature1dType quadrature_1d_;
  const ProductQuadratureType product_quadrature_;
}; // class LinearReconstructionOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RECONSTRUCTION_LINEAR_HH
