// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_BASE_HH
#define DUNE_GDT_OPERATORS_BASE_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/local/assembler.hh>
#include <dune/gdt/assembler/wrapper.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required for the traits
template <class M,
          class RS,
          class GL = typename RS::GridLayerType,
          class SS = RS,
          class F = typename M::RealType,
          ChoosePattern pt = ChoosePattern::face_and_volume,
          class ORS = RS,
          class OSS = SS>
class MatrixOperatorBase;


namespace internal {


template <class MatrixImp,
          class RangeSpaceImp,
          class GridLayerImp,
          class SourceSpaceImp,
          class FieldImp,
          ChoosePattern pt,
          class OuterRangeSpaceImp,
          class OuterSourceSpaceImp>
class MatrixOperatorBaseTraits
{
  static_assert(XT::LA::is_matrix<MatrixImp>::value, "");
  static_assert(is_space<RangeSpaceImp>::value, "");
  static_assert(is_space<SourceSpaceImp>::value, "");
  static_assert(is_space<OuterRangeSpaceImp>::value, "");
  static_assert(is_space<OuterSourceSpaceImp>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<typename RangeSpaceImp::GridLayerType>,
                             XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "RangeSpaceType and GridLayerType have to match!");
  static_assert(std::is_same<XT::Grid::extract_entity_t<typename SourceSpaceImp::GridLayerType>,
                             XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "SourceSpaceType and GridLayerType have to match!");
  static_assert(std::is_same<XT::Grid::extract_entity_t<typename OuterRangeSpaceImp::GridLayerType>,
                             XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "OuterRangeSpaceImp and GridLayerType have to match!");
  static_assert(std::is_same<XT::Grid::extract_entity_t<typename OuterSourceSpaceImp::GridLayerType>,
                             XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "OuterSourceSpaceImp and GridLayerType have to match!");

public:
  typedef MatrixOperatorBase<MatrixImp,
                             RangeSpaceImp,
                             GridLayerImp,
                             SourceSpaceImp,
                             FieldImp,
                             pt,
                             OuterRangeSpaceImp,
                             OuterSourceSpaceImp>
      derived_type;
  typedef FieldImp FieldType;
  typedef std::unique_ptr<derived_type> JacobianType;
};


} // namespace internal


/**
 * \todo add static checks of dimensions
 * \note Does a const_cast in apply() and apply2(), not sure yet if this is fine.
 * \warning: only apply2(DiscreteFunction, DiscreteFunction) automagically sums over mpi processes
 */
template <class MatrixImp,
          class RangeSpaceImp,
          class GridLayerImp,
          class SourceSpaceImp,
          class FieldImp,
          ChoosePattern pt,
          class OuterRangeSpaceImp,
          class OuterSourceSpaceImp>
class MatrixOperatorBase
    : public OperatorInterface<internal::MatrixOperatorBaseTraits<MatrixImp,
                                                                  RangeSpaceImp,
                                                                  GridLayerImp,
                                                                  SourceSpaceImp,
                                                                  FieldImp,
                                                                  pt,
                                                                  OuterRangeSpaceImp,
                                                                  OuterSourceSpaceImp>>,
      public SystemAssembler<RangeSpaceImp, GridLayerImp, SourceSpaceImp, OuterRangeSpaceImp, OuterSourceSpaceImp>
{
  typedef OperatorInterface<internal::MatrixOperatorBaseTraits<MatrixImp,
                                                               RangeSpaceImp,
                                                               GridLayerImp,
                                                               SourceSpaceImp,
                                                               FieldImp,
                                                               pt,
                                                               OuterRangeSpaceImp,
                                                               OuterSourceSpaceImp>>
      BaseOperatorType;
  typedef SystemAssembler<RangeSpaceImp, GridLayerImp, SourceSpaceImp, OuterRangeSpaceImp, OuterSourceSpaceImp>
      BaseAssemblerType;
  typedef MatrixOperatorBase<MatrixImp,
                             RangeSpaceImp,
                             GridLayerImp,
                             SourceSpaceImp,
                             FieldImp,
                             pt,
                             OuterRangeSpaceImp,
                             OuterSourceSpaceImp>
      ThisType;

public:
  typedef internal::MatrixOperatorBaseTraits<MatrixImp,
                                             RangeSpaceImp,
                                             GridLayerImp,
                                             SourceSpaceImp,
                                             FieldImp,
                                             pt,
                                             OuterRangeSpaceImp,
                                             OuterSourceSpaceImp>
      Traits;
  typedef typename BaseAssemblerType::TestSpaceType RangeSpaceType;
  typedef typename BaseAssemblerType::AnsatzSpaceType SourceSpaceType;
  typedef typename BaseAssemblerType::OuterTestSpaceType OuterRangeSpaceType;
  typedef typename BaseAssemblerType::OuterAnsatzSpaceType OuterSourceSpaceType;
  typedef typename RangeSpaceType::BaseFunctionSetType RangeBaseType;
  typedef typename SourceSpaceType::BaseFunctionSetType SourceBaseType;
  typedef typename OuterRangeSpaceType::BaseFunctionSetType OuterRangeBaseType;
  typedef typename OuterSourceSpaceType::BaseFunctionSetType OuterSourceBaseType;
  typedef XT::LA::SparsityPatternDefault PatternType;
  typedef MatrixImp MatrixType;
  using typename BaseOperatorType::FieldType;
  using typename BaseOperatorType::JacobianType;
  using typename BaseOperatorType::derived_type;
  using typename BaseAssemblerType::GridLayerType;
  using typename BaseAssemblerType::IntersectionType;
  static const constexpr ChoosePattern pattern_type = pt;

private:
  typedef XT::LA::Solver<MatrixType, typename SourceSpaceType::DofCommunicatorType> LinearSolverType;

  template <ChoosePattern pp = ChoosePattern::face_and_volume, bool anything = true>
  struct Compute
  {
    static PatternType
    pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridLayerType& grd_layr)
    {
      return rng_spc.compute_face_and_volume_pattern(grd_layr, src_spc);
    }
  };

  template <bool anything>
  struct Compute<ChoosePattern::volume, anything>
  {
    static PatternType
    pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridLayerType& grd_layr)
    {
      return rng_spc.compute_volume_pattern(grd_layr, src_spc);
    }
  };

  template <bool anything>
  struct Compute<ChoosePattern::face, anything>
  {
    static PatternType
    pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridLayerType& grd_layr)
    {
      return rng_spc.compute_face_pattern(grd_layr, src_spc);
    }
  };

public:
  static PatternType
  pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridLayerType& grd_layr)
  {
    return Compute<pt>::pattern(rng_spc, src_spc, grd_layr);
  }

  template <class R>
  static typename std::enable_if<std::is_same<R, RangeSpaceType>::value && std::is_same<R, SourceSpaceType>::value
                                     && std::is_same<typename R::GridLayerType, GridLayerType>::value,
                                 PatternType>::type
  pattern(const R& rng_spc)
  {
    return pattern(rng_spc, rng_spc);
  }

  template <class R, class S>
  static typename std::enable_if<std::is_same<R, RangeSpaceType>::value && std::is_same<S, SourceSpaceType>::value
                                     && std::is_same<typename R::GridLayerType, GridLayerType>::value,
                                 PatternType>::type
  pattern(const R& rng_spc, const S& src_spc)
  {
    return pattern(rng_spc, src_spc, rng_spc.grid_layer());
  }

  template <class R, class G>
  static typename std::enable_if<std::is_same<R, RangeSpaceType>::value && std::is_same<R, SourceSpaceType>::value
                                     && std::is_same<G, GridLayerType>::value,
                                 PatternType>::type
  pattern(const R& rng_spc, const G& grd_layr)
  {
    return pattern(rng_spc, rng_spc, grd_layr);
  }

  template <class... Args>
  explicit MatrixOperatorBase(MatrixType& mtrx, Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_in_in_(mtrx)
    , matrix_out_out_(matrix_in_in_.access())
    , matrix_in_out_(matrix_in_in_.access())
    , matrix_out_in_(matrix_in_in_.access())
  {
    if (matrix_in_in_.access().rows() != this->range_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix.rows(): " << matrix_in_in_.access().rows() << "\n"
                                   << "range_space().mapper().size(): "
                                   << this->range_space().mapper().size());
    if (matrix_in_in_.access().cols() != this->source_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix.cols(): " << matrix_in_in_.access().cols() << "\n"
                                   << "source_space().mapper().size(): "
                                   << this->source_space().mapper().size());
  } // MatrixOperatorBase(...)

  template <class... Args>
  explicit MatrixOperatorBase(MatrixType& mtrx_in_in,
                              MatrixType& mtrx_out_out,
                              MatrixType& mtrx_in_out,
                              MatrixType& mtrx_out_in,
                              Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_in_in_(mtrx_in_in)
    , matrix_out_out_(mtrx_out_out)
    , matrix_in_out_(mtrx_in_out)
    , matrix_out_in_(mtrx_out_in)
  {
    if (matrix_in_in_.access().rows() != this->range_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_in_in.rows(): " << matrix_in_in_.access().rows() << "\n"
                                         << "range_space().mapper().size(): "
                                         << this->range_space().mapper().size());
    if (matrix_in_in_.access().cols() != this->source_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_in_in.cols(): " << matrix_in_in_.access().cols() << "\n"
                                         << "source_space().mapper().size(): "
                                         << this->source_space().mapper().size());
    if (matrix_out_out_.access().rows() != this->outer_range_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_out_out.rows(): " << matrix_out_out_.access().rows() << "\n"
                                           << "outer_range_space().mapper().size(): "
                                           << this->outer_range_space().mapper().size());
    if (matrix_out_out_.access().cols() != this->outer_source_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_out_out.cols(): " << matrix_out_out_.access().cols() << "\n"
                                           << "outer_source_space().mapper().size(): "
                                           << this->outer_source_space().mapper().size());
    if (matrix_in_out_.access().rows() != this->range_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_in_out.rows(): " << matrix_in_out_.access().rows() << "\n"
                                          << "range_space().mapper().size(): "
                                          << this->range_space().mapper().size());
    if (matrix_in_out_.access().cols() != this->outer_source_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_in_out.cols(): " << matrix_in_out_.access().cols() << "\n"
                                          << "outer_source_space().mapper().size(): "
                                          << this->outer_source_space().mapper().size());
    if (matrix_out_in_.access().rows() != this->outer_range_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_out_in.rows(): " << matrix_out_in_.access().rows() << "\n"
                                          << "outer_range_space().mapper().size(): "
                                          << this->outer_range_space().mapper().size());
    if (matrix_out_in_.access().cols() != this->source_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix_out_in.cols(): " << matrix_out_in_.access().cols() << "\n"
                                          << "source_space().mapper().size(): "
                                          << this->source_space().mapper().size());
  } // MatrixOperatorBase(...)

  /// \todo Guard against copy and move ctor (Args = ThisType)!
  template <class... Args>
  explicit MatrixOperatorBase(Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_in_in_(new MatrixType(this->range_space().mapper().size(),
                                   this->source_space().mapper().size(),
                                   pattern(this->range_space(), this->source_space(), this->grid_layer())))
    , matrix_out_out_(matrix_in_in_.access())
    , matrix_in_out_(matrix_in_in_.access())
    , matrix_out_in_(matrix_in_in_.access())
  {
  }

  /// \sa SystemAssembler
  MatrixOperatorBase(const ThisType& other) = delete;
  MatrixOperatorBase(ThisType&& source) = delete;
  MatrixOperatorBase(ThisType& other) = delete; // <- b.c. of the too perfect forwarding ctor

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const MatrixType& matrix() const
  {
    return matrix_in_in_.access();
  }

  MatrixType& matrix()
  {
    return matrix_in_in_.access();
  }

  const SourceSpaceType& source_space() const
  {
    return this->ansatz_space();
  }

  const RangeSpaceType& range_space() const
  {
    return this->test_space();
  }

  const OuterSourceSpaceType& outer_source_space() const
  {
    return this->outer_ansatz_space();
  }

  const OuterRangeSpaceType& outer_range_space() const
  {
    return this->outer_test_space();
  }

  using BaseAssemblerType::append;

  ThisType& append(
      const LocalVolumeTwoFormInterface<RangeBaseType, SourceBaseType, FieldType>& local_volume_twoform,
      const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    this->append(local_volume_twoform, matrix_in_in_.access(), where);
    return *this;
  }

  ThisType& append(const LocalCouplingTwoFormInterface<RangeBaseType,
                                                       IntersectionType,
                                                       SourceBaseType,
                                                       OuterRangeBaseType,
                                                       OuterSourceBaseType,
                                                       FieldType>& local_coupling_twoform,
                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where =
                       new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridLayerType>())
  {
    this->append(local_coupling_twoform,
                 matrix_in_in_.access(),
                 matrix_out_out_.access(),
                 matrix_in_out_.access(),
                 matrix_out_in_.access(),
                 where);
    return *this;
  }

  ThisType& append(const LocalBoundaryTwoFormInterface<RangeBaseType, IntersectionType, SourceBaseType, FieldType>&
                       local_boundary_twoform,
                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where =
                       new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridLayerType>())
  {
    this->append(local_boundary_twoform, matrix_in_in_.access(), where);
    return *this;
  }

  template <class S, class R>
  void apply(const XT::LA::VectorInterface<S>& source,
             XT::LA::VectorInterface<R>& range,
             const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    const_cast<ThisType&>(*this).assemble();
    matrix().mv(source.as_imp(), range.as_imp());
  }

  template <class S, class R>
  void apply(const ConstDiscreteFunction<SourceSpaceType, S>& source,
             DiscreteFunction<RangeSpaceType, R>& range,
             const Dune::XT::Common::Parameter& param = {}) const
  {
    apply(source.vector(), range.vector(), param);
  }

  template <class R, class S>
  FieldType apply2(const XT::LA::VectorInterface<R>& range,
                   const XT::LA::VectorInterface<S>& source,
                   const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    const_cast<ThisType&>(*this).assemble();
    auto tmp = range.copy();
    matrix().mv(source.as_imp(source), tmp);
    return range.dot(tmp);
  }

  template <class R, class S>
  FieldType apply2(const ConstDiscreteFunction<RangeSpaceType, R>& range,
                   const ConstDiscreteFunction<SourceSpaceType, S>& source,
                   const Dune::XT::Common::Parameter& param = {}) const
  {
    auto ret = apply2(range.vector(), source.vector(), param);
    return range.space().grid_layer().grid().comm().sum(ret);
  }

  template <class SourceType>
  JacobianType jacobian(const SourceType& /*source*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    return JacobianType(matrix(), range_space(), source_space());
  }

  template <class SourceType>
  void
  jacobian(const SourceType& /*source*/, JacobianType& jac, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    jac->matrix() = matrix();
  }

  using BaseOperatorType::apply_inverse;

  template <class R, class S>
  void apply_inverse(const XT::LA::VectorInterface<R>& range,
                     XT::LA::VectorInterface<S>& source,
                     const XT::Common::Configuration& opts,
                     const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    this->assemble();
    LinearSolverType(matrix(), source_space().dof_communicator()).apply(range.as_imp(), source.as_imp(), opts);
  }

  template <class R, class S>
  void apply_inverse(const ConstDiscreteFunction<RangeSpaceType, R>& range,
                     ConstDiscreteFunction<SourceSpaceType, S>& source,
                     const XT::Common::Configuration& opts,
                     const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    apply_inverse(range.vector(), source.vector(), opts);
  }

  std::vector<std::string> invert_options() const
  {
    return LinearSolverType::types();
  }

  XT::Common::Configuration invert_options(const std::string& type) const
  {
    return LinearSolverType::options(type);
  }

protected:
  using BaseAssemblerType::codim0_functors_;
  using BaseAssemblerType::codim1_functors_;

private:
  Dune::XT::Common::StorageProvider<MatrixType> matrix_in_in_;
  Dune::XT::Common::StorageProvider<MatrixType> matrix_out_out_;
  Dune::XT::Common::StorageProvider<MatrixType> matrix_in_out_;
  Dune::XT::Common::StorageProvider<MatrixType> matrix_out_in_;
}; // class MatrixOperatorBase


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_BASE_HH