// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_ELLIPTIC_HH
#define DUNE_GDT_OPERATOR_ELLIPTIC_HH

#include <type_traits>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/space/interface.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operator {


// forward, to be used in the traits
template <class DiffusionImp, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp                                                              = typename SourceSpaceImp::GridViewType>
class Elliptic;


template <class DiffusionImp, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp                                                              = typename SourceSpaceImp::GridViewType>
class EllipticTraits
{
  static_assert(std::is_base_of<Stuff::LocalizableFunctionInterface<
                                    typename DiffusionImp::EntityType, typename DiffusionImp::DomainFieldType,
                                    DiffusionImp::dimDomain, typename DiffusionImp::RangeFieldType,
                                    DiffusionImp::dimRange, DiffusionImp::dimRangeCols>,
                                DiffusionImp>::value,
                "DiffusionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::LA::MatrixInterface<typename MatrixImp::Traits>, MatrixImp>::value,
                "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename SourceSpaceImp::Traits>, SourceSpaceImp>::value,
                "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename RangeSpaceImp::Traits>, RangeSpaceImp>::value,
                "RangeSpaceImp has to be derived from SpaceInterface!");

public:
  typedef Elliptic<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp> derived_type;
  typedef MatrixImp MatrixType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef GridViewImp GridViewType;

private:
  typedef DiffusionImp DiffusionType;

public:
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> LocalOperatorType;

private:
  friend class Elliptic<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp>;
}; // class EllipticTraits


template <class DiffusionImp, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp>
class Elliptic : public AssemblableVolumeOperatorBase<EllipticTraits<DiffusionImp, MatrixImp, SourceSpaceImp,
                                                                     RangeSpaceImp, GridViewImp>>
{
public:
  typedef EllipticTraits<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp> Traits;

private:
  typedef AssemblableVolumeOperatorBase<Traits> BaseType;

  typedef typename Traits::DiffusionType DiffusionType;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::GridViewType GridViewType;

  using BaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  Elliptic(const DiffusionType& diffusion, MatrixType& matrix, const SourceSpaceType& source_space,
           const RangeSpaceType& range_space, const GridViewType& grid_view)
    : BaseType(matrix, source_space, range_space, grid_view)
    , local_operator_(diffusion)
  {
  }

  Elliptic(const DiffusionType& diffusion, MatrixType& matrix, const SourceSpaceType& source_space,
           const RangeSpaceType& range_space)
    : BaseType(matrix, source_space, range_space)
    , local_operator_(diffusion)
  {
  }

  Elliptic(const DiffusionType& diffusion, MatrixType& matrix, const SourceSpaceType& source_space)
    : BaseType(matrix, source_space)
    , local_operator_(diffusion)
  {
  }

  virtual const LocalOperatorType& local_operator() const
  {
    return local_operator_;
  }

private:
  const LocalOperatorType local_operator_;
}; // class Elliptic


} // namespace Operator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_ELLIPTIC_HH
