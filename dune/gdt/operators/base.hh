// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_BASE_HH
#define DUNE_GDT_OPERATORS_BASE_HH

#include <dune/stuff/la/solver.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


template <class Traits>
class MatrixBased : public AssemblableOperatorInterface<Traits>
{
  typedef AssemblableOperatorInterface<Traits> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::MatrixType;

private:
  typedef Stuff::LA::Solver<MatrixType> LinearSolverType;

public:
  MatrixBased(MatrixType& matrix, const SourceSpaceType& source_space, const RangeSpaceType& range_space,
              const GridViewType& grid_view)
    : matrix_(matrix)
    , source_space_(source_space)
    , range_space_(range_space)
    , grid_view_(grid_view)
    , assembled_(false)
  {
  }

  MatrixBased(MatrixType& matrix, const SourceSpaceType& source_space, const RangeSpaceType& range_space)
    : matrix_(matrix)
    , source_space_(source_space)
    , range_space_(range_space)
    , grid_view_(*(source_space_.grid_view()))
    , assembled_(false)
  {
  }

  MatrixBased(MatrixType& matrix, const SourceSpaceType& source_space)
    : matrix_(matrix)
    , source_space_(source_space)
    , range_space_(source_space_)
    , grid_view_(*(source_space_.grid_view()))
    , assembled_(false)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const RangeSpaceType& range_space() const
  {
    return range_space_;
  }

  const SourceSpaceType& source_space() const
  {
    return source_space_;
  }

  MatrixType& matrix()
  {
    return matrix_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

  virtual void assemble() = 0;

  template <class S, class R>
  void apply(const Stuff::LA::VectorInterface<S>& source, Stuff::LA::VectorInterface<R>& range)
  {
    assemble();
    matrix_.mv(source.as_imp(), range.as_imp());
  } // ... apply(...)

  static std::vector<std::string> invert_options()
  {
    return LinearSolverType::options();
  }

  static Stuff::Common::ConfigTree invert_options(const std::string& type)
  {
    return LinearSolverType::options(type);
  }

  template <class R, class S>
  void apply_inverse(const Stuff::LA::VectorInterface<R>& range, Stuff::LA::VectorInterface<S>& source,
                     const Stuff::Common::ConfigTree& opts)
  {
    assemble();
    LinearSolverType(matrix).apply(range.as_imp(), source.as_imp(), opts);
  } // ... apply_inverse(...)

private:
  MatrixType& matrix_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const GridViewType& grid_view_;
  bool assembled_;
}; // class MatrixBased


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_BASE_HH
