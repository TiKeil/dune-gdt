// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_RAVIARTTHOMAS_PDELAB_HH
#define DUNE_GDT_PLAYGROUND_SPACES_RAVIARTTHOMAS_PDELAB_HH

#warning This header is deprecated, include <dune/gdt/playground/spaces/rt/pdelab.hh> instead (21.11.2014)!
#include <dune/gdt/playground/spaces/rt/pdelab.hh>

namespace RaviartThomas {


template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class DUNE_DEPRECATED_MSG("Use RT::PdelabBased instead (21.11.2014)!") PdelabBased
    : public RT::PdelabBased<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>
{
public:
  template <class... Args>
  PdelabBased(Args&&... args)
    : RT::PdelabBased<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>(std::forward<Args>(args)...)
  {
  }
};


} // namespace RaviartThomas

#endif // DUNE_GDT_PLAYGROUND_SPACES_RAVIARTTHOMAS_PDELAB_HH
