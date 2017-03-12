// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_DG_BINDINGS_HH
#define DUNE_GDT_SPACES_DG_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/xt/grid/grids.bindings.hh>

#include "dg.hh"
#include "interface.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


#define _DEFINE_DG_SPACE(_G, _layer, _backend, _p, _r, _rC)                                                            \
  typedef                                                                                                              \
      typename DgSpaceProvider<_G, XT::Grid::Layers::_layer, ChooseSpaceBackend::_backend, _p, double, _r, _rC>::type  \
          Dg_##_G##_##_layer##_to_##_r##x##_rC##_##p##_p##_##_backend##_Space

#if HAVE_DUNE_FEM
_DEFINE_DG_SPACE(YASP_2D_EQUIDISTANT_OFFSET, leaf, fem, 1, 1, 1);
_DEFINE_DG_SPACE(YASP_2D_EQUIDISTANT_OFFSET, leaf, fem, 2, 1, 1);
_DEFINE_DG_SPACE(YASP_2D_EQUIDISTANT_OFFSET, level, fem, 1, 1, 1);
_DEFINE_DG_SPACE(YASP_2D_EQUIDISTANT_OFFSET, level, fem, 2, 1, 1);
#endif
#if HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
_DEFINE_DG_SPACE(ALU_2D_SIMPLEX_CONFORMING, leaf, fem, 1, 1, 1);
_DEFINE_DG_SPACE(ALU_2D_SIMPLEX_CONFORMING, leaf, fem, 2, 1, 1);
_DEFINE_DG_SPACE(ALU_2D_SIMPLEX_CONFORMING, level, fem, 1, 1, 1);
_DEFINE_DG_SPACE(ALU_2D_SIMPLEX_CONFORMING, level, fem, 2, 1, 1);
#endif
#endif // HAVE_DUNE_ALUGRID

#undef _DEFINE_DG_SPACE


// this is used by other headers
#define DG_SPACE(_G, _layer, _backend, _p, _r, _rC) Dg_##_G##_##_layer##_to_##_r##x##_rC##_##p##_p##_##_backend##_Space


#define DUNE_GDT_SPACES_DG_BIND_FEM(_prefix, _GRID)                                                                    \
  _prefix class SpaceInterface<DG_SPACE(_GRID, leaf, fem, 1, 1, 1)>;                                                   \
  _prefix class SpaceInterface<DG_SPACE(_GRID, level, fem, 1, 1, 1)>;                                                  \
  _prefix class SpaceInterface<DG_SPACE(_GRID, leaf, fem, 2, 1, 1)>;                                                   \
  _prefix class SpaceInterface<DG_SPACE(_GRID, level, fem, 2, 1, 1)>


// these lines have to match the corresponding ones in the .cc source file
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_DG_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET);
#endif

#if HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_DG_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif // HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_DG_BINDINGS_HH