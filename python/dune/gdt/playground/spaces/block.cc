// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/gdt/shared.hh>

#include "block.hh"


#define DUNE_GDT_SPACES_BLOCK_BIND(_m, _GRID, _s_type, _s_backend, _p)                                                 \
  Dune::GDT::bindings::BlockSpace<Dune::GDT::SpaceProvider<_GRID,                                                      \
                                                           Dune::XT::Grid::Layers::dd_subdomain,                       \
                                                           Dune::GDT::SpaceType::_s_type,                              \
                                                           Dune::GDT::Backends::_s_backend,                            \
                                                           _p,                                                         \
                                                           double,                                                     \
                                                           1,                                                          \
                                                           1>>::bind(_m)


PYBIND11_MODULE(__spaces_block, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  DUNE_GDT_SPACES_BLOCK_BIND(m, GDT_BINDINGS_GRID, dg, gdt, 1);
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.spaces.block");
}
