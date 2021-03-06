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

#include <python/dune/gdt/prolongations/bindings.hh>
#include <python/dune/gdt/shared.hh>


#define DUNE_GDT_PROLONGATIONS_BIND(                                                                                   \
    _m, _s_G, _s_gl, _s_backend, _s_type, _s_p, _r, _s_la, _r_G, _r_gl, _r_backend, _r_type, _r_p, _r_la)              \
  Dune::GDT::bindings::prolong<                                                                                        \
      typename Dune::GDT::SpaceProvider<_s_G,                                                                          \
                                        Dune::XT::Grid::Layers::_s_gl,                                                 \
                                        Dune::GDT::SpaceType::_s_type,                                                 \
                                        Dune::GDT::Backends::_s_backend,                                               \
                                        _s_p,                                                                          \
                                        double,                                                                        \
                                        _r,                                                                            \
                                        1>::type,                                                                      \
      typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_s_la>::VectorType,                             \
      typename Dune::GDT::SpaceProvider<_r_G,                                                                          \
                                        Dune::XT::Grid::Layers::_r_gl,                                                 \
                                        Dune::GDT::SpaceType::_r_type,                                                 \
                                        Dune::GDT::Backends::_r_backend,                                               \
                                        _r_p,                                                                          \
                                        double,                                                                        \
                                        _r,                                                                            \
                                        1>::type,                                                                      \
      typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_r_la>::VectorType>::bind(_m)


PYBIND11_MODULE(__prolongations, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.prolongations");

  using G = GDT_BINDINGS_GRID;

  DUNE_GDT_PROLONGATIONS_BIND(m, G, leaf, gdt, dg, 1, 1, istl_dense, G, leaf, gdt, dg, 1, istl_dense);
  DUNE_GDT_PROLONGATIONS_BIND(m, G, leaf, gdt, dg, 1, 1, istl_dense, G, leaf, gdt, dg, 2, istl_dense);
  DUNE_GDT_PROLONGATIONS_BIND(m, G, leaf, gdt, dg, 1, 1, istl_dense, G, leaf, gdt, dg, 3, istl_dense);
  DUNE_GDT_PROLONGATIONS_BIND(m, G, dd_subdomain, gdt, block_dg, 1, 1, istl_dense, G, leaf, gdt, dg, 1, istl_dense);
  DUNE_GDT_PROLONGATIONS_BIND(m, G, dd_subdomain, gdt, block_dg, 1, 1, istl_dense, G, leaf, gdt, dg, 2, istl_dense);
  DUNE_GDT_PROLONGATIONS_BIND(m, G, dd_subdomain, gdt, block_dg, 1, 1, istl_dense, G, leaf, gdt, dg, 3, istl_dense);
}
