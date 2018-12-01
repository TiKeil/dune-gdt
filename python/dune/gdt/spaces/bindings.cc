// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#  include <dune/common/parallel/mpihelper.hh>

#  include <dune/pybindxi/pybind11.h>
#  include <dune/pybindxi/stl.h>

#  include <python/dune/xt/common/bindings.hh>
#  include <python/dune/gdt/shared.hh>

#  include "bindings.hh"


PYBIND11_MODULE(__spaces, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

  DUNE_GDT_SPACES_CG_BIND(m);
  DUNE_GDT_SPACES_DG_BIND(m);
  DUNE_GDT_SPACES_FV_BIND(m);
  DUNE_GDT_SPACES_RT_BIND(m);

  add_initialization(m, "dune.gdt.spaces");
}

#endif // HAVE_DUNE_PYBINDXI
