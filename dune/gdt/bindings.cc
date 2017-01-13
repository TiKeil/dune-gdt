// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include <config.h>

//#if HAVE_DUNE_PYBINXI

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/xt/common/string.hh>
#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/spaces/fv.hh>
#include <dune/gdt/spaces/rt.hh>

#include <dune/gdt/spaces.pbh>

namespace py = pybind11;
using namespace pybind11::literals;


template <class SP>
void bind_space(py::module& m, const std::string& grid_id, const std::string& space_id, const std::string& backend)
{
  using namespace Dune::XT;
  using Common::to_string;
  using Common::to_lower;
  typedef typename SP::type S; //            without the size_t(...) we get linker errors on module import
  const std::string r = to_string(size_t(S::dimRange));
  const std::string rC = to_string(size_t(S::dimRangeCols));
  const std::string p = to_string(int(S::polOrder));
  const std::string space_suffix = grid_id + "_to_" + r + "x" + rC + "__p" + p + backend;
  auto s = Dune::GDT::bind_space<S>(m, space_id + "Space__" + space_suffix);

  typedef typename Grid::extract_grid<typename S::GridViewType>::type G;
  m.def(std::string("make_" + to_lower(space_id) + "_space__" + space_suffix).c_str(),
        [](Grid::GridProvider<G>& grid_provider, const int level = 0) { return SP::create(grid_provider, level); },
        "grid_provider"_a,
        "level"_a = 0);
} // ... bind_space(...)


template <class G>
void addbind_for_Grid(py::module& m, const std::string& grid_id)
{
  // FV
  bind_space<Dune::GDT::
                 FvSpaceProvider<G, Dune::XT::Grid::Layers::leaf, Dune::GDT::ChooseSpaceBackend::gdt, double, 1, 1>>(
      m, grid_id + "_leaf", "Fv", "");
  bind_space<Dune::GDT::FvSpaceProvider<G,
                                        Dune::XT::Grid::Layers::leaf,
                                        Dune::GDT::ChooseSpaceBackend::gdt,
                                        double,
                                        G::dimension,
                                        1>>(m, grid_id + "_leaf", "Fv", "");
// CG
#if HAVE_DUNE_FEM
  bind_space<Dune::GDT::
                 CgSpaceProvider<G, Dune::XT::Grid::Layers::leaf, Dune::GDT::ChooseSpaceBackend::fem, 1, double, 1, 1>>(
      m, grid_id + "_leaf", "Cg", "__fem");
  bind_space<Dune::GDT::CgSpaceProvider<G,
                                        Dune::XT::Grid::Layers::level,
                                        Dune::GDT::ChooseSpaceBackend::fem,
                                        1,
                                        double,
                                        1,
                                        1>>(m, grid_id + "_level", "Cg", "__fem");
#endif
#if HAVE_DUNE_PDELAB
  bind_space<Dune::GDT::CgSpaceProvider<G,
                                        Dune::XT::Grid::Layers::leaf,
                                        Dune::GDT::ChooseSpaceBackend::pdelab,
                                        1,
                                        double,
                                        1,
                                        1>>(m, grid_id + "_leaf", "Cg", "__pdelab");
  bind_space<Dune::GDT::CgSpaceProvider<G,
                                        Dune::XT::Grid::Layers::level,
                                        Dune::GDT::ChooseSpaceBackend::pdelab,
                                        1,
                                        double,
                                        1,
                                        1>>(m, grid_id + "_level", "Cg", "__pdelab");
#endif
// DG
#if HAVE_DUNE_FEM
  bind_space<Dune::GDT::
                 DgSpaceProvider<G, Dune::XT::Grid::Layers::leaf, Dune::GDT::ChooseSpaceBackend::fem, 1, double, 1, 1>>(
      m, grid_id + "_leaf", "Dg", "__fem");
  bind_space<Dune::GDT::DgSpaceProvider<G,
                                        Dune::XT::Grid::Layers::level,
                                        Dune::GDT::ChooseSpaceBackend::fem,
                                        1,
                                        double,
                                        1,
                                        1>>(m, grid_id + "_level", "Dg", "__fem");
#endif
//#if HAVE_DUNE_PDELAB // these need to be guarded against simplicial grids
//  bind_space<Dune::GDT::DgSpaceProvider<G,
//                                        Dune::XT::Grid::Layers::leaf,
//                                        Dune::GDT::ChooseSpaceBackend::pdelab,
//                                        1,
//                                        double,
//                                        1,
//                                        1>>(m, grid_id + "_leaf", "Dg", "__pdelab");
//  bind_space<Dune::GDT::DgSpaceProvider<G,
//                                        Dune::XT::Grid::Layers::level,
//                                        Dune::GDT::ChooseSpaceBackend::pdelab,
//                                        1,
//                                        double,
//                                        1,
//                                        1>>(m, grid_id + "_level", "Dg", "__pdelab");
//#endif
// RT
#if HAVE_DUNE_PDELAB
  bind_space<Dune::GDT::RtSpaceProvider<G,
                                        Dune::XT::Grid::Layers::leaf,
                                        Dune::GDT::ChooseSpaceBackend::pdelab,
                                        0,
                                        double,
                                        G::dimension,
                                        1>>(m, grid_id + "_leaf", "Rt", "__pdelab");
  bind_space<Dune::GDT::RtSpaceProvider<G,
                                        Dune::XT::Grid::Layers::level,
                                        Dune::GDT::ChooseSpaceBackend::pdelab,
                                        0,
                                        double,
                                        G::dimension,
                                        1>>(m, grid_id + "_level", "Rt", "__pdelab");
#endif
} // ... addbind_for_Grid(...)


PYBIND11_PLUGIN(gdt)
{
  py::module m("gdt", "dune-gdt");

  py::module::import("common");
  py::module::import("grid");
  py::module::import("functions");
  py::module::import("la");

  m.def("init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = args.size();
          char** argv = Dune::XT::Common::vector_to_main_args(args);
#if HAVE_DUNE_FEM
          Dune::Fem::MPIManager::initialize(argc, argv);
#else
          Dune::MPIHelper::instance(argc, argv);
#endif
        },
        "args"_a = std::vector<std::string>());

  typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> G;
  const std::string grid_id = "2d_simplex_aluconform";

  addbind_for_Grid<G>(m, grid_id);

  return m.ptr();
}


//#endif // HAVE_DUNE_PYBINXI
