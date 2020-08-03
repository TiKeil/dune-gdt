// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/interpolations/default.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>
#include <python/dune/gdt/operators/interfaces.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class VectorTag, class GV, size_t r = 1>
class default_interpolation
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;
  using E = XT::Grid::extract_entity_t<GV>;
  using GF = XT::Functions::GridFunction<E, r>;
  using DF = DiscreteFunction<V, GV, r>;
  using S = typename DF::SpaceType;

public:
  static void bind(pybind11::module& m, const std::string& function_id = "default_interpolation")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        function_id.c_str(),
        [](GF source, DF& target) { GDT::default_interpolation(source, target); },
        "source"_a,
        "target"_a,
        py::call_guard<py::gil_scoped_release>());
    if (std::is_same<VectorTag, XT::LA::bindings::Istl>::value)
      m.def(
          function_id.c_str(),
          [](GF source, const S& target_space, const VectorTag&) {
            return GDT::default_interpolation<V>(source, target_space);
          },
          "source"_a,
          "target"_a,
          "la_backend"_a = XT::LA::bindings::Istl(),
          py::call_guard<py::gil_scoped_release>());
    else
      m.def(
          function_id.c_str(),
          [](GF source, const S& target_space, const VectorTag&) {
            return GDT::default_interpolation<V>(source, target_space);
          },
          "source"_a,
          "target"_a,
          "la_backend"_a,
          py::call_guard<py::gil_scoped_release>());
  } // ... bind(...)
}; // class default_interpolation


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class VT, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct default_interpolation_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::default_interpolation;

    default_interpolation<V, VT, GV>::bind(m);
    if (d > 1)
      default_interpolation<V, VT, GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    default_interpolation_for_all_grids<V, VT, typename GridTypes::tail_type>::bind(m);
  }
};

template <class V, class VT>
struct default_interpolation_for_all_grids<V, VT, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_interpolations_default, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._discretefunction_discretefunction");
  py::module::import("dune.gdt._spaces_interface");

#if 0
  // bindings for all but dune-istl disabled for the moment
  default_interpolation_for_all_grids<LA::CommonDenseVector<double>,
                                      LA::bindings::Common,
                                      XT::Grid::AvailableGridTypes>::bind(m);
#  if HAVE_EIGEN
  default_interpolation_for_all_grids<LA::EigenDenseVector<double>, LA::bindings::Eigen, XT::Grid::AvailableGridTypes>::
      bind(m);
#  endif
#endif // 0
  default_interpolation_for_all_grids<LA::IstlDenseVector<double>, LA::bindings::Istl, XT::Grid::AvailableGridTypes>::
      bind(m);
}