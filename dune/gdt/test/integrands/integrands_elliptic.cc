// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/matrix.hh>

#include <dune/xt/la/container/eye-matrix.hh>

#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/local/integrands/elliptic.hh>

#include <dune/gdt/test/integrands/integrands.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct EllipticIntegrandTest : public IntegrandTest<G>
{
  using BaseType = IntegrandTest<G>;
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::GV;
  using typename BaseType::VectorJacobianType;
  using ScalarIntegrandType = LocalEllipticIntegrand<E, 1>;
  using VectorIntegrandType = LocalEllipticIntegrand<E, d>;

  virtual void SetUp() override
  {
    BaseType::SetUp();
    diffusion_factor_ = std::make_shared<XT::Functions::GenericGridFunction<E, 1>>(
        2, [](const E&) {}, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    diffusion_tensor_ = std::make_shared<XT::Functions::GenericGridFunction<E, 2, 2>>(
        1,
        [](const E&) {},
        [](const DomainType& x, const XT::Common::Parameter&) {
          return VectorJacobianType{{x[0], x[1]}, {1., 2.}};
        });
  }

  virtual void is_constructable() override final
  {
    ScalarIntegrandType scalar_integrand1;
    ScalarIntegrandType scalar_integrand2(1., XT::LA::eye_matrix<FieldMatrix<D, d, d>>(d, d));
    const XT::Functions::GenericFunction<d, 1> scalar_function(
        2, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    const XT::Functions::GenericFunction<d, 2, 2> matrix_function(
        1, [](const DomainType& x, const XT::Common::Parameter&) {
          return VectorJacobianType{{x[0], x[1]}, {1., 2.}};
        });
    ScalarIntegrandType scalar_integrand3(scalar_function, matrix_function);
    ScalarIntegrandType scalar_integrand4(*diffusion_factor_, *diffusion_tensor_);
    DUNE_UNUSED_PARAMETER(scalar_integrand1);
    DUNE_UNUSED_PARAMETER(scalar_integrand2);
    DUNE_UNUSED_PARAMETER(scalar_integrand3);
    DUNE_UNUSED_PARAMETER(scalar_integrand4);
    VectorIntegrandType vector_integrand1;
    VectorIntegrandType vector_integrand2(1., XT::LA::eye_matrix<FieldMatrix<D, d, d>>(d, d));
    VectorIntegrandType vector_integrand3(scalar_function, matrix_function);
    VectorIntegrandType vector_integrand4(*diffusion_factor_, *diffusion_tensor_);
    DUNE_UNUSED_PARAMETER(vector_integrand1);
    DUNE_UNUSED_PARAMETER(vector_integrand2);
    DUNE_UNUSED_PARAMETER(vector_integrand3);
    DUNE_UNUSED_PARAMETER(vector_integrand4);
  }

  virtual void evaluates_correctly_for_scalar_bases()
  {
    ScalarIntegrandType scalar_integrand(*diffusion_factor_, *diffusion_tensor_);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    scalar_integrand.bind(element);
    const auto integrand_order = scalar_integrand.order(*scalar_test_, *scalar_ansatz_);
    EXPECT_EQ(8, integrand_order);
    DynamicMatrix<D> result(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      scalar_integrand.evaluate(*scalar_test_, *scalar_ansatz_, x, result);
      DynamicMatrix<D> expected_result{
          {1, 2 * (x[0] * x[1] + std::pow(x[0], 2))},
          {x[0] * std::pow(x[1], 3) + 3 * x[0] * std::pow(x[1], 2),
           3 * std::pow(x[0], 2) * std::pow(x[1], 4)
               + 6 * (std::pow(x[0], 2) * std::pow(x[1], 3) + std::pow(x[0], 3) * std::pow(x[1], 2))}};
      expected_result *= x[0] * x[1];
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < 2; ++jj)
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
    }
  }

  virtual void evaluates_correctly_for_vector_bases()
  {
    VectorIntegrandType integrand(*diffusion_factor_, *diffusion_tensor_);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    integrand.bind(element);
    const auto integrand_order = integrand.order(*vector_test_, *vector_ansatz_);
    EXPECT_EQ(5, integrand_order);
    DynamicMatrix<D> result(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      integrand.evaluate(*vector_test_, *vector_ansatz_, x, result);
      DynamicMatrix<D> expected_result{
          {0., 0.},
          {x[0] * x[1] + x[0] + x[1] + 2,
           2 * (std::pow(x[0], 2) * x[1] + std::pow(x[1], 2) + std::pow(x[0], 2) + 2 * x[1])}};
      expected_result *= x[0] * x[1];
      for (size_t ii = 1; ii < 2; ++ii)
        for (size_t jj = 0; jj < 2; ++jj) {
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
        }
    }
  }

  using BaseType::grid_provider_;
  using BaseType::scalar_ansatz_;
  using BaseType::scalar_test_;
  using BaseType::vector_ansatz_;
  using BaseType::vector_test_;
  std::shared_ptr<XT::Functions::GenericGridFunction<E, 1>> diffusion_factor_;
  std::shared_ptr<XT::Functions::GenericGridFunction<E, 2, 2>> diffusion_tensor_;
}; // struct EllipticIntegrandTest


} // namespace Test
} // namespace GDT
} // namespace Dune


template <class G>
using EllipticIntegrandTest = Dune::GDT::Test::EllipticIntegrandTest<G>;
TYPED_TEST_CASE(EllipticIntegrandTest, Grids2D);

TYPED_TEST(EllipticIntegrandTest, is_constructable)
{
  this->is_constructable();
}
TYPED_TEST(EllipticIntegrandTest, evaluates_correctly_for_scalar_bases)
{
  this->evaluates_correctly_for_scalar_bases();
}

TYPED_TEST(EllipticIntegrandTest, evaluates_correctly_for_vector_bases)
{
  this->evaluates_correctly_for_vector_bases();
}
