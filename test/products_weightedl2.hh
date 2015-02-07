// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PRODUCTS_WEIGHTEDL2_HH
#define DUNE_GDT_TEST_PRODUCTS_WEIGHTEDL2_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/la/container/common.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/weightedl2.hh>

#include "products.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class SpaceType>
struct WeightedL2ProductBase : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const size_t dimRange = SpaceType::dimRange;
  typedef Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange> FunctionType;

  WeightedL2ProductBase()
    : grid_(GridProviderType(0.0, 1.0, boost::numeric_cast<size_t>(dsc_grid_elements())).grid_ptr())
    , space_(Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(*grid_))
    , one_("x", "1.0", 0)
  {
  }

  virtual RangeFieldType compute(const FunctionType& function) const = 0;

  void constant_arguments() const
  {
    check(compute(one_), 1.0);
  } // ... constant_arguments(...)

  void linear_arguments() const
  {
    const FunctionType linear("x", "x[0] - 1.0", 1);
    check(compute(linear), 1.0 / 3.0);
  }

  void quadratic_arguments() const
  {
    const FunctionType quadratic("x", "x[0]*x[0]", 2);
    check(compute(quadratic), 1.0 / 5.0);
  }

  void check(const RangeFieldType& result, const RangeFieldType& expected, const RangeFieldType epsilon = 1e-14) const
  {
    const auto error = std::abs(expected - result);
    EXPECT_LE(error, epsilon) << "result:   " << result << "\n"
                              << "expected: " << expected << "\n"
                              << "difference: " << std::scientific << error;
  } // ... check(...)

  std::shared_ptr<GridType> grid_;
  const SpaceType space_;
  const FunctionType one_;
}; // struct WeightedL2ProductBase


template <class SpaceType>
struct WeightedL2LocalizableProduct : public WeightedL2ProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> BaseType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  virtual RangeFieldType compute(const FunctionType& function) const override
  {
    return Products::WeightedL2Localizable<GridViewType, FunctionType, FunctionType, FunctionType>(
               this->space_.grid_view(), function, function, this->one_)
        .apply2();
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::WeightedL2Localizable<GridViewType, FunctionType, FunctionType, FunctionType> ProductType;
    ProductType product(this->space_.grid_view(), this->one_, this->one_, this->one_);
    LocalizableProductBase<SpaceType, ProductType>::fulfills_interface(product);
  }
}; // struct WeightedL2LocalizableProduct


template <class SpaceType>
struct WeightedL2AssemblableProduct : public WeightedL2ProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename Dune::Stuff::LA::CommonDenseVector<RangeFieldType> VectorType;
  typedef typename Dune::Stuff::LA::CommonDenseMatrix<RangeFieldType> MatrixType;
  typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
  typedef Dune::GDT::Operators::Projection<GridViewType> ProjectionOperatorType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    // create the product
    typedef Products::WeightedL2Assemblable<MatrixType, FunctionType, SpaceType, GridViewType, SpaceType> Product;
    Product product(this->space_, this->one_);
    product.assemble(false);
    // project the function
    DiscreteFunctionType discrete_function(this->space_);
    ProjectionOperatorType(this->space_.grid_view()).apply(function, discrete_function);
    // compute the product
    const auto result = product.apply2(discrete_function, discrete_function);
    Product product_tbb(this->space_, this->one_);
    product_tbb.assemble(true);
    const auto result_tbb = product_tbb.apply2(discrete_function, discrete_function);
    EXPECT_DOUBLE_EQ(result_tbb, result);
    return result;
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::WeightedL2Assemblable<MatrixType, FunctionType, SpaceType, GridViewType, SpaceType> ProductType;
    ProductType product(this->space_, this->one_);
    AssemblableProductBase<SpaceType, ProductType, VectorType>::fulfills_interface(product);
  }
}; // struct WeightedL2AssemblableProduct


template <class SpaceType>
struct WeightedL2Product : public WeightedL2ProductBase<SpaceType>
{
  typedef WeightedL2ProductBase<SpaceType> BaseType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::GridViewType GridViewType;

  virtual RangeFieldType compute(const FunctionType& function) const override final
  {
    Products::WeightedL2<GridViewType, FunctionType> product(this->space_.grid_view(), this->one_);
    return product.apply2(function, function);
  } // ... compute(...)

  void fulfills_interface() const
  {
    typedef Products::WeightedL2<GridViewType, FunctionType> ProductType;
    ProductBase<SpaceType, ProductType>::fulfills_interface(ProductType(this->space_.grid_view(), this->one_));
  }
}; // struct WeightedL2Product


#endif // DUNE_GDT_TEST_PRODUCTS_WEIGHTEDL2_HH
