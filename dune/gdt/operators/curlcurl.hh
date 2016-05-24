// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_CURLCURL_HH
#define DUNE_GDT_OPERATORS_CURLCURL_HH

#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/local/integrands/curlcurl.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/integrals.hh>
//#include <dune/gdt/assembler/local/codim0.hh>

//new
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/spaces/interface.hh>

#include <dune/gdt/assembler/system.hh>

#include "base.hh"

/*
#ifndef DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_OPERATORS_WEIGHTED_L2_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/spaces/interface.hh>
*/

namespace Dune {
namespace GDT {


// ////////////////////////// //
// CurlCurlLocalizableProduct //
// ////////////////////////// //

template <class FunctionType, class GridView, class Range, class Source = Range,
          class Field                                                         = typename Range::RangeFieldType>
class CurlCurlLocalizableProduct : public LocalizableProductBase<GridView, Range, Source, Field>
{
  typedef LocalizableProductBase<GridView, Range, Source, Field> BaseType;
  typedef LocalVolumeIntegralOperator<LocalCurlCurlIntegrand<FunctionType>> CurlCurlLocalizableOperatorType;

public:
  template <class... Args>
  CurlCurlLocalizableProduct(const FunctionType& function, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_curlcurl_operator_(function)
  {
    this->add(local_curlcurl_operator_);
  }

  template <class... Args>
  WeightedL2LocalizableProduct(const size_t over_integrate, const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(over_integrate, weight)
  {
    this->add(local_weighted_l2_operator_);
  }

private:
  const LocalWeightedL2OperatorType local_weighted_l2_operator_;
}; // class WeightedL2LocalizableProduct


// //////////////////////////////////// //
// make_weighted_l2_localizable_product //
// //////////////////////////////////// //

/**
 * \sa WeightedL2LocalizableProduct
 */
template <class WeightFunctionType, class GridViewType, class RangeType, class SourceType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value
                            && Stuff::is_localizable_function<RangeType>::value
                            && Stuff::is_localizable_function<SourceType>::value,
                        std::unique_ptr<WeightedL2LocalizableProduct<WeightFunctionType, GridViewType, RangeType,
                                                                     SourceType>>>::type
make_weighted_l2_localizable_product(const WeightFunctionType& weight, const GridViewType& grid_view,
                                     const RangeType& range, const SourceType& source, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2LocalizableProduct<WeightFunctionType, GridViewType, RangeType, SourceType>>(
      over_integrate, weight, grid_view, range, source);
}


// //////////////////////// //
// WeightedL2MatrixOperator //
// //////////////////////// //

template <class WeightFunctionType, class RangeSpace,
          class Matrix   = typename Stuff::LA::Container<typename RangeSpace::RangeFieldType>::MatrixType,
          class GridView = typename RangeSpace::GridViewType, class SourceSpace = RangeSpace,
          class Field = typename RangeSpace::RangeFieldType>
class WeightedL2MatrixOperator
    : public MatrixOperatorBase<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume>
{
  typedef MatrixOperatorBase<Matrix, RangeSpace, GridView, SourceSpace, Field, ChoosePattern::volume> BaseType;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<WeightFunctionType>> LocalWeightedL2OperatorType;

public:
  template <class... Args>
  explicit WeightedL2MatrixOperator(const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(weight)
  {
    this->add(local_weighted_l2_operator_);
  }

  template <class... Args>
  explicit WeightedL2MatrixOperator(const size_t over_integrate, const WeightFunctionType& weight, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_weighted_l2_operator_(over_integrate, weight)
  {
    this->add(local_weighted_l2_operator_);
  }

private:
  const LocalWeightedL2OperatorType local_weighted_l2_operator_;
}; // class WeightedL2MatrixOperator


// //////////////////////////////// //
// make_weighted_l2_matrix_operator //
// //////////////////////////////// //

// without matrix

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically,
 *        source and range space are given by space, grid_view of the space is used).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value
                            && Stuff::is_localizable_function<WeightFunctionType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically,
 *        source and range space are given by space).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, space, grid_view);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class SpaceType, class GridViewType>
typename std::
    enable_if<Stuff::LA::is_matrix<MatrixType>::value && Stuff::is_localizable_function<WeightFunctionType>::value
                  && is_space<SpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
              std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>>::type
    make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const SpaceType& space,
                                     const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>(
      over_integrate, weight, space, grid_view);
}

/**
 * \brief Creates a weighted L2 matrix operator (MatrixType has to be supllied, a matrix is created automatically).
 * \note  MatrixType has to be supplied, i.e., use like
\code
auto op = make_weighted_l2_matrix_operator< MatrixType >(weight, range_space, source_space, grid_view);
\endcode
 */
template <class MatrixType, class WeightFunctionType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::enable_if<Stuff::LA::is_matrix<MatrixType>::value
                            && Stuff::is_localizable_function<WeightFunctionType>::value
                            && is_space<RangeSpaceType>::value && is_space<SourceSpaceType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType,
                                                                 GridViewType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, const RangeSpaceType& range_space,
                                 const SourceSpaceType& source_space, const GridViewType& grid_view,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType,
                                                   RangeSpaceType,
                                                   MatrixType,
                                                   GridViewType,
                                                   SourceSpaceType>>(
      over_integrate, weight, range_space, source_space, grid_view);
}

// with matrix

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space, grid_view of the space is
 *        used).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::LA::is_matrix<MatrixType>::value && is_space<SpaceType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix, const SpaceType& space,
                                 const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType>>(
      over_integrate, weight, matrix, space);
}

/**
 * \brief Creates a weighted L2 matrix operator (source and range space are given by space).
 */
template <class WeightFunctionType, class MatrixType, class SpaceType, class GridViewType>
typename std::
    enable_if<Stuff::is_localizable_function<WeightFunctionType>::value && Stuff::LA::is_matrix<MatrixType>::value
                  && is_space<SpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
              std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>>::type
    make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix, const SpaceType& space,
                                     const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType, SpaceType, MatrixType, GridViewType>>(
      over_integrate, weight, matrix, space, grid_view);
}

/**
 * \brief Creates a weighted L2 matrix operator.
 */
template <class WeightFunctionType, class MatrixType, class RangeSpaceType, class SourceSpaceType, class GridViewType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::LA::is_matrix<MatrixType>::value && is_space<RangeSpaceType>::value
                            && is_space<SourceSpaceType>::value && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2MatrixOperator<WeightFunctionType, RangeSpaceType, MatrixType,
                                                                 GridViewType, SourceSpaceType>>>::type
make_weighted_l2_matrix_operator(const WeightFunctionType& weight, MatrixType& matrix,
                                 const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                                 const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<WeightedL2MatrixOperator<WeightFunctionType,
                                                   RangeSpaceType,
                                                   MatrixType,
                                                   GridViewType,
                                                   SourceSpaceType>>(
      over_integrate, weight, matrix, range_space, source_space, grid_view);
}


// ////////////////// //
// WeightedL2Operator //
// ////////////////// //

// forward, needed for the traits
template <class WeightFunctionType, class GridView, class Field = typename WeightFunctionType::RangeFieldType>
class WeightedL2Operator;


namespace internal {


template <class WeightFunctionType, class GridViewType, class Field>
class WeightedL2OperatorTraits
{
public:
  typedef WeightedL2Operator<WeightFunctionType, GridViewType, Field> derived_type;
  typedef Field FieldType;
};


} // namespace internal


template <class WeightFunctionType, class GridViewType, class Field>
class WeightedL2Operator
    : public OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridViewType, Field>>
{
  typedef OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridViewType, Field>> BaseType;

public:
  using typename BaseType::FieldType;

  WeightedL2Operator(const WeightFunctionType& weight, GridViewType grid_view, const size_t over_integrate = 0)
    : weight_(weight)
    , grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {
  }

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range) const
  {
    typedef typename Stuff::LA::Container<typename VectorType::ScalarType, VectorType::sparse_matrix_type>::MatrixType
        MatrixType;
    auto op = make_weighted_l2_matrix_operator<MatrixType>(
        weight_, source.space(), range.space(), grid_view_, over_integrate_);
    op->apply(source, range);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& range,
                   const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& source) const
  {
    auto product = make_weighted_l2_localizable_product(weight_, grid_view_, range, source, over_integrate_);
    return product->apply2();
  }

  using BaseType::apply_inverse;

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/, SourceType& /*source*/,
                     const Stuff::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "yet");
    return {"depends_on_the_vector_type_of_the_discrete_function"};
  }

  Stuff::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

private:
  const WeightFunctionType& weight_;
  GridViewType grid_view_;
  const size_t over_integrate_;
}; // class WeightedL2Operator


// ///////////////////////// //
// make_weighted_l2_operator //
// ///////////////////////// //

template <class GridViewType, class WeightFunctionType>
typename std::enable_if<Stuff::is_localizable_function<WeightFunctionType>::value
                            && Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<WeightedL2Operator<WeightFunctionType, GridViewType,
                                                           typename WeightFunctionType::RangeFieldType>>>::type
make_weighted_l2_operator(const GridViewType& grid_view, const WeightFunctionType& weight,
                          const size_t over_integrate = 0)
{
  return DSC::
      make_unique<WeightedL2Operator<WeightFunctionType, GridViewType, typename WeightFunctionType::RangeFieldType>>(
          weight, grid_view, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_WEIGHTED_L2_HH


namespace Dune {
namespace GDT {


//forward
template< class FunctionType, class MatrixImp, class SourceSpaceImp,
          class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp = typename SourceSpaceImp::GridViewType >
class CurlCurl;



namespace internal {

/** \brief Traits for the CurlCurl-operator
 *
 *\sa CurlCurl
 */
template< class FunctionType, class MatrixImp, class SourceSpaceImp,
          class RangeSpaceImp, class GridViewImp >
class CurlCurlTraits
{
  static_assert(Stuff::is_localizable_function< FunctionType >::value,
                "FunctionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::LA::is_matrix< MatrixImp >::value,
                "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(is_space< SourceSpaceImp >::value, "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(is_space< RangeSpaceImp >::value, "RangeSpaceImp has to be derived from SpaceInterface!");
public:
  typedef CurlCurl< FunctionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp >
          derived_type;
  typedef MatrixImp      MatrixType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef RangeSpaceImp  RangeSpaceType;
  typedef GridViewImp    GridViewType;
}; //class CurlCurlTraits

} //namespace internal


/** \brief Implements a (global) curlcurl operator
 * \note only scalar parameters are supported for sure
 *
 * \tparam FunctionType Type of the parameter function
 * \tparam MatrixImp Type for the system matrix everything is assembled in
 * \tparam SourceSpaceImp Type of the ansatz space
 * \tparam RangeSpaceImp Type of the test space
 * \tparam GridViewImp Type of the grid
 */
template< class FunctionType, class MatrixImp, class SourceSpaceImp,
          class RangeSpaceImp, class GridViewImp >
class CurlCurl
  : Stuff::Common::StorageProvider< MatrixImp >
  , public Operators::MatrixBased< internal::CurlCurlTraits< FunctionType, MatrixImp,
                                                             SourceSpaceImp, RangeSpaceImp, GridViewImp > >
  , public SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp >
{
  typedef Stuff::Common::StorageProvider< MatrixImp >                                        StorageProvider;
  typedef SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp >                      AssemblerBaseType;
  typedef Operators::MatrixBased< internal::CurlCurlTraits< FunctionType, MatrixImp, SourceSpaceImp,
                                                            RangeSpaceImp, GridViewImp > >   OperatorBaseType;
  typedef LocalOperator::Codim0Integral< LocalEvaluation::CurlCurl< FunctionType > >         LocalCurlOperatorType;
  typedef LocalAssembler::Codim0Matrix< LocalCurlOperatorType >                              LocalCurlAssemblerType;
  public:
    typedef internal::CurlCurlTraits< FunctionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp>
        Traits;

  typedef typename Traits::MatrixType      MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType  RangeSpaceType;
  typedef typename Traits::GridViewType    GridViewType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    return range_space.compute_face_and_volume_pattern(grid_view, source_space);
  }

  CurlCurl(const FunctionType& mu,
           MatrixType& mtrx,
           const SourceSpaceType& src_space,
           const RangeSpaceType& rng_space,
           const GridViewType& grid_view)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_space, rng_space, grid_view)
    , AssemblerBaseType(rng_space, grid_view, src_space)
    , mu_(mu)
    , local_curl_operator_(mu_)
    , local_curl_assembler_(local_curl_operator_)
    , assembled_(false)
  {
    setup();
  }

  CurlCurl(const FunctionType& mu,
           const SourceSpaceType& src_space,
           const RangeSpaceType& rng_space,
           const GridViewType& grid_view)
    : StorageProvider(new MatrixType(rng_space.mapper().size(),
                                     src_space.mapper().size(),
                                     pattern(rng_space, src_space, grid_view)))
    , OperatorBaseType(this->storage_access(), src_space, rng_space, grid_view)
    , AssemblerBaseType(rng_space, grid_view, src_space)
    , mu_(mu)
    , local_curl_operator_(mu_)
    , local_curl_assembler_(local_curl_operator_)
    , assembled_(false)
  {
    setup();
  }

  CurlCurl(const FunctionType& mu,
           MatrixType& mtrx,
           const SourceSpaceType& src_space,
           const RangeSpaceType& rng_space)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_space, rng_space)
    , AssemblerBaseType(rng_space, src_space)
    , mu_(mu)
    , local_curl_operator_(mu_)
    , local_curl_assembler_(local_curl_operator_)
    , assembled_(false)
  {
    setup();
  }

  CurlCurl(const FunctionType& mu,
           const SourceSpaceType& src_space,
           const RangeSpaceType& rng_space)
    : StorageProvider(new MatrixType(rng_space.mapper().size(),
                                     src_space.mapper().size(),
                                     pattern(rng_space, src_space)))
    , OperatorBaseType(this->storage_access(), src_space, rng_space)
    , AssemblerBaseType(rng_space, src_space)
    , mu_(mu)
    , local_curl_operator_(mu_)
    , local_curl_assembler_(local_curl_operator_)
    , assembled_(false)
  {
    setup();
  }

  CurlCurl(const FunctionType& mu,
           MatrixType& mtrx,
           const SourceSpaceType& src_space)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_space)
    , AssemblerBaseType(src_space)
    , mu_(mu)
    , local_curl_operator_(mu_)
    , local_curl_assembler_(local_curl_operator_)
    , assembled_(false)
  {
    setup();
  }

  CurlCurl(const FunctionType& mu,
           const SourceSpaceType& src_space)
    : StorageProvider(new MatrixType(src_space.mapper().size(), src_space.mapper().size(),
                                     pattern(src_space)))
    , OperatorBaseType(this->storage_access(), src_space)
    , AssemblerBaseType(src_space)
    , mu_(mu)
    , local_curl_operator_(mu_)
    , local_curl_assembler_(local_curl_operator_)
    , assembled_(false)
  {
    setup();
  }

  virtual ~CurlCurl() {}

  virtual void assemble() override final
  {
    if(!assembled_) {
      AssemblerBaseType::assemble(true);
      assembled_ = true;
    }
  } // ... assemble(...)


private:
  void setup()
  {
    this->add(local_curl_assembler_, this->matrix());
  } //... setup()

  const FunctionType& mu_;
  const LocalCurlOperatorType local_curl_operator_;
  const LocalCurlAssemblerType local_curl_assembler_;
  bool assembled_;
};  //class CurlCurl


} //namespace GDT
} //namespace Dune




#endif // DUNE_GDT_OPERATORS_CURLCURL_HH
