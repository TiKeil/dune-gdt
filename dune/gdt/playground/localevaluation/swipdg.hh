// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_LOCALEVALUATION_SWIPDG_HH
#define DUNE_GDT_PLAYGROUND_LOCALEVALUATION_SWIPDG_HH

#include <dune/stuff/common/fmatrix.hh>

#include <dune/gdt/localevaluation/swipdg.hh>

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace SWIPDG {


template <class DiffusionFactorType, class DiffusionTensorType>
class Inner
    : public LocalEvaluation::Codim1Interface<internal::InnerTraits<DiffusionFactorType, DiffusionTensorType>, 4>
{
public:
  typedef internal::InnerTraits<DiffusionFactorType, DiffusionTensorType> Traits;

  Inner(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& inducingFunction,
        const double beta = internal::default_beta(DiffusionTensorType::dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(inducingFunction)
    , beta_(beta)
  {
  }

  template <class EntityType>
  class LocalfunctionTuple
  {
    typedef typename DiffusionFactorType::LocalfunctionType LocalDiffusionFactorType;
    typedef typename DiffusionTensorType::LocalfunctionType LocalDiffusionTensorType;

  public:
    typedef std::tuple<std::shared_ptr<LocalDiffusionFactorType>, std::shared_ptr<LocalDiffusionTensorType>> Type;
  };

  template <class EntityType>
  typename LocalfunctionTuple<EntityType>::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity), diffusion_tensor_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class E, class N, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  size_t order(const typename LocalfunctionTuple<E>::Type& localFunctionsEntity,
               const typename LocalfunctionTuple<N>::Type& localFunctionsNeighbor,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBaseEntity,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBaseEntity,
               const Stuff::LocalfunctionSetInterface<N, D, d, R, rT, rCT>& testBaseNeighbor,
               const Stuff::LocalfunctionSetInterface<N, D, d, R, rA, rCA>& ansatzBaseNeighbor) const
  {
    const auto local_diffusion_factor_entity   = std::get<0>(localFunctionsEntity);
    const auto local_diffusion_tensor_entity   = std::get<1>(localFunctionsEntity);
    const auto local_diffusion_factor_neighbor = std::get<0>(localFunctionsNeighbor);
    const auto local_diffusion_tensor_neighbor = std::get<1>(localFunctionsNeighbor);
    return order(*local_diffusion_factor_entity,
                 *local_diffusion_tensor_entity,
                 *local_diffusion_factor_neighbor,
                 *local_diffusion_tensor_neighbor,
                 testBaseEntity,
                 ansatzBaseEntity,
                 testBaseNeighbor,
                 ansatzBaseNeighbor);
  }

  template <class E, class N, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rT, int rCT, int rA,
            int rCA>
  size_t order(const Stuff::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& localDiffusionFactorEntity,
               const Stuff::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& localDiffusionTensorEntity,
               const Stuff::LocalfunctionInterface<N, D, d, R, rDF, rCDF>& localDiffusionFactorNeighbor,
               const Stuff::LocalfunctionInterface<N, D, d, R, rDT, rCDT>& localDiffusionTensorNeighbor,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBaseEntity,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBaseEntity,
               const Stuff::LocalfunctionSetInterface<N, D, d, R, rT, rCT>& testBaseNeighbor,
               const Stuff::LocalfunctionSetInterface<N, D, d, R, rA, rCA>& ansatzBaseNeighbor) const
  {
    return std::max(localDiffusionFactorEntity.order(), localDiffusionFactorNeighbor.order())
           + std::max(localDiffusionTensorEntity.order(), localDiffusionTensorNeighbor.order())
           + std::max(testBaseEntity.order(), testBaseNeighbor.order())
           + std::max(ansatzBaseEntity.order(), ansatzBaseNeighbor.order());
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class E, class N, class IntersectionType, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const typename LocalfunctionTuple<E>::Type& localFunctionsEntity,
                const typename LocalfunctionTuple<N>::Type& localFunctionsNeighbor,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBaseEntity,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface<N, D, d, R, rT, rCT>& testBaseNeighbor,
                const Stuff::LocalfunctionSetInterface<N, D, d, R, rA, rCA>& ansatzBaseNeighbor,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
                Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    const auto local_diffusion_factor_entity   = std::get<0>(localFunctionsEntity);
    const auto local_diffusion_tensor_entity   = std::get<1>(localFunctionsEntity);
    const auto local_diffusion_factor_neighbor = std::get<0>(localFunctionsNeighbor);
    const auto local_diffusion_tensor_neighbor = std::get<1>(localFunctionsNeighbor);
    evaluate(*local_diffusion_factor_entity,
             *local_diffusion_tensor_entity,
             *local_diffusion_factor_neighbor,
             *local_diffusion_tensor_neighbor,
             testBaseEntity,
             ansatzBaseEntity,
             testBaseNeighbor,
             ansatzBaseNeighbor,
             intersection,
             localPoint,
             entityEntityRet,
             neighborNeighborRet,
             entityNeighborRet,
             neighborEntityRet);
  }

  template <class E, class N, class IntersectionType, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT,
            int rT, int rCT, int rA, int rCA>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& /*localDiffusionFactorEntity*/,
                const Stuff::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& /*localDiffusionTensorEntity*/,
                const Stuff::LocalfunctionInterface<N, D, d, R, rDF, rCDF>& /*localDiffusionFactorNeighbor*/,
                const Stuff::LocalfunctionInterface<N, D, d, R, rDT, rCDT>& /*localDiffusionTensorNeighbor*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& /*ansatzBaseEntity*/,
                const Stuff::LocalfunctionSetInterface<N, D, d, R, rT, rCT>& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface<N, D, d, R, rA, rCA>& /*ansatzBaseNeighbor*/,
                const IntersectionType& /*intersection*/, const Dune::FieldVector<D, d - 1>& /*localPoint*/,
                Dune::DynamicMatrix<R>& /*entityEntityRet*/, Dune::DynamicMatrix<R>& /*neighborNeighborRet*/,
                Dune::DynamicMatrix<R>& /*entityNeighborRet*/, Dune::DynamicMatrix<R>& /*neighborEntityRet*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  }

  template <class E, class D, int d, class R, class N, class IntersectionType>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, 1, 1>& localDiffusionFactorEntity,
                const Stuff::LocalfunctionInterface<E, D, d, R, d, d>& localDiffusionTensorEntity,
                const Stuff::LocalfunctionInterface<N, D, d, R, 1, 1>& localDiffusionFactorNeighbor,
                const Stuff::LocalfunctionInterface<N, D, d, R, d, d>& localDiffusionTensorNeighbor,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>& testBaseEntity,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface<N, D, d, R, 1, 1>& testBaseNeighbor,
                const Stuff::LocalfunctionSetInterface<N, D, d, R, 1, 1>& ansatzBaseNeighbor,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
                Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    // clear ret
    Stuff::Common::clear(entityEntityRet);
    Stuff::Common::clear(neighborNeighborRet);
    Stuff::Common::clear(entityNeighborRet);
    Stuff::Common::clear(neighborEntityRet);
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::DomainType DomainType;
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::RangeType RangeType;
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::JacobianRangeType JacobianRangeType;
    // convert local point (which is in intersection coordinates) to entity/neighbor coordinates
    const DomainType localPointEn    = intersection.geometryInInside().global(localPoint);
    const DomainType localPointNe    = intersection.geometryInOutside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    typedef Stuff::Common::FieldMatrix<R, d, d> TensorType;
    const auto local_diffusion_factor_en       = localDiffusionFactorEntity.evaluate(localPointEn);
    const TensorType local_diffusion_tensor_en = localDiffusionTensorEntity.evaluate(localPointEn);
    const auto local_diffusion_factor_ne       = localDiffusionFactorNeighbor.evaluate(localPointNe);
    const TensorType local_diffusion_tensor_ne = localDiffusionTensorNeighbor.evaluate(localPointNe);
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(testBaseEntity.order(),
                 std::max(ansatzBaseEntity.order(), std::max(testBaseNeighbor.order(), ansatzBaseNeighbor.order())));
    const R sigma = internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R delta_plus   = unitOuterNormal * (local_diffusion_tensor_ne * unitOuterNormal);
    const R delta_minus  = unitOuterNormal * (local_diffusion_tensor_en * unitOuterNormal);
    const R gamma        = (delta_plus * delta_minus) / (delta_plus + delta_minus);
    const R penalty      = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    const R weight_plus  = delta_minus / (delta_plus + delta_minus);
    const R weight_minus = delta_plus / (delta_plus + delta_minus);
    // compute diffusion value (should be factor * tensor, but this is the same)
    const TensorType diffusion_value_en = local_diffusion_tensor_en * local_diffusion_factor_en[0];
    const TensorType diffusion_value_ne = local_diffusion_tensor_ne * local_diffusion_factor_ne[0];
    // evaluate bases
    // * entity
    //   * test
    const size_t rowsEn = testBaseEntity.size();
    std::vector<RangeType> testValuesEn(rowsEn, RangeType(0));
    std::vector<JacobianRangeType> testGradientsEn(rowsEn, JacobianRangeType(0));
    testBaseEntity.evaluate(localPointEn, testValuesEn);
    testBaseEntity.jacobian(localPointEn, testGradientsEn);
    //   * ansatz
    const size_t colsEn = ansatzBaseEntity.size();
    std::vector<RangeType> ansatzValuesEn(colsEn, RangeType(0));
    std::vector<JacobianRangeType> ansatzGradientsEn(colsEn, JacobianRangeType(0));
    ansatzBaseEntity.evaluate(localPointEn, ansatzValuesEn);
    ansatzBaseEntity.jacobian(localPointEn, ansatzGradientsEn);
    // * neighbor
    //   * test
    const size_t rowsNe = testBaseNeighbor.size();
    std::vector<RangeType> testValuesNe(rowsNe, RangeType(0));
    std::vector<JacobianRangeType> testGradientsNe(rowsNe, JacobianRangeType(0));
    testBaseNeighbor.evaluate(localPointNe, testValuesNe);
    testBaseNeighbor.jacobian(localPointNe, testGradientsNe);
    //   * ansatz
    const size_t colsNe = ansatzBaseNeighbor.size();
    std::vector<RangeType> ansatzValuesNe(colsNe, RangeType(0));
    std::vector<JacobianRangeType> ansatzGradientsNe(colsNe, JacobianRangeType(0));
    ansatzBaseNeighbor.evaluate(localPointNe, ansatzValuesNe);
    ansatzBaseNeighbor.jacobian(localPointNe, ansatzGradientsNe);
    // compute the evaluations
    assert(entityEntityRet.rows() >= rowsEn);
    assert(entityEntityRet.cols() >= colsEn);
    assert(entityNeighborRet.rows() >= rowsEn);
    assert(entityNeighborRet.cols() >= colsNe);
    assert(neighborEntityRet.rows() >= rowsNe);
    assert(neighborEntityRet.cols() >= colsEn);
    assert(neighborNeighborRet.rows() >= rowsNe);
    assert(neighborNeighborRet.cols() >= colsNe);
    // loop over all entity test basis functions
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      auto& entityEntityRetRow   = entityEntityRet[ii];
      auto& entityNeighborRetRow = entityNeighborRet[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        // consistency term
        entityEntityRetRow[jj] +=
            -weight_minus * ((diffusion_value_en * ansatzGradientsEn[jj][0]) * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityEntityRetRow[jj] +=
            -weight_minus * ansatzValuesEn[jj] * ((diffusion_value_en * testGradientsEn[ii][0]) * unitOuterNormal);
        // penalty term
        entityEntityRetRow[jj] += penalty * ansatzValuesEn[jj] * testValuesEn[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        entityNeighborRetRow[jj] +=
            -weight_plus * ((diffusion_value_ne * ansatzGradientsNe[jj][0]) * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityNeighborRetRow[jj] +=
            weight_minus * ansatzValuesNe[jj] * ((diffusion_value_en * testGradientsEn[ii][0]) * unitOuterNormal);
        // penalty term
        entityNeighborRetRow[jj] += -1.0 * penalty * ansatzValuesNe[jj] * testValuesEn[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all entity test basis functions
    // loop over all neighbor test basis functions
    for (size_t ii = 0; ii < rowsNe; ++ii) {
      auto& neighborEntityRetRow   = neighborEntityRet[ii];
      auto& neighborNeighborRetRow = neighborNeighborRet[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        // consistency term
        neighborEntityRetRow[jj] +=
            weight_minus * ((diffusion_value_en * ansatzGradientsEn[jj][0]) * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborEntityRetRow[jj] +=
            -weight_plus * ansatzValuesEn[jj] * ((diffusion_value_ne * testGradientsNe[ii][0]) * unitOuterNormal);
        // penalty term
        neighborEntityRetRow[jj] += -1.0 * penalty * ansatzValuesEn[jj] * testValuesNe[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        neighborNeighborRetRow[jj] +=
            weight_plus * ((diffusion_value_ne * ansatzGradientsNe[jj][0]) * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborNeighborRetRow[jj] +=
            weight_plus * ansatzValuesNe[jj] * ((diffusion_value_ne * testGradientsNe[ii][0]) * unitOuterNormal);
        // penalty term
        neighborNeighborRetRow[jj] += penalty * ansatzValuesNe[jj] * testValuesNe[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // void evaluate< ..., 1, 1 >(...) const

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const double beta_;
}; // Inner


template <class DiffusionFactorType, class DiffusionTensorType>
class BoundaryLHS
    : public LocalEvaluation::Codim1Interface<internal::BoundaryLHSTraits<DiffusionFactorType, DiffusionTensorType>, 2>
{
public:
  typedef internal::BoundaryLHSTraits<DiffusionFactorType, DiffusionTensorType> Traits;

  BoundaryLHS(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
              const double beta = internal::default_beta(DiffusionFactorType::dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , beta_(beta)
  {
  }

  template <class EntityType>
  class LocalfunctionTuple
  {
    typedef typename DiffusionFactorType::LocalfunctionType LocalDiffusionFactorType;
    typedef typename DiffusionTensorType::LocalfunctionType LocalDiffusionTensorType;

  public:
    typedef std::tuple<std::shared_ptr<LocalDiffusionFactorType>, std::shared_ptr<LocalDiffusionTensorType>> Type;
  };

  template <class EntityType>
  typename LocalfunctionTuple<EntityType>::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity), diffusion_tensor_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class E, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  size_t order(const typename LocalfunctionTuple<E>::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBase) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    return order(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase);
  }

  template <class E, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rT, int rCT, int rA, int rCA>
  size_t order(const Stuff::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& localDiffusionFactor,
               const Stuff::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& localDiffusionTensor,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBase) const
  {
    return localDiffusionFactor.order() + localDiffusionTensor.order() + testBase.order() + ansatzBase.order();
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class E, class IntersectionType, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const typename LocalfunctionTuple<E>::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    evaluate(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase, intersection, localPoint, ret);
  }

  template <class E, class IntersectionType, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rT,
            int rCT, int rA, int rCA>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& /*localDiffusionFactor*/,
                const Stuff::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& /*localDiffusionTensor*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& /*testBase*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& /*ansatzBase*/,
                const IntersectionType& /*intersection*/, const Dune::FieldVector<D, d - 1>& /*localPoint*/,
                Dune::DynamicMatrix<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  }

  template <class E, class D, int d, class R, class IntersectionType>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, 1, 1>& localDiffusionFactor,
                const Stuff::LocalfunctionInterface<E, D, d, R, d, d>& localDiffusionTensor,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>& ansatzBase,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& ret) const
  {
    // clear ret
    Stuff::Common::clear(ret);
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::DomainType DomainType;
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::RangeType RangeType;
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::JacobianRangeType JacobianRangeType;
    // get local point (which is in intersection coordinates) in entity coordinates
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    typedef Stuff::Common::FieldMatrix<R, d, d> TensorType;
    const auto diffusion_factor_value       = localDiffusionFactor.evaluate(localPointEntity);
    const TensorType diffusion_tensor_value = localDiffusionTensor.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(testBase.order(), ansatzBase.order());
    const R sigma             = internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = unitOuterNormal * (diffusion_tensor_value * unitOuterNormal);
    const R penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // compute diffusion value (should be factor * tensor, but this is the same)
    const TensorType diffusion_value = diffusion_tensor_value * diffusion_factor_value;
    // evaluate bases
    // * test
    const size_t rows = testBase.size();
    std::vector<RangeType> testValues(rows, RangeType(0));
    std::vector<JacobianRangeType> testGradients(rows, JacobianRangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    testBase.jacobian(localPointEntity, testGradients);
    // * ansatz
    const size_t cols = ansatzBase.size();
    std::vector<RangeType> ansatzValues(cols, RangeType(0));
    std::vector<JacobianRangeType> ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.evaluate(localPointEntity, ansatzValues);
    ansatzBase.jacobian(localPointEntity, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    // loop over all test basis functions
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      // loop over all ansatz basis functions
      for (size_t jj = 0; jj < cols; ++jj) {
        // consistency term
        retRow[jj] += -1.0 * ((diffusion_value * ansatzGradients[jj][0]) * unitOuterNormal) * testValues[ii];
        // symmetry term
        retRow[jj] += -1.0 * ansatzValues[jj] * ((diffusion_value * testGradients[ii][0]) * unitOuterNormal);
        // penalty term
        retRow[jj] += penalty * ansatzValues[jj] * testValues[ii];
      } // loop over all ansatz basis functions
    } // loop over all test basis functions
  } // void evaluate(...) const

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const double beta_;
}; // class BoundaryLHS


template <class DiffusionFactorType, class DirichletType, class DiffusionTensorType>
class BoundaryRHS
    : public LocalEvaluation::
          Codim1Interface<internal::BoundaryRHSTraits<DiffusionFactorType, DirichletType, DiffusionTensorType>, 1>
{
public:
  typedef internal::BoundaryRHSTraits<DiffusionFactorType, DirichletType, DiffusionTensorType> Traits;

  BoundaryRHS(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
              const DirichletType& dirichlet,
              const double beta = internal::default_beta(DiffusionFactorType::dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , dirichlet_(dirichlet)
    , beta_(beta)
  {
  }

  template <class EntityType>
  class LocalfunctionTuple
  {
    typedef typename DiffusionFactorType::LocalfunctionType LocalDiffusionFactorType;
    typedef typename DiffusionTensorType::LocalfunctionType LocalDiffusionTensorType;
    typedef typename DirichletType::LocalfunctionType LocalDirichletType;

  public:
    typedef std::tuple<std::shared_ptr<LocalDiffusionFactorType>, std::shared_ptr<LocalDiffusionTensorType>,
                       std::shared_ptr<LocalDirichletType>> Type;
  };

  template <class EntityType>
  typename LocalfunctionTuple<EntityType>::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           dirichlet_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class E, class D, int d, class R, int r, int rC>
  size_t order(const typename LocalfunctionTuple<E>::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>& testBase) const
  {
    const auto localDiffusionFactor = std::get<0>(localFuncs);
    const auto localDiffusionTensor = std::get<1>(localFuncs);
    const auto localDirichlet = std::get<2>(localFuncs);
    return order(*localDiffusionFactor, *localDiffusionTensor, *localDirichlet, testBase);
  }

  template <class E, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rLR, int rCLR, int rT, int rCT>
  size_t order(const Stuff::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& localDiffusionFactor,
               const Stuff::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& localDiffusionTensor,
               const Stuff::LocalfunctionInterface<E, D, d, R, rLR, rCLR>& localDirichlet,
               const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase) const
  {
    const size_t testOrder         = testBase.order();
    const size_t testGradientOrder = std::max(ssize_t(testOrder) - 1, ssize_t(0));
    const size_t diffusionOrder    = localDiffusionFactor.order() + localDiffusionTensor.order();
    const size_t dirichletOrder = localDirichlet.order();
    return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
  } // static int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class E, class IntersectionType, class D, int d, class R, int r, int rC>
  void evaluate(const typename LocalfunctionTuple<E>::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, r, rC>& testBase,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicVector<R>& ret) const
  {
    const auto localDiffusionFactor = std::get<0>(localFuncs);
    const auto localDiffusionTensor = std::get<1>(localFuncs);
    const auto localDirichlet = std::get<2>(localFuncs);
    evaluate(*localDiffusionFactor, *localDiffusionTensor, *localDirichlet, testBase, intersection, localPoint, ret);
  }

  template <class E, class IntersectionType, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rLDR,
            int rCLDR, int rT, int rCT>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, rDF, rCDF>& /*localDiffusionFactor*/,
                const Stuff::LocalfunctionInterface<E, D, d, R, rDT, rCDT>& /*localDiffusionTensor*/,
                const Stuff::LocalfunctionInterface<E, D, d, R, rLDR, rCLDR>& /*localDirichlet*/,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& /*testBase*/,
                const IntersectionType& /*intersection*/, const Dune::FieldVector<D, d - 1>& /*localPoint*/,
                Dune::DynamicVector<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  }

  template <class E, class D, int d, class R, class IntersectionType>
  void evaluate(const Stuff::LocalfunctionInterface<E, D, d, R, 1, 1>& localDiffusionFactor,
                const Stuff::LocalfunctionInterface<E, D, d, R, d, d>& localDiffusionTensor,
                const Stuff::LocalfunctionInterface<E, D, d, R, 1, 1>& localDirichlet,
                const Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>& testBase,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicVector<R>& ret) const
  {
    // clear ret
    Stuff::Common::clear(ret);
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::DomainType DomainType;
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::RangeType RangeType;
    typedef typename Stuff::LocalfunctionSetInterface<E, D, d, R, 1, 1>::JacobianRangeType JacobianRangeType;
    // get local point (which is in intersection coordinates) in entity coordinates
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local functions
    typedef Stuff::Common::FieldMatrix<R, d, d> TensorType;
    const auto diffusionFactorValue       = localDiffusionFactor.evaluate(localPointEntity);
    const TensorType diffusionTensorValue = localDiffusionTensor.evaluate(localPointEntity);
    const RangeType dirichletValue        = localDirichlet.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t polorder = testBase.order();
    const R sigma         = internal::boundary_sigma(polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = unitOuterNormal * (diffusionTensorValue * unitOuterNormal);
    const R penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // compute diffusion value (should be factor * tensor, but this is the same)
    const TensorType diffusionValue = diffusionTensorValue * diffusionFactorValue;
    // evaluate basis
    const size_t size = testBase.size();
    std::vector<RangeType> testValues(size, RangeType(0));
    std::vector<JacobianRangeType> testGradients(size, JacobianRangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    testBase.jacobian(localPointEntity, testGradients);
    // compute
    assert(ret.size() >= size);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      // symmetry term
      ret[ii] += -1.0 * dirichletValue * ((diffusionValue * testGradients[ii][0]) * unitOuterNormal);
      // penalty term
      ret[ii] += penalty * dirichletValue * testValues[ii];
    } // loop over all test basis functions
  } // void evaluate(...) const

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const DirichletType& dirichlet_;
  const double beta_;
}; // class BoundaryRHS


} // namespace SWIPDG
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_LOCALEVALUATION_SWIPDG_HH
