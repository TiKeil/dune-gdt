#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTEGRATOR_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTEGRATOR_HH

// dune fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>

// dune-fem-functionals includes
#include <dune/fem/localoperation/interface.hh>

namespace Dune {

namespace Functionals {

namespace LocalOperation {

template <class FunctionSpaceImp, class LocalOperationImp>
class VolumeIntegrator
    : public Dune::Functionals::LocalOperation::Interface<FunctionSpaceImp,
                                                          VolumeIntegrator<FunctionSpaceImp, LocalOperationImp>>
{
public:
  typedef typename Dune::Functionals::LocalOperation::
      TypeSelector<VolumeIntegrator<FunctionSpaceImp, LocalOperationImp>>::Select InterfaceType;

  typedef typename Dune::Functionals::LocalOperation::TypeSelector<LocalOperationImp>::Select LocalOperationType;

  typedef FunctionSpaceImp FunctionSpaceType;

  typedef typename InterfaceType::RangeFieldType RangeFieldType;

  typedef typename InterfaceType::RangeType RangeType;

  typedef typename InterfaceType::DomainType DomainType;

  VolumeIntegrator(const LocalOperationType& localOperation)
    : InterfaceType()
    , localOperation_(localOperation)
  {
    std::cout << "VolumeIntegrator::VolumeIntegrator()" << std::endl;
  }

  InterfaceType localOperation() const
  {
    return localOperation_;
  }

  /**
    * \tparam LocalFunctionType
    *         Should comply to the Dune::LocalFunction interface.
    **/
  template <class LocalFunctionType>
  RangeFieldType operateLocal(const LocalFunctionType& localFunction) const
  {
    // init return value
    RangeFieldType ret = 0.0;

    // some types we will need
    typedef typename LocalFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef typename LocalFunctionType::EntityType EntityType;

    typedef typename EntityType::Geometry EntityGeometryType;

    typedef CachingQuadrature<GridPartType, 0> VolumeQuadratureType;

    // entity and geometry
    const EntityType& entity                 = localFunction.entity();
    const EntityGeometryType& entityGeometry = entity.geometry();

    // quadrature
    const unsigned int quadratureOrder = localFunction.order();
    const VolumeQuadratureType volumeQuadrature(entity, quadratureOrder);
    const unsigned int numberOfQuadraturePoints = volumeQuadrature.nop();

    for (unsigned int q = 0; q < numberOfQuadraturePoints; ++q) {
      // local coordinate
      const DomainType x = volumeQuadrature.point(q);

      // integration factors
      const double integrationFactor = entityGeometry.integrationElement(x);
      const double quadratureWeight  = volumeQuadrature.weight(q);

      // evaluate the local operation
      const RangeFieldType localOperationEvalauted = localOperation_.evaluateLocal(localFunction, x);

      // compute integral
      ret += integrationFactor * quadratureWeight * localOperationEvalauted;
    }

    // return
    return ret;
  }

private:
  const LocalOperationType localOperation_;

}; // end class VolumeIntegrator

} // end namespace LocalOperation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTEGRATOR_HH
