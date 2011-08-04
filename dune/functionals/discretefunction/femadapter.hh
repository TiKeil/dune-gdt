#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTION_FEMADAPTER_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTION_FEMADAPTER_HH

// local includes
#include "local.hh"

namespace Dune {

namespace Functionals {

namespace DiscreteFunction {

template <class FucntionalsDiscreteFunctionImp>
class FemAdapter
{
public:
  typedef FucntionalsDiscreteFunctionImp HostDiscreteFunctionType;

  typedef typename HostDiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  enum
  {
    polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder
  };

  typedef FemAdapter<HostDiscreteFunctionType> ThisType;

  typedef Dune::Functionals::DiscreteFunction::LocalConst<ThisType> LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

  typedef typename HostDiscreteFunctionType::StorageType StorageType;

  FemAdapter(HostDiscreteFunctionType& hostDiscreteFunction)
    : hostDiscreteFunction_(hostDiscreteFunction)
  {
  }

  //! copy constructor
  FemAdapter(const ThisType& other)
    : hostDiscreteFunction_(other.hostDiscreteFunction())
  {
  }

  //! assignment operator
  ThisType& operator=(const ThisType& other)
  {
    if (this != other) {
      hostDiscreteFunction_ = other.hostDiscreteFunction();
    }
    return *this;
  }

  HostDiscreteFunctionType& hostDiscreteFunction()
  {
    return hostDiscreteFunction();
  }

  const HostDiscreteFunctionType& hostDiscreteFunction() const
  {
    return hostDiscreteFunction();
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return hostDiscreteFunction_.space();
  }

  const std::string name() const
  {
    return hostDiscreteFunction_.name();
  }

  void setName(const std::string& newName = "")
  {
    hostDiscreteFunction_.setName(newName);
  }

  void clear()
  {
    std::cout << "DiscreteFunction::FemAdapter::clear() does not do anything!" << std::endl;
  }

  const StorageType& storage() const
  {
    return hostDiscreteFunction_.storage();
  }

  StorageType& storage()
  {
    return hostDiscreteFunction_.storage();
  }

  //  RangeFieldType& operator[] ( const unsigned int globalDofNumber )
  //  {
  //    return hostDiscreteFunction_.operator[]( globalDofNumber );
  //  }

  const RangeFieldType& operator[](const unsigned int globalDofNumber) const
  {
    return hostDiscreteFunction_.operator[](globalDofNumber);
  }

  //  template< class EntityType >
  //  LocalFunctionType localFunction( const EntityType& entity )
  //  {
  //    return LocalFunctionType( (*this), entity );
  //  }

  template <class EntityType>
  const LocalFunctionType localFunction(const EntityType& entity) const
  {
    return LocalFunctionType((*this), entity);
  }

  bool continuous() const
  {
    return hostDiscreteFunction_.continuous();
  }

  int oder() const
  {
    return hostDiscreteFunction_.order();
  }

private:
  HostDiscreteFunctionType& hostDiscreteFunction_;

}; // end class FemAdapter

} // end namespace DiscreteFunction

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTION_FEMADAPTER_HH
