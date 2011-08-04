#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTION_CONTINUOUS_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTION_CONTINUOUS_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-istl includes
#include <dune/istl/bvector.hh>

// dune-fem-tools includes
#include <dune/fem-tools/space/projection.hh>

// local includes
#include "local.hh"

namespace Dune {

namespace Functionals {

namespace DiscreteFunction {

namespace Continuous {

template <class ContinuousDiscreteFunctionSpaceImp>
class BlockVector
{
public:
  typedef ContinuousDiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  enum
  {
    polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder
  };

  typedef BlockVector<DiscreteFunctionSpaceType> ThisType;

  typedef Dune::Functionals::DiscreteFunction::Local<ThisType> LocalFunctionType;

  typedef Dune::Functionals::DiscreteFunction::LocalConst<ThisType> ConstLocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

  typedef Dune::BlockVector<Dune::FieldVector<RangeFieldType, dimRange>> StorageType;

  BlockVector(const DiscreteFunctionSpaceType& discreteFunctionSpace,
              const std::string name = "continuousBlockVectorFunction")
    : space_(discreteFunctionSpace)
    , storage_(discreteFunctionSpace.map().size())
    , name_(name)
  {
  }

  BlockVector(const DiscreteFunctionSpaceType& discreteFunctionSpace, const StorageType& storage,
              const std::string name = "continuousBlockVectorFunction")
    : space_(discreteFunctionSpace)
    , storage_(space_.map().size())
    , name_(name)
  {
    assert(storage.size() == storage_.size());
    for (unsigned int i = 0; i < storage_.size(); ++i) {
      storage_[i] = storage[i];
    }
  }

  template <class FunctionType>
  BlockVector(const DiscreteFunctionSpaceType& discreteFunctionSpace, const std::string name,
              const FunctionType& function, const std::string projectionType)
    : space_(discreteFunctionSpace)
    , storage_(space_.map().size())
    , name_(name)
  {
    if (projectionType.compare("dirichlet") == 0) {
      Dune::FemTools::Projection::Dirichlet::project(function, *this);
    } else {
      throw Dune::NotImplemented();
    }
  }

  //! copy constructor
  BlockVector(const ThisType& other)
    : space_(other.space())
    , storage_(space_.map().size())
    , name_("copyOF" + other.name())
  {
    for (unsigned int i = 0; i < storage_.size(); ++i) {
      operator[](i) = other[i];
    }
  }

  //! assignment operator
  ThisType& operator=(const ThisType& other)
  {
    if (this != other) {
      // we should do something like
      //      assert( this->space() == other.space() );
      assert(other.space().map().size() == this->space().map().size());
      for (unsigned int i = 0; i < storage_.size(); ++i) {
        operator[](i) = other[i];
      }
    }
    return *this;
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  const std::string name() const
  {
    return name_;
  }

  void setName(const std::string& newName = "")
  {
    name_ = newName;
  }

  void clear()
  {
    std::cout << "DiscreteFunction::Continuous::BlockVector::clear() does not do anything!" << std::endl;
  }

  const StorageType& storage() const
  {
    return storage_;
  }

  StorageType& storage()
  {
    return storage_;
  }

  RangeFieldType& operator[](const unsigned int globalDofNumber)
  {
    return storage_[globalDofNumber][0];
  }

  const RangeFieldType& operator[](const unsigned int globalDofNumber) const
  {
    return storage_[globalDofNumber][0];
  }

  template <class EntityType>
  LocalFunctionType localFunction(const EntityType& entity)
  {
    return LocalFunctionType((*this), entity);
  }

  template <class EntityType>
  const ConstLocalFunctionType localFunction(const EntityType& entity) const
  {
    return ConstLocalFunctionType((*this), entity);
  }

  /**
    \attention  This is not correct for order 0
    \todo       fix me
    **/
  bool continuous() const
  {
    return true;
  }

  int oder() const
  {
    return space_.order();
  }

private:
  const DiscreteFunctionSpaceType& space_;
  StorageType storage_;
  std::string name_;

}; // end class BlockVector

} // end namespace Continuous

} // end namespace DiscreteFunction

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTION_CONTINUOUS_HH
