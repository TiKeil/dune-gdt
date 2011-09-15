#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/basefunctionset/local/lagrange.hh>
#include <dune/detailed-discretizations/constraints/dirichlet.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace DiscreteFunctionSpace {

namespace Subspace {

namespace Linear {

template <class SuperSpaceImp>
class Dirichlet
{
public:
  typedef SuperSpaceImp SuperSpaceType;

  typedef Dirichlet<SuperSpaceType> ThisType;

  typedef Dune::DetailedDiscretizations::Constraints::DirichletZero<SuperSpaceType> ConstraintsType;

  typedef typename SuperSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename SuperSpaceType::GridPartType GridPartType;

  enum
  {
    polynomialOrder = SuperSpaceType::polynomialOrder
  };

  typedef typename SuperSpaceType::MapperType MapperType;

  typedef typename SuperSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  static const unsigned int dimDomain = SuperSpaceType::dimDomain;

  static const unsigned int dimRange = SuperSpaceType::dimRange;

  /**
      @name Convenience
      @{
   **/
  typedef typename SuperSpaceType::IteratorType IteratorType;

  typedef typename SuperSpaceType::EntityType EntityType;
  /**
      @}
   **/

  Dirichlet(const SuperSpaceType& superSpace)
    : superSpace_(superSpace)
    , constraints_(superSpace_)
  {
  }

  const SuperSpaceType& superSpace() const
  {
    return superSpace_;
  }

  const ConstraintsType& constraints() const
  {
    return constraints_;
  }

  const GridPartType& gridPart() const
  {
    return superSpace_.gridPart();
  }

  const MapperType& map() const
  {
    return superSpace_.map();
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return superSpace_.baseFunctionSet();
  }

  int order() const
  {
    return superSpace_.order();
  }

  bool continuous() const
  {
    return superSpace_.continuous();
  }

  /**
      @name Convenience methods
      @{
   **/
  IteratorType begin() const
  {
    return superSpace_.gridPart().template begin<0>();
  }

  IteratorType end() const
  {
    return superSpace_.gridPart().template end<0>();
  }
  /**
      @}
   **/

private:
  //! copy constructor
  Dirichlet(const ThisType&);

  //! assignment operator
  ThisType& operator=(const ThisType&);

  const SuperSpaceType& superSpace_;
  const ConstraintsType constraints_;

}; // end class Dirichlet

} // end namespace Linear

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
