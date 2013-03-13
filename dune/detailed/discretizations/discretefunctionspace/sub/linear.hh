#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

#include <memory>

#include <dune/detailed/discretizations/constraints/dirichlet.hh>

#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace DiscreteFunctionSpace {

namespace Sub {

namespace Linear {

template <class SuperSpaceImp>
class Dirichlet
{
public:
  typedef SuperSpaceImp SuperSpaceType;

  typedef Dirichlet<SuperSpaceType> ThisType;

  typedef typename SuperSpaceType::GridPartType GridPartType;

  typedef typename SuperSpaceType::GridViewType GridViewType;

  typedef typename Dune::Stuff::Grid::BoundaryInfo::Interface<GridViewType> BoundaryInfoType;

  typedef Dune::Detailed::Discretizations::Constraints::DirichletZero<SuperSpaceType> ConstraintsType;

  typedef typename SuperSpaceType::FunctionSpaceType FunctionSpaceType;

  static const int polynomialOrder = SuperSpaceType::polynomialOrder;

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

  typedef typename SuperSpaceType::PatternType PatternType;

  typedef typename Dune::Stuff::Grid::BoundaryInfo::AllDirichlet<GridViewType> DefaultBoundaryInfoType;

  Dirichlet(const SuperSpaceType& superSpace,
            const std::shared_ptr<const BoundaryInfoType> boundaryInfo =
                std::shared_ptr<const DefaultBoundaryInfoType>(new DefaultBoundaryInfoType(DefaultBoundaryInfoType())))
    : superSpace_(superSpace)
    , boundaryInfo_(boundaryInfo)
    , constraints_(superSpace_, *boundaryInfo_)
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

  const BoundaryInfoType& boundaryInfo() const
  {
    return *boundaryInfo_;
  }

  const GridPartType& gridPart() const
  {
    return superSpace_.gridPart();
  }

  const GridViewType& gridView() const
  {
    return superSpace_.gridView();
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

  template <class LocalGridPartType, class OtherDiscreteFunctionSpaceType>
  std::shared_ptr<const PatternType> computeLocalPattern(const LocalGridPartType& localGridPart,
                                                         const OtherDiscreteFunctionSpaceType& other) const
  {
    return superSpace_.computeLocalPattern(localGridPart, other);
  }

  template <class LocalGridPartType>
  std::shared_ptr<const PatternType> computeLocalPattern(const LocalGridPartType& localGridPart) const
  {
    return superSpace_.computeLocalPattern(localGridPart);
  }

  template <class CouplingGridPartType, class OutsideDiscreteFunctionSpaceType>
  std::shared_ptr<const PatternType> computeCouplingPattern(const CouplingGridPartType& couplingGridPart,
                                                            const OutsideDiscreteFunctionSpaceType& outerSpace) const
  {
    return superSpace_.computeCouplingPattern(couplingGridPart, outerSpace);
  }

  template <class OtherDiscreteFunctionSpaceType>
  std::shared_ptr<const PatternType> computePattern(const OtherDiscreteFunctionSpaceType& other) const
  {
    return superSpace_.computePattern(other);
  }

  std::shared_ptr<const PatternType> computePattern() const
  {
    return superSpace_.computePattern();
  }

private:
  //! copy constructor
  Dirichlet(const ThisType&);

  //! assignment operator
  ThisType& operator=(const ThisType&);

  const SuperSpaceType& superSpace_;
  const std::shared_ptr<const BoundaryInfoType> boundaryInfo_;
  const ConstraintsType constraints_;

}; // end class Dirichlet

} // end namespace Linear

} // end namespace Sub

} // end namespace DiscreteFunctionSpace

} // namespace Discretizations

} // end namespace Discretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
