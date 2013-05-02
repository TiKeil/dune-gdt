#ifndef DUNE_DETAILED_DISCRETIZATIONS_SPACE_LAGRANGE_CONTINUOUS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_SPACE_LAGRANGE_CONTINUOUS_HH

#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/fem_localfunctions/localfunctions/transformations.hh>
#include <dune/fem_localfunctions/basefunctions/genericbasefunctionsetstorage.hh>
#include <dune/fem_localfunctions/basefunctionsetmap/basefunctionsetmap.hh>
#include <dune/fem_localfunctions/space/genericdiscretefunctionspace.hh>

#include <dune/stuff/common/color.hh>

#include "../mapper/fem.hh"
#include "../basefunctionset/fem-localfunctions.hh"
#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDimRows = 1, int rangeDimCols = 1>
class ContinuousLagrangeSpace;


// forward, to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDimRows = 1, int rangeDimCols = 1>
class ContinuousLagrangeSpaceTraits;


/**
 *  \brief Traits class for ContinuousLagrangeSpace for dimRange 1x1.
 */
template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class ContinuousLagrangeSpaceTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
{
public:
  typedef GridPartImp GridPartType;
  static const int polOrder = polynomialOrder;
  dune_static_assert((polOrder >= 1), "ERROR: wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = 1;
  static const unsigned int dimRangeCols = 1;
  typedef ContinuousLagrangeSpace<GridPartType, polOrder, RangeFieldType, dimRangeRows, dimRangeCols> derived_type;

private:
  typedef ContinuousLagrangeSpaceTraits<GridPartType, polOrder, RangeFieldType, dimRangeRows, dimRangeCols> ThisType;
  typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, dimDomain, DomainFieldType, RangeFieldType>
      FiniteElementType;
  typedef Dune::FemLocalFunctions::BaseFunctionSetMap<GridPartType, FiniteElementType,
                                                      Dune::FemLocalFunctions::NoTransformation,
                                                      Dune::FemLocalFunctions::SimpleStorage, polOrder,
                                                      polOrder> BaseFunctionSetMapType;

public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace<BaseFunctionSetMapType> BackendType;
  typedef MapperWrappedFemDofMapper<typename BackendType::MapperType> MapperType;
  typedef BaseFunctionSetFemLocalfunctionsWrapper<BaseFunctionSetMapType> BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType EntityType;

private:
  friend class ContinuousLagrangeSpace<GridPartType, polOrder, RangeFieldType, dimRangeRows, dimRangeCols>;
}; // class ContinuousLagrangeSpaceTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >


template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class ContinuousLagrangeSpace<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public SpaceInterface<ContinuousLagrangeSpaceTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>>
{
public:
  typedef ContinuousLagrangeSpaceTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridPartType GridPartType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const int polOrder           = Traits::polOrder;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = Traits::dimRangeRows;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  ContinuousLagrangeSpace(const GridPartType& gridP)
    : gridPart_(checkGridPart(gridP))
    , baseFunctionSetMap_(gridPart_)
    , backend_(const_cast<GridPartType&>(gridPart_), baseFunctionSetMap_)
    , mapper_(backend_.mapper())
  {
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(baseFunctionSetMap_, entity);
  }

private:
  static const GridPartType& checkGridPart(const GridPartType& gP)
  {
    // check
    typedef typename Dune::Fem::AllGeomTypes<typename GridPartType::IndexSetType, typename GridPartType::GridType>
        AllGeometryTypes;
    const AllGeometryTypes allGeometryTypes(gP.indexSet());
    const std::vector<Dune::GeometryType>& geometryTypes = allGeometryTypes.geomTypes(0);
    if (!(geometryTypes.size() == 1 && geometryTypes[0].isSimplex()))
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                      << " this space is only implemented for simplicial grids!");
    return gP;
  }

  const GridPartType& gridPart_;
  BaseFunctionSetMapType baseFunctionSetMap_;
  const BackendType backend_;
  const MapperType mapper_;
}; // class ContinuousLagrangeSpace< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_SPACE_LAGRANGE_CONTINUOUS_HH
