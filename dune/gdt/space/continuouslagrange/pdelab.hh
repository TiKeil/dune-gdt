// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_CONTINUOUSLAGRANGE_PDELAB_HH
#define DUNE_GDT_SPACE_CONTINUOUSLAGRANGE_PDELAB_HH

#include <memory>

#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

#include "../../mapper/pdelab.hh"
#include "../../basefunctionset/pdelab.hh"

#include "../interface.hh"
#include "../constraints.hh"

namespace Dune {
namespace GDT {
namespace ContinuousLagrangeSpace {


// forward, to be used in the traits and to allow for specialization
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class PdelabWrapper;


/**
 *  \brief Traits class for ContinuousLagrangeSpace::PdelabWrapper.
 */
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class PdelabWrapperTraits
{
public:
  typedef PdelabWrapper<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridPartImp GridPartType;
  static const int polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;

private:
  template <class G, bool single_geom, bool is_simplex, bool is_cube>
  struct Choose
  {
    static_assert(Dune::AlwaysFalse<G>::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template <class G>
  struct Choose<G, true, true, false>
  {
    typedef PDELab::PkLocalFiniteElementMap<typename GridPartType::GridViewType, DomainFieldType, RangeFieldType,
                                            polOrder> FEMapType;
  };
  template <class G>
  struct Choose<G, true, false, true>
  {
    typedef PDELab::QkLocalFiniteElementMap<typename GridPartType::GridViewType, DomainFieldType, RangeFieldType,
                                            polOrder> FEMapType;
  };
  typedef typename GridPartType::GridType GridType;
  typedef typename Choose<GridType, Dune::Capabilities::hasSingleGeometryType<GridType>::v,
                          Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                              == GenericGeometry::SimplexTopology<dimDomain>::type::id,
                          Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                              == GenericGeometry::CubeTopology<dimDomain>::type::id>::FEMapType FEMapType;

public:
  typedef PDELab::GridFunctionSpace<typename GridPartType::GridViewType, FEMapType> BackendType;
  typedef Mapper::PdelabWrapper<BackendType> MapperType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef BaseFunctionSet::PdelabWrapper<BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange,
                                         dimRangeCols> BaseFunctionSetType;

private:
  friend class PdelabWrapper<GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>;
}; // class SpaceWrappedFemContinuousLagrangeTraits


// unspecialized version, to give better compile errors
template <class GP, int p, class R, int r, int rC>
class PdelabWrapper : public SpaceInterface<PdelabWrapperTraits<GP, p, R, r, rC>>
{
public:
  typedef PdelabWrapperTraits<GP, p, R, r, rC> Traits;

  typedef typename Traits::GridPartType GridPartType;
  static const int polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = Traits::dimRange;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  PdelabWrapper(const std::shared_ptr<const GridPartType>& /*grid_prt*/)
  {
    static_assert((Dune::AlwaysFalse<GP>::value), "Not yet implemented for this combination of dimensions!");
  }
}; // PdelabWrapper


template <class GridPartImp, int polynomialOrder, class RangeFieldImp>
class PdelabWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public SpaceInterface<PdelabWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>>
{
  typedef SpaceInterface<PdelabWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1>> BaseType;
  typedef PdelabWrapper<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;

public:
  typedef PdelabWrapperTraits<GridPartImp, polynomialOrder, RangeFieldImp, 1, 1> Traits;

  typedef typename Traits::GridPartType GridPartType;
  static const int polOrder = Traits::polOrder;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = Traits::dimRange;
  static const unsigned int dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::MapperType MapperType;
  typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

  PdelabWrapper(const std::shared_ptr<const GridPartType>& gridP)
    : gridPart_(gridP)
    , fe_map_(std::make_shared<FEMapType>(gridPart_->gridView()))
    , backend_(std::make_shared<BackendType>(const_cast<GridPartType&>(*gridPart_).gridView(), *fe_map_))
    , mapper_(std::make_shared<MapperType>(*backend_))
    , tmpMappedRows_(mapper_->maxNumDofs())
    , tmpMappedCols_(mapper_->maxNumDofs())
  {
  }

  PdelabWrapper(const ThisType& other)
    : gridPart_(other.gridPart_)
    , fe_map_(other.fe_map_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
    , tmpMappedRows_(mapper_->maxNumDofs())
    , tmpMappedCols_(mapper_->maxNumDofs())
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridPart_ = other.gridPart_;
      fe_map_   = other.fe_map_;
      backend_  = other.backend_;
      mapper_ = other.mapper_;
      tmpMappedRows_.resize(mapper_->maxNumDofs());
      tmpMappedCols_.resize(mapper_->maxNumDofs());
    }
    return *this;
  }

  ~PdelabWrapper()
  {
  }

  const std::shared_ptr<const GridPartType>& gridPart() const
  {
    return gridPart_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(*backend_, entity);
  }

  template <class R>
  void localConstraints(const EntityType& /*entity*/, Constraints::LocalDefault<R>& /*ret*/) const
  {
    static_assert((Dune::AlwaysFalse<R>::value), "Not implemented for arbitrary constraints!");
  }

  void
  localConstraints(const EntityType& /*entity*/,
                   Constraints::Dirichlet<typename GridPartType::IntersectionType, RangeFieldType, true>& /*ret*/) const
  {
    //    std::set< size_t > localDirichletDofs;
    //    const auto& gridBoundary = ret.gridBoundary();
    //    if (polOrder == 1) {
    //      localDirichletDofs = this->findLocalDirichletDoFs(entity, gridBoundary);
    //    } else {
    //      const auto& lagrangePointSet = backend_->lagrangePointSet(entity);
    //      // loop over all intersections
    //      const auto intersectionEndIt = gridPart_->iend(entity);
    //      for (auto intersectionIt = gridPart_->ibegin(entity); intersectionIt !=intersectionEndIt; ++intersectionIt)
    //      {
    //        const auto& intersection = *intersectionIt;
    //        // only work on dirichlet intersections
    //        if (gridBoundary.dirichlet(intersection)) {
    //          // get local face number of boundary intersection
    //          const int intersectionIndex = intersection.indexInInside();
    //          // iterate over face dofs and set unit row
    //          const auto faceDofEndIt = lagrangePointSet.template endSubEntity< 1 >(intersectionIndex);
    //          for (auto faceDofIt = lagrangePointSet.template beginSubEntity< 1 >(intersectionIndex);
    //               faceDofIt != faceDofEndIt;
    //               ++faceDofIt) {
    //            const size_t localDofIndex = *faceDofIt;
    //            localDirichletDofs.insert(localDofIndex);
    //          } // iterate over face dofs and set unit row
    //        } // only work on dirichlet intersections
    //      } // loop over all intersections
    //    }
    //    const size_t numRows = localDirichletDofs.size();
    //    if (numRows > 0) {
    //      const size_t numCols = mapper_->numDofs(entity);
    //      ret.setSize(numRows, numCols);
    //      mapper_->globalIndices(entity, tmpMappedRows_);
    //      mapper_->globalIndices(entity, tmpMappedCols_);
    //      size_t localRow = 0;
    //      const RangeFieldType zero(0);
    //      const RangeFieldType one(1);
    //      for (auto localDirichletDofIt = localDirichletDofs.begin();
    //           localDirichletDofIt != localDirichletDofs.end();
    //           ++localDirichletDofIt) {
    //        const size_t& localDirichletDofIndex = * localDirichletDofIt;
    //        ret.globalRow(localRow) = tmpMappedRows_[localDirichletDofIndex];
    //        for (size_t jj = 0; jj < ret.cols(); ++jj) {
    //          ret.globalCol(jj) = tmpMappedCols_[jj];
    //          if (tmpMappedCols_[jj] == tmpMappedRows_[localDirichletDofIndex])
    //            ret.value(localRow, jj) = one;
    //          else
    //            ret.value(localRow, jj) = zero;
    //        }
    //        ++localRow;
    //      }
    //    } else {
    //      ret.setSize(0, 0);
    //    }
  } // ... localConstraints(..., Dirichlet< ..., true >)

  void localConstraints(
      const EntityType& /*entity*/,
      Constraints::Dirichlet<typename GridPartType::IntersectionType, RangeFieldType, false>& /*ret*/) const
  {
    //    std::set< size_t > localDirichletDofs;
    //    const auto& gridBoundary = ret.gridBoundary();
    //    if (polOrder == 1) {
    //      localDirichletDofs = this->findLocalDirichletDoFs(entity, gridBoundary);
    //    } else {
    //      const auto& lagrangePointSet = backend_->lagrangePointSet(entity);
    //      // loop over all intersections
    //      const auto intersectionEndIt = gridPart_->iend(entity);
    //      for (auto intersectionIt = gridPart_->ibegin(entity); intersectionIt !=intersectionEndIt; ++intersectionIt)
    //      {
    //        const auto& intersection = *intersectionIt;
    //        // only work on dirichlet intersections
    //        if (gridBoundary.dirichlet(intersection)) {
    //          // get local face number of boundary intersection
    //          const int intersectionIndex = intersection.indexInInside();
    //          // iterate over face dofs and set unit row
    //          const auto faceDofEndIt = lagrangePointSet.template endSubEntity< 1 >(intersectionIndex);
    //          for (auto faceDofIt = lagrangePointSet.template beginSubEntity< 1 >(intersectionIndex);
    //               faceDofIt != faceDofEndIt;
    //               ++faceDofIt) {
    //            const size_t localDofIndex = *faceDofIt;
    //            localDirichletDofs.insert(localDofIndex);
    //          } // iterate over face dofs and set unit row
    //        } // only work on dirichlet intersections
    //      } // loop over all intersections
    //    }
    //    const size_t numRows = localDirichletDofs.size();
    //    if (numRows > 0) {
    //      const size_t numCols = mapper_->numDofs(entity);
    //      ret.setSize(numRows, numCols);
    //      mapper_->globalIndices(entity, tmpMappedRows_);
    //      mapper_->globalIndices(entity, tmpMappedCols_);
    //      size_t localRow = 0;
    //      const RangeFieldType zero(0);
    //      for (auto localDirichletDofIt = localDirichletDofs.begin();
    //           localDirichletDofIt != localDirichletDofs.end();
    //           ++localDirichletDofIt) {
    //        const size_t& localDirichletDofIndex = * localDirichletDofIt;
    //        ret.globalRow(localRow) = tmpMappedRows_[localDirichletDofIndex];
    //        for (size_t jj = 0; jj < ret.cols(); ++jj) {
    //          ret.globalCol(jj) = tmpMappedCols_[jj];
    //          ret.value(localRow, jj) = zero;
    //        }
    //        ++localRow;
    //      }
    //    } else {
    //      ret.setSize(0, 0);
    //    }
  } // ... localConstraints(..., Dirichlet< ..., false >)

  using BaseType::computePattern;

  template <class LocalGridPartType, class OtherSpaceType>
  PatternType* computePattern(const LocalGridPartType& localGridPart, const OtherSpaceType& otherSpace) const
  {
    return BaseType::computeCodim0Pattern(localGridPart, otherSpace);
  }

private:
  std::shared_ptr<const GridPartType> gridPart_;
  std::shared_ptr<const FEMapType> fe_map_;
  std::shared_ptr<const BackendType> backend_;
  std::shared_ptr<const MapperType> mapper_;
  mutable Dune::DynamicVector<size_t> tmpMappedRows_;
  mutable Dune::DynamicVector<size_t> tmpMappedCols_;
}; // class PdelabWrapper< ..., 1 >


} // namespace ContinuousLagrangeSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_CONTINUOUSLAGRANGE_PDELAB_HH
