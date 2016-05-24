// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_DUNE_PDELAB_WRAPPER_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_DUNE_PDELAB_WRAPPER_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#endif

#include <dune/stuff/common/type_utils.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {

#if HAVE_DUNE_PDELAB


// forwards, to be used in the traits and to allow for specialization
template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "Untested for arbitrary dimension!");
};


template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class PiolaTransformedDunePdelabWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "Untested for these dimensions!");
};

//new for Nedelec
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class Edges05DunePdelabWrapper
{
  static_assert(Dune::AlwaysFalse< PdelabSpaceImp >::value, "Untested for these dimensions!");
};

namespace internal {


// forward, to allow for specialization
template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols>
class DunePdelabWrapperTraits;


//! Specialization for rangeDimCols = 1
template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim>
class DunePdelabWrapperTraits<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef DunePdelabWrapper<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      derived_type;

private:
  typedef PDELab::LocalFunctionSpace<PdelabSpaceImp, PDELab::TrialSpaceTag> PdelabLFSType;
  typedef FiniteElementInterfaceSwitch<typename PdelabSpaceImp::Traits::FiniteElementType> FESwitchType;

public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;

private:
  friend class DunePdelabWrapper<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>;
};


template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim>
class PiolaTransformedDunePdelabWrapperTraits
{
  static_assert(domainDim == rangeDim, "Untested!");

public:
  typedef PiolaTransformedDunePdelabWrapper<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                            rangeDim> derived_type;

private:
  typedef PDELab::LocalFunctionSpace<PdelabSpaceImp, PDELab::TrialSpaceTag> PdelabLFSType;
  typedef FiniteElementInterfaceSwitch<typename PdelabSpaceImp::Traits::FiniteElementType> FESwitchType;

public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;

private:
  friend class PiolaTransformedDunePdelabWrapper<PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                                 rangeDim, 1>;
};

//new for nedelec
/**  Traits for the class Edges05DunePdelabWrapper
 *
 * Edges05DunePdelabWrapper
*/
template< class PdelabSpaceImp, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim >
class Edges05DunePdelabWrapperTraits
{
  static_assert(domainDim == rangeDim, "Untested!");
  static_assert(domainDim == 3, "Untested!");
public:
  typedef Edges05DunePdelabWrapper< PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim >
      derived_type;
private:
  typedef PDELab::LocalFunctionSpace< PdelabSpaceImp, PDELab::TrialSpaceTag > PdelabLFSType;
  typedef FiniteElementInterfaceSwitch< typename PdelabSpaceImp::Traits::FiniteElementType > FESwitchType;
public:
  typedef typename FESwitchType::Basis BackendType;
  typedef EntityImp EntityType;
private:
  friend class Edges05DunePdelabWrapper< PdelabSpaceImp, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >;
};

} // namespace internal


//! Specialization for dimRange = 1, dimRangeRows = 1
template <class PdelabSpaceType, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp>
class DunePdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public BaseFunctionSetInterface<internal::DunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp,
                                                                        domainDim, RangeFieldImp, 1, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
{
  typedef DunePdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::DunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp,
                                                                     domainDim, RangeFieldImp, 1, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, 1, 1> BaseType;

public:
  typedef internal::DunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
      Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType FESwitchType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  DunePdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , tmp_domain_(0)
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_     = std::unique_ptr<PdelabLFSType>(lfs_ptr);
    backend_ = std::unique_ptr<BackendType>(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
  } // DunePdelabWrapper(...)

  DunePdelabWrapper(ThisType&& source) = default;
  DunePdelabWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override final
  {
    return backend_->size();
  }

  virtual size_t order() const override final
  {
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const override final
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, ret);
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const override final
  {
    assert(ret.size() >= backend_->size());
    backend_->evaluateJacobian(xx, ret);
    const auto jacobian_inverse_transposed = this->entity().geometry().jacobianInverseTransposed(xx);
    for (size_t ii = 0; ii < ret.size(); ++ii) {
      jacobian_inverse_transposed.mv(ret[ii][0], tmp_domain_);
      ret[ii][0] = tmp_domain_;
    }
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  mutable DomainType tmp_domain_;
  std::unique_ptr<const PdelabLFSType> lfs_;
  std::unique_ptr<const BackendType> backend_;
}; // class DunePdelabWrapper


template <class PdelabSpaceType, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim>
class DunePdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public BaseFunctionSetInterface<internal::DunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp,
                                                                        domainDim, RangeFieldImp, rangeDim, 1>,
                                      DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef DunePdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::DunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp,
                                                                     domainDim, RangeFieldImp, rangeDim, 1>,
                                   DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

public:
  typedef internal::DunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                            rangeDim, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType FESwitchType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange  = rangeDim;

  DunePdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , tmp_domain_(0)
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_     = std::unique_ptr<PdelabLFSType>(lfs_ptr);
    backend_ = std::unique_ptr<BackendType>(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
  } // DunePdelabWrapper(...)

  DunePdelabWrapper(ThisType&& source) = default;
  DunePdelabWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override final
  {
    return backend_->size() * dimRange;
  }

  virtual size_t order() const override final
  {
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const override final
  {
    assert(ret.size() >= size());
    const size_t size_of_factor_basefunctionset = backend_->size();
    std::vector<Dune::FieldVector<RangeFieldImp, 1>> factor_ret(size_of_factor_basefunctionset);
    backend_->evaluateFunction(xx, factor_ret);
    // if factor_ret is [1 2] and we have two factors, we want to return [[1 0] [2 0] [0 1] [0 2]]
    for (size_t jj = 0; jj < size_of_factor_basefunctionset; ++jj) {
      for (size_t ii = 0; ii < dimRange; ++ii) {
        ret[ii * size_of_factor_basefunctionset + jj] *= 0;
        ret[ii * size_of_factor_basefunctionset + jj][ii] = factor_ret[jj][0];
      }
    }
  }

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const override final
  {
    assert(ret.size() >= backend_->size());
    const size_t size_of_factor_basefunctionset = backend_->size();
    std::vector<FieldMatrix<RangeFieldImp, 1, dimDomain>> factor_ret(size_of_factor_basefunctionset);
    backend_->evaluateJacobian(xx, factor_ret);
    const auto jacobian_inverse_transposed = this->entity().geometry().jacobianInverseTransposed(xx);
    for (size_t jj = 0; jj < factor_ret.size(); ++jj) {
      jacobian_inverse_transposed.mv(factor_ret[jj][0], tmp_domain_);
      factor_ret[jj][0] = tmp_domain_;
    }
    for (size_t ii = 0; ii < dimRange; ++ii) {
      for (size_t jj = 0; jj < size_of_factor_basefunctionset; ++jj) {
        ret[ii * size_of_factor_basefunctionset + jj] *= 0;
        ret[ii * size_of_factor_basefunctionset + jj][ii] = factor_ret[jj][0];
      }
    }
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  mutable DomainType tmp_domain_;
  std::unique_ptr<const PdelabLFSType> lfs_;
  std::unique_ptr<const BackendType> backend_;
}; // class DunePdelabWrapper


template <class PdelabSpaceType, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim>
class PiolaTransformedDunePdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                        1>
    : public BaseFunctionSetInterface<internal::PiolaTransformedDunePdelabWrapperTraits<PdelabSpaceType, EntityImp,
                                                                                        DomainFieldImp, domainDim,
                                                                                        RangeFieldImp, rangeDim>,
                                      DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef PiolaTransformedDunePdelabWrapper<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                            rangeDim, 1> ThisType;
  typedef BaseFunctionSetInterface<internal::PiolaTransformedDunePdelabWrapperTraits<PdelabSpaceType, EntityImp,
                                                                                     DomainFieldImp, domainDim,
                                                                                     RangeFieldImp, rangeDim>,
                                   DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

public:
  typedef internal::PiolaTransformedDunePdelabWrapperTraits<PdelabSpaceType, EntityImp, DomainFieldImp, domainDim,
                                                            RangeFieldImp, rangeDim> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType FESwitchType;

public:
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  PiolaTransformedDunePdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , tmp_domain_(DomainFieldType(0))
    , tmp_jacobian_transposed_(DomainFieldType(0))
    , tmp_jacobian_inverse_transposed_(DomainFieldType(0))
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_                 = std::unique_ptr<PdelabLFSType>(lfs_ptr);
    backend_             = std::unique_ptr<BackendType>(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
    tmp_ranges_          = std::vector<RangeType>(backend_->size(), RangeType(0));
    tmp_jacobian_ranges_ = std::vector<JacobianRangeType>(backend_->size(), JacobianRangeType(0));
  } // DunePdelabWrapper(...)

  PiolaTransformedDunePdelabWrapper(ThisType&& source) = default;

  PiolaTransformedDunePdelabWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override final
  {
    return backend_->size();
  }

  virtual size_t order() const override final
  {
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const override final
  {
    assert(lfs_);
    assert(backend_);
    assert(tmp_ranges_.size() >= backend_->size());
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, tmp_ranges_);
    const auto geometry                       = this->entity().geometry();
    tmp_jacobian_transposed_                  = geometry.jacobianTransposed(xx);
    const DomainFieldType integration_element = geometry.integrationElement(xx);
    for (size_t ii = 0; ii < backend_->size(); ++ii) {
      tmp_jacobian_transposed_.mtv(tmp_ranges_[ii], ret[ii]);
      ret[ii] /= integration_element;
    }
  } // ... evaluate(...)

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const override final
  {
    assert(lfs_);
    assert(backend_);
    assert(ret.size() >= backend_->size());
    backend_->evaluateJacobian(xx, tmp_jacobian_ranges_);
    const auto geometry                       = this->entity().geometry();
    tmp_jacobian_transposed_                  = geometry.jacobianTransposed(xx);
    tmp_jacobian_inverse_transposed_          = geometry.jacobianInverseTransposed(xx);
    const DomainFieldType integration_element = geometry.integrationElement(xx);
    for (size_t ii = 0; ii < backend_->size(); ++ii) {
      for (size_t jj = 0; jj < dimDomain; ++jj) {
        tmp_jacobian_inverse_transposed_.mv(tmp_jacobian_ranges_[ii][jj], ret[ii][jj]);
        tmp_jacobian_transposed_.mv(ret[ii][jj], tmp_jacobian_ranges_[ii][jj]);
        tmp_jacobian_ranges_[ii][jj] /= integration_element;
        ret[ii][jj] = tmp_jacobian_ranges_[ii][jj];
      }
    }
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  mutable DomainType tmp_domain_;
  mutable typename EntityType::Geometry::JacobianTransposed tmp_jacobian_transposed_;
  mutable typename EntityType::Geometry::JacobianInverseTransposed tmp_jacobian_inverse_transposed_;
  std::unique_ptr<const PdelabLFSType> lfs_;
  std::unique_ptr<const BackendType> backend_;
  mutable std::vector<RangeType> tmp_ranges_;
  mutable std::vector<JacobianRangeType> tmp_jacobian_ranges_;
}; // class PiolaTransformedDunePdelabWrapper


// New for Nedelec
/** \brief Wrapper-class for the Hcurl-conforming transformation, e.g. for member functions of the Nedelec spaces
 *
 * \note As the curl is only defined in dimension 3, the domainDim and the rangeDim have to equal 3
 *
 * \tparam PdelabSpaceType Type of function space, has to be implemented in dune-pdelab/finiteelementmap
 * \tparam EntityImp Type of Entity the transformation maps from
 * \tparam DomainFieldImp Type of the domain field
 * \tparam domainDim Dimension of the domain, has to be 3
 * \tparam RangeFieldImp Type of the range field
 * \tparam rangeDim Dimension of the range, has to be 3
 */
template< class PdelabSpaceType, class EntityImp,
          class DomainFieldImp, size_t domainDim,
          class RangeFieldImp, size_t rangeDim >
class Edges05DunePdelabWrapper< PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface< internal::Edges05DunePdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                                                       DomainFieldImp, domainDim,
                                                                       RangeFieldImp, rangeDim >,
                                     DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef Edges05DunePdelabWrapper
    < PdelabSpaceType, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > ThisType;
  typedef BaseFunctionSetInterface< internal::Edges05DunePdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                                                      DomainFieldImp, domainDim,
                                                                      RangeFieldImp, rangeDim >,
                                    DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
      BaseType;
public:
  typedef internal::Edges05DunePdelabWrapperTraits< PdelabSpaceType, EntityImp,
                                            DomainFieldImp, domainDim,
                                            RangeFieldImp, rangeDim >
      Traits;
  typedef typename Traits::BackendType   BackendType;
  typedef typename Traits::EntityType    EntityType;
private:
  typedef typename Traits::PdelabLFSType PdelabLFSType;
  typedef typename Traits::FESwitchType  FESwitchType;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using BaseType::dimDomain;

  Edges05DunePdelabWrapper(const PdelabSpaceType& space, const EntityType& ent)
    : BaseType(ent)
    , tmp_domain_(DomainFieldType(0))
    , tmp_jacobian_transposed_(DomainFieldType(0))
    , tmp_jacobian_inverse_transposed_(DomainFieldType(0))
  {
    PdelabLFSType* lfs_ptr = new PdelabLFSType(space);
    lfs_ptr->bind(this->entity());
    lfs_ = std::unique_ptr< PdelabLFSType >(lfs_ptr);
    backend_ = std::unique_ptr< BackendType >(new BackendType(FESwitchType::basis(lfs_->finiteElement())));
    tmp_ranges_ = std::vector< RangeType >(backend_->size(), RangeType(0));
    tmp_jacobian_ranges_ = std::vector< JacobianRangeType >(backend_->size(), JacobianRangeType(0));
  } // DunePdelabWrapper(...)

  Edges05DunePdelabWrapper(ThisType&& source) = default;

  Edges05DunePdelabWrapper(const ThisType& /*other*/) = delete;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  virtual size_t size() const override final
  {
    return backend_->size();
  }

  virtual size_t order() const override final
  {
    return backend_->order();
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const override final
  {
    assert(lfs_);
    assert(backend_);
    assert(ret.size() >= backend_->size());
    backend_->evaluateFunction(xx, ret);
    //evaluateFunction for edges05 already gives the value on the element and not in the reference configuration, so no mapping here
  } // ... evaluate(...)

  using BaseType::evaluate;

  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const override final
  {
    assert(lfs_);
    assert(backend_);
    assert(ret.size() >= backend_->size());
    backend_->evaluateJacobian(xx, ret);
    //again, evaluateJacobian for edge05 already gives the value on the local element and not on the reference element, so no mapping
  } // ... jacobian(...)

  using BaseType::jacobian;

private:
  mutable DomainType tmp_domain_;
  mutable typename EntityType::Geometry::JacobianTransposed        tmp_jacobian_transposed_;
  mutable typename EntityType::Geometry::JacobianInverseTransposed tmp_jacobian_inverse_transposed_;
  std::unique_ptr< const PdelabLFSType > lfs_;
  std::unique_ptr< const BackendType >   backend_;
  mutable std::vector< RangeType >         tmp_ranges_;
  mutable std::vector< JacobianRangeType > tmp_jacobian_ranges_;
}; // class Edges05DunePdelabWrapper


#else // HAVE_DUNE_PDELAB


template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabWrapper
{
  static_assert(AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};


template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class PiolaTransformedDunePdelabWrapper
{
  static_assert(AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};

//new for nedelec

template <class PdelabSpaceImp, class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp,
          size_t rangeDim, size_t rangeDimCols = 1>
class Edges05DunePdelabWrapper
{
  static_assert(AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace BaseFunctionSet
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_DUNE_PDELAB_WRAPPER_HH
