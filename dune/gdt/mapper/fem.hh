// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_MAPPER_FEM_HH
#define DUNE_GDT_MAPPER_FEM_HH

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/mapper/nonblockmapper.hh>
#endif

#include <dune/stuff/common/type_utils.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Mapper {

#if HAVE_DUNE_FEM


// forward, to be used in the traits and to allow for specialization
template <class FemDofMapperImp, int block_size = 1>
class FemDofWrapper;


template <class FemDofMapperImp, int block_size>
class FemDofWrapperTraits
{
public:
  typedef FemDofWrapper<FemDofMapperImp, block_size> derived_type;
  typedef Fem::NonBlockMapper<FemDofMapperImp, block_size> BackendType;
  typedef typename BackendType::ElementType EntityType;
};


template <class FemDofMapperImp>
class FemDofWrapperTraits<FemDofMapperImp, 1>
{
public:
  typedef FemDofWrapper<FemDofMapperImp, 1> derived_type;
  typedef FemDofMapperImp BackendType;
  typedef typename BackendType::ElementType EntityType;
};


template <class FemDofMapperImp, int block_size>
class FemDofWrapper : public MapperInterface<FemDofWrapperTraits<FemDofMapperImp, block_size>>
{
  typedef MapperInterface<FemDofWrapperTraits<FemDofMapperImp, block_size>> InterfaceType;

public:
  typedef FemDofWrapperTraits<FemDofMapperImp, block_size> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  explicit FemDofWrapper(FemDofMapperImp& femNonBlockMapper)
    : backend_(femNonBlockMapper)
  {
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size();
  }

  size_t numDofs(const EntityType& entity) const
  {
    return backend_.numDofs(entity);
  }

  size_t maxNumDofs() const
  {
    return backend_.maxNumDofs();
  }

private:
  class Functor
  {
  public:
    explicit Functor(Dune::DynamicVector<size_t>& globalIndices)
      : globalIndices_(globalIndices)
    {
    }

    void operator()(int localDoF, int globalDoF)
    {
      assert(localDoF < int(globalIndices_.size()));
      globalIndices_[localDoF] = globalDoF;
    }

  private:
    Dune::DynamicVector<size_t>& globalIndices_;
  };

public:
  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    // some checks
    const size_t numLocalDofs = numDofs(entity);
    if (ret.size() < numLocalDofs)
      ret.resize(numLocalDofs);
    // compute
    Functor functor(ret);
    backend_.mapEach(entity, functor);
  }

  using InterfaceType::globalIndices;

  /**
   *  \attention  This method is implemented using globalIndices() and thus not optimal!
   */
  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    const size_t numLocalDofs = numDofs(entity);
    assert(localIndex < numLocalDofs);
    Dune::DynamicVector<size_t> tmpGlobalIndices(numLocalDofs);
    globalIndices(entity, tmpGlobalIndices);
    return tmpGlobalIndices[localIndex];
  }

private:
  const BackendType backend_;
}; // class FemDofWrapper


template <class FemDofMapperImp>
class FemDofWrapper<FemDofMapperImp, 1> : public MapperInterface<FemDofWrapperTraits<FemDofMapperImp, 1>>
{
  typedef MapperInterface<FemDofWrapperTraits<FemDofMapperImp, 1>> InterfaceType;

public:
  typedef FemDofWrapperTraits<FemDofMapperImp, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  explicit FemDofWrapper(const BackendType& femMapper)
    : backend_(femMapper)
  {
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size();
  }

  size_t numDofs(const EntityType& entity) const
  {
    return backend_.numDofs(entity);
  }

  size_t maxNumDofs() const
  {
    return backend_.maxNumDofs();
  }

private:
  class Functor
  {
  public:
    explicit Functor(Dune::DynamicVector<size_t>& globalIndices)
      : globalIndices_(globalIndices)
    {
    }

    void operator()(int localDoF, int globalDoF)
    {
      assert(localDoF < int(globalIndices_.size()));
      globalIndices_[localDoF] = globalDoF;
    }

  private:
    Dune::DynamicVector<size_t>& globalIndices_;
  };

public:
  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    // some checks
    const size_t numLocalDofs = numDofs(entity);
    if (ret.size() < numLocalDofs)
      ret.resize(numLocalDofs);
    // compute
    Functor functor(ret);
    backend_.mapEach(entity, functor);
  }

  using InterfaceType::globalIndices;

  /**
   *  \attention  This method is implemented using globalIndices() and thus not optimal!
   */
  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    const size_t numLocalDofs = numDofs(entity);
    assert(localIndex < numLocalDofs);
    Dune::DynamicVector<size_t> tmpGlobalIndices(numLocalDofs);
    globalIndices(entity, tmpGlobalIndices);
    return tmpGlobalIndices[localIndex];
  }

private:
  const BackendType& backend_;
}; // class FemDofWrapper< ..., 1 >


#else // HAVE_DUNE_FEM


template <class FemDofMapperImp>
class FemDofWrapper
{
  static_assert(Dune::AlwaysFalse<FemDofMapperImp>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MAPPER_FEM_HH
