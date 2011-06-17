#ifndef DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
#define DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH

// dune-functionals includes
#include <dune/functionals/common/localvector.hh>

namespace Dune {

namespace Functionals {

namespace Assembler {

namespace Local {

namespace Codim0 {

template <class LocalFunctionalImp>
class Vector
{
public:
  typedef LocalFunctionalImp LocalFunctionalType;

  typedef typename LocalFunctionalType::RangeFieldType RangeFieldType;

  //! constructor
  Vector(const LocalFunctionalType& localFunctional)
    : localFunctional_(localFunctional)
  {
  }

  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

  template <class TestSpaceType, class EntityType, class VectorType,
            class LocalVectorType = Dune::Functionals::Common::LocalVector<RangeFieldType>>
  void assembleLocal(const TestSpaceType& testSpace, const EntityType& entity, VectorType& vector,
                     LocalVectorType localVector) const
  {
    // write local operator application to tmpLocalMatrix
    localFunctional_.applyLocal(testSpace.localBaseFunctionSet(entity), localVector);

    // write local matrix to global
    addToVector(testSpace, entity, localVector, vector);
  }

private:
  template <class TestSpaceType, class EntityType, class LocalVectorType, class VectorType>
  void addToVector(const TestSpaceType& testSpace, const EntityType& entity, const LocalVectorType& localVector,
                   VectorType& vector) const
  {
    for (int j = 0; j < testSpace.baseFunctionSet(entity).numBaseFunctions(); ++j) {
      const int globalJ = testSpace.mapToGlobal(entity, j);

      vector[globalJ] += localVector[j];
    }
  } // end method addToVector

  const LocalFunctionalType& localFunctional_;

}; // end class Vector

} // end namespace Codim0

} // end namespace Local

} // end namespace Assembler

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
