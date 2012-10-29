#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH

// system
#include <vector>

// dune-stuff
#include <dune/stuff/common/vector.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Local {

namespace Codim0 {

template <class LocalFunctionalImp>
class Vector
{
public:
  typedef LocalFunctionalImp LocalFunctionalType;

  typedef Vector<LocalFunctionalType> ThisType;

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

  std::vector<unsigned int> numTmpObjectsRequired() const
  {
    std::vector<unsigned int> ret(2, 0);
    // we require 1 tmp vector in this local assembler
    ret[0] = 1;
    // the functional itself requires that much local matrices
    ret[1] = localFunctional_.numTmpObjectsRequired();
    return ret;
  }

  template <class TestSpaceType, class EntityType, class VectorType, class LocalVectorType>
  void assembleLocal(const TestSpaceType& testSpace, const EntityType& entity, VectorType& vector,
                     std::vector<std::vector<LocalVectorType>>& tmpLocalVectorsContainer) const
  {
    // get the basefunctionset
    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType LocalTestBaseFunctionSetType;

    const LocalTestBaseFunctionSetType localTestBaseFunctionSet = testSpace.baseFunctionSet().local(entity);

    // check tmp local vectors
    assert(tmpLocalVectorsContainer.size() > 1);
    std::vector<LocalVectorType>& tmpLocalVectors = tmpLocalVectorsContainer[0];
    if (tmpLocalVectors.size() < 1) {
      tmpLocalVectors.resize(1, LocalVectorType(testSpace.map().maxLocalSize(), RangeFieldType(0.0)));
    }

    // clear target vector
    Dune::Stuff::Common::clear(tmpLocalVectors[0]);

    // write local functional application to tmpLocalVector
    localFunctional_.applyLocal(localTestBaseFunctionSet, tmpLocalVectors[0], tmpLocalVectorsContainer[1]);

    // write local vector to global
    addToVector(testSpace, entity, tmpLocalVectors[0], vector);
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  //! copy constructor
  Vector(const ThisType&);

  template <class TestSpaceType, class EntityType, class LocalVectorType, class VectorType>
  void addToVector(const TestSpaceType& testSpace, const EntityType& entity, const LocalVectorType& localVector,
                   VectorType& vector) const
  {
    for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
      const unsigned int globalJ = testSpace.map().toGlobal(entity, j);

      vector.add(globalJ, localVector[j]);
    }
  } // end method addToVector

  const LocalFunctionalType& localFunctional_;

}; // end class Vector

} // end namespace Codim0

} // end namespace Local

} // end namespace Assembler

} // namespace Discretization

} // end namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM0_VECTOR_HH
