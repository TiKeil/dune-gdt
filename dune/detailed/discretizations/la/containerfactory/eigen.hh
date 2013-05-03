#ifndef DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH

//#ifdef HAVE_EIGEN

#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/detailed/discretizations/space/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


template <class ElementType>
class ContainerFactoryEigen
{
public:
  typedef Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementType> RowMajorSparseMatrixType;
  typedef Dune::Stuff::LA::Container::EigenDenseMatrix<ElementType> DenseMatrixType;
  typedef Dune::Stuff::LA::Container::EigenDenseVector<ElementType> DenseVectorType;
  typedef Dune::Stuff::LA::Container::SparsityPatternDefault PatternType;

  template <class T, class A>
  static RowMajorSparseMatrixType* createRowMajorSparseMatrix(const SpaceInterface<T>& testSpace,
                                                              const SpaceInterface<A>& ansatzSpace)
  {
    const std::shared_ptr<const PatternType> pattern(testSpace.computePattern(ansatzSpace));
    return createRowMajorSparseMatrix(testSpace, ansatzSpace, *pattern);
  } // static ... createRowMajorSparseMatrix(...)

  template <class T, class A>
  static RowMajorSparseMatrixType* createRowMajorSparseMatrix(const SpaceInterface<T>& testSpace,
                                                              const SpaceInterface<A>& ansatzSpace,
                                                              const PatternType& pattern)
  {
    return new RowMajorSparseMatrixType(testSpace.mapper().size(), ansatzSpace.mapper().size(), pattern);
  } // static ... createRowMajorSparseMatrix(...)

  template <class T, class A>
  static DenseMatrixType createDenseMatrix(const SpaceInterface<T>& testSpace, const SpaceInterface<A>& ansatzSpace)
  {
    return new DenseMatrixType(testSpace.mapper().size(), ansatzSpace.mapper().size());
  } // static ... createDenseMatrix(...)

  template <class S>
  static DenseVectorType* createDenseVector(const SpaceInterface<S>& space)
  {
    return new DenseVectorType(space.mapper().size());
  } // static ... createDenseVector(const SpaceType& space)
}; // class ContainerFactoryEigen


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

//#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH
