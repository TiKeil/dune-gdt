#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_BINARY_ELLIPTIC_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_BINARY_ELLIPTIC_HH

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>

// dune-helper-tools includes
#include <dune/helper-tools/function/expression.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace Evaluation {

namespace Local {

namespace Binary {

/**
  \brief  This represents the operation \f$a\nabla u \nabla v\f$.

          \f$a\f$ is a given scalar function (in this case 1) and \f$u\f$ and \f$u\f$ may be local functions, i.e.
          ansatz- and testfunctions.
  \tparam FunctionSpaceImp
          Type of the function space, where \f$f\f$, \f$u\f$ and \f$v\f$ live in.
  **/
template <class FunctionSpaceImp, class InducingFunctionImp = Dune::HelperTools::Function::Expression<
                                      typename FunctionSpaceImp::DomainFieldType, FunctionSpaceImp::DimDomain,
                                      typename FunctionSpaceImp::RangeFieldType, FunctionSpaceImp::DimRange>>
class Elliptic
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef Elliptic<FunctionSpaceType> ThisType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef InducingFunctionImp InducingFunctionType;

  //! constructor, takes the inducing function as a runtime parameter
  Elliptic(const Dune::shared_ptr<const InducingFunctionType>& inducingFunction, int order = 0)
    : inducingFunction_(inducingFunction)
    , order_(std::max(0, order))
  {
  }

  Elliptic(const Dune::ParameterTree& paramTree)
    : inducingFunction_(new InducingFunctionType(paramTree))
    , order_(paramTree.get("order", 0))
  {
  }

  //! returns the inducing function
  InducingFunctionType& inducingFunction() const
  {
    return *inducingFunction_;
  }

  unsigned int order() const
  {
    return order_;
  }

  /**
    * \brief      Evaluates \f$a(x)\nabla u(x) \nabla v(x)\f$ for a given local point \f$x\f$.
    *
    * \tparam     LocalAnsatzFunctionType
    *             Type of the local ansatz function \f$u\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalTestFunctionType
    *             Type of the local test function \f$v\f$, i.e. Dune::LocalFunction.
    * \tparam     LocalPointType
    *             Type of the local point \f$x\f$, i.e. Dune::FieldVector.
    * \param[in]  localAnsatzFunction
    *             The local function \f$u\f$.
    * \param[in]  localTestFunction
    *             The local function \f$v\f$.
    * \param[in]  localPoint
    *             The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
    *             element.
    * \return     \f$a(x)\nabla u(x) \nabla v(x)\f$
    **/
  template <class LocalAnsatzBaseFunctionSetType, class LocalTestBaseFunctionSetType, class LocalMatrixType>
  void evaluateLocal(const LocalAnsatzBaseFunctionSetType& localAnsatzBaseFunctionSet,
                     const LocalTestBaseFunctionSetType& localTestBaseFunctionSet, const DomainType& localPoint,
                     LocalMatrixType& ret) const
  {
    // get global point
    const DomainType globalPoint = localAnsatzBaseFunctionSet.entity().geometry().global(localPoint);

    // evaluate first gradient
    const unsigned int rows = localAnsatzBaseFunctionSet.size();
    std::vector<JacobianRangeType> gradientLocalAnsatzBaseFunctionSet(rows, JacobianRangeType(0.0));
    localAnsatzBaseFunctionSet.jacobian(localPoint, gradientLocalAnsatzBaseFunctionSet);

    // evaluate second gradient
    const unsigned int cols = localTestBaseFunctionSet.size();
    std::vector<JacobianRangeType> gradientLocalTestBaseFunctionSet(cols, JacobianRangeType(0.0));
    localTestBaseFunctionSet.jacobian(localPoint, gradientLocalTestBaseFunctionSet);

    // evaluate inducing function
    RangeType functionValue(0.0);
    inducingFunction_->evaluate(globalPoint, functionValue);

    // do loop over all ansatz and test basefunctions
    assert(ret.rows() == rows);
    assert(ret.cols() == cols);
    for (unsigned int i = 0; i < rows; ++i) {
      for (unsigned int j = 0; j < cols; ++j) {
        const RangeFieldType gradientProduct =
            gradientLocalAnsatzBaseFunctionSet[i][0] * gradientLocalTestBaseFunctionSet[j][0];
        ret[i][j] = functionValue * gradientProduct;
      }
    }
  } // end method evaluateLocal

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  //! copy constructor
  Elliptic(const Elliptic& other);

  const Dune::shared_ptr<const InducingFunctionType> inducingFunction_;
  unsigned int order_;
}; // end class Elliptic

} // end namespace Binary

} // end namespace Local

} // end namespace Evaluation

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_LOCAL_BINARY_ELLIPTIC_HH
