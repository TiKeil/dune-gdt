/**
  \file   main.cc
  \brief  Main file fir the finite element example.
  **/
// disable warnings about problems in dune headers
#include <dune/fem-tools/header/disablewarnings.hh>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

// dune-grid includes
#include <dune/grid/utility/gridtype.hh>

// dune-istl includes
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/bvector.hh>

// dune-fem includes
#include <dune/fem/misc/mpimanager.hh>
//#include <dune/fem/gridpart/gridpart.hh>
//#include <dune/fem/storage/vector.hh>
//#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// reenable warnings
#include <dune/fem-tools/header/enablewarnings.hh>

// dune-functionals includes
#include <dune/functionals/discretefunctionspace/finiteelement.hh>
#include <dune/functionals/discretefunctionspace/subspace/linear.hh>
#include <dune/functionals/discretefunctionspace/subspace/affine.hh>
//#include <dune/fem/localoperation/interface.hh>
//#include <dune/fem/localoperation/integrator.hh>
//#include <dune/fem/subspace/subspaces.hh>
//#include <dune/fem/operator/linear.hh>
//#include <dune/fem/functional/finiteelement.hh>
//#include <dune/fem/container/factory.hh>
//#include <dune/fem/solver/femassembler.hh>

// dune-fem-tools includes
//#include <dune/fem-tools/function/functiontools.hh>
//#include <dune/fem-tools/space/projection.hh>

using namespace Dune::Functionals;

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif

///**
//  \brief  This represents the operation \f$fv\f$.

//          \f$f\f$ is a given right hand side (in this case 1) and \f$v\f$ may be a local function, i.e. a
//          testfunction.
//  \tparam FunctionSpaceImp
//          Type of the function space, where \f$f\f$ and \f$v\f$ live in.
//  **/
// template< class FunctionSpaceImp >
// class ProductOperation
//  : public Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
//{
// public:

//  typedef FunctionSpaceImp
//    FunctionSpaceType;

//  typedef typename FunctionSpaceType::RangeFieldType
//    RangeFieldType;

//  typedef typename FunctionSpaceType::RangeType
//    RangeType;

//  /**
//    \brief      Evaluates \f$f(x)v(x)\f$ for a given local point \f$x\f$.

//    \tparam     LocalTestFunctionType
//                Type of the local function \f$v\f$, i.e. Dune::LocalFunction.
//    \tparam     LocalPointType
//                Type of the local point \f$x\f$, i.e. Dune::FieldVector.
//    \param[in]  localTestFunction
//                The local function \f$v\f$.
//    \param[in]  localPoint
//                The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
//                element.
//    \return     \f$f(x)v(x)\f$
//    **/
//  template< class LocalTestFunctionType, class LocalPointType >
//  RangeFieldType evaluate(  const LocalTestFunctionType& localTestFunction,
//                            const LocalPointType& localPoint ) const
//  {
//    // init return value
//    RangeFieldType ret = 0.0;

//    // evaluate local function
//    RangeType localTestFunctionValue( 0.0 );
//    localTestFunction.evaluate( localPoint, localTestFunctionValue );

//    // 1.0 * v(x)
//    ret = 1.0 * localTestFunctionValue;

//    // return
//    return ret;
//  }

//}; // end class ProductOperation


///**
//  \brief  This represents the operation \f$a\nabla u \nabla v\f$.

//          \f$a\f$ is a given scalar function (in this case 1) and \f$u\f$ and \f$u\f$ may be local functions, i.e.
//          ansatz- and testfunctions.
//  \tparam FunctionSpaceImp
//          Type of the function space, where \f$f\f$, \f$u\f$ and \f$v\f$ live in.
//  **/
// template< class FunctionSpaceImp >
// class EllipticOperation
//  : public Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
//{
// public:

//  typedef FunctionSpaceImp
//    FunctionSpaceType;

//  typedef Dune::Functionals::LocalOperation::Interface< FunctionSpaceImp >
//    BaseType;

//  typedef typename FunctionSpaceType::RangeFieldType
//    RangeFieldType;

//  typedef typename FunctionSpaceType::RangeType
//    RangeType;

//  typedef typename FunctionSpaceType::JacobianRangeType
//    JacobianRangeType;

//  /**
//    * \brief      Evaluates \f$a(x)\nabla u(x) \nabla v(x)\f$ for a given local point \f$x\f$.
//    *
//    * \tparam     LocalAnsatzFunctionType
//    *             Type of the local ansatz function \f$u\f$, i.e. Dune::LocalFunction.
//    * \tparam     LocalTestFunctionType
//    *             Type of the local test function \f$v\f$, i.e. Dune::LocalFunction.
//    * \tparam     LocalPointType
//    *             Type of the local point \f$x\f$, i.e. Dune::FieldVector.
//    * \param[in]  localAnsatzFunction
//    *             The local function \f$u\f$.
//    * \param[in]  localTestFunction
//    *             The local function \f$v\f$.
//    * \param[in]  localPoint
//    *             The local point \f$x\f$. This point is local in the sense, that this is a point on a reference
//    *             element.
//    * \return     \f$a(x)\nabla u(x) \nabla v(x)\f$
//    **/
//  template< class LocalAnsatzFunctionType, class LocalTestFunctionType, class LocalPointType >
//  RangeFieldType evaluate(  const LocalAnsatzFunctionType& localAnsatzFunction,
//                            const LocalTestFunctionType& localTestFunction,
//                            const LocalPointType& localPoint ) const
//  {
//    // init return value
//    RangeFieldType ret = 0.0;

//    // evaluate first gradient
//    JacobianRangeType gradientLocalAnsatzFunction( 0.0 );
//    localAnsatzFunction.jacobian( localPoint, gradientLocalAnsatzFunction );

//    // evaluate second gradient
//    JacobianRangeType gradientLocalTestFunction( 0.0 );
//    localTestFunction.jacobian( localPoint, gradientLocalTestFunction );

//    const RangeFieldType product = gradientLocalAnsatzFunction[0] * gradientLocalTestFunction[0];

//    // 1.0 * \gradient u(x) \gradient v(x)
//    ret = 1.0 * product;

//    // return
//    return ret;
//  }

//}; // end class EllipticOperation

int main(int argc, char** argv)
{
  try {

    // MPI manager
    Dune::MPIManager::initialize(argc, argv);

    // grid
    static const unsigned int dimRange = 1;

    typedef Dune::GridSelector::GridType GridType;

    typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;

    Dune::GridPtr<GridType> gridPtr("../macrogrids/unitcube2.dgf");

    GridPartType gridPart(*gridPtr);

    // function spaces and functions
    typedef Dune::FunctionSpace<double, double, GridType::dimension, dimRange> FunctionSpaceType;

    //    typedef Dune::Function< double, double >
    //      FunctionType;


    //    // local operations
    //    typedef ProductOperation< FunctionSpaceType >
    //      ProductOperationType;

    //    ProductOperationType productOperation;

    //    typedef EllipticOperation< FunctionSpaceType >
    //      EllipticOperationType;

    //    EllipticOperationType ellipticOperation;


    //    // integration
    //    typedef LocalOperation::Integrator::Codim0< FunctionSpaceType, ProductOperationType >
    //      ProductIntegratorType;

    //    ProductIntegratorType productIntegrator( productOperation );

    //    typedef LocalOperation::Integrator::Codim0< FunctionSpaceType, EllipticOperationType >
    //      EllipticIntegratorType;

    //    EllipticIntegratorType ellipticIntegrator( ellipticOperation );


    // discrete function space
    typedef DiscreteFunctionSpace::ContinuousFiniteElement<FunctionSpaceType, GridPartType, polOrder> DiscreteH1Type;

    const DiscreteH1Type discreteH1(gridPart);

    typedef DiscreteFunctionSpace::Subspace::Linear::DirichletZero<DiscreteH1Type> DiscreteH10Type;

    const DiscreteH10Type discreteH10(discreteH1);

    typedef DiscreteFunctionSpace::Subspace::Affine::Dirichlet<DiscreteH10Type> DiscreteH1GType;

    const DiscreteH1GType discreteH1G(discreteH10, "[1.0;1.0;1.0]");


    //    // operator and functional
    //    typedef Operator::Linear< EllipticIntegratorType, DiscreteH1Type >
    //      FEMellipticOperatorType;

    //    FEMellipticOperatorType femEllipticOperator( ellipticIntegrator, discreteH1 );

    //    typedef Functional::FiniteElementLOP< DiscreteH1Type, ProductIntegratorType >
    //      FEMrhsFunctionalType;

    //    FEMrhsFunctionalType femRhsFunctional( discreteH1, productIntegrator );


    //    // matrix, rhs and solution storage
    //    typedef Dune::FieldMatrix< double, dimRange, dimRange >
    //      FieldMatrixType;

    //    typedef Container::MatrixFactory< Dune::BCRSMatrix< FieldMatrixType > >
    //      MatrixFactoryType;

    //    typedef MatrixFactoryType::ContainerType
    //      MatrixContainer;

    //    typedef MatrixFactoryType::AutoPtrType
    //      MatrixContainerPtr;

    //    typedef Container::VectorFactory< Dune::BlockVector< Dune::FieldVector< double, 1 > > >
    //      VectorFactoryType;

    //    typedef VectorFactoryType::ContainerType
    //      VectorContainer;

    //    typedef VectorFactoryType::AutoPtrType
    //      VectorContainerPtr;

    //    MatrixContainerPtr A  = MatrixFactoryType::create( discreteH10 );
    //    VectorContainerPtr F  = VectorFactoryType::create( discreteH10 );
    //    VectorContainerPtr u0 = VectorFactoryType::create( discreteH10 );


    //    // assembler
    //    typedef Assembler::FiniteElement< MatrixContainer, VectorContainer >
    //      Assembler;

    //    Assembler::assembleMatrix( femEllipticOperator, *A );
    //    Assembler::applyMatrixConstraints( discreteH10, *A );

    //    Assembler::assembleVector( femRhsFunctional, *F );
    //    Assembler::applyVectorConstraints( discreteH10, *F );


    //    // preconditioner and solver
    //    typedef Dune::MatrixAdapter< MatrixContainer, VectorContainer, VectorContainer >
    //      MatrixAdapterType;

    //    MatrixAdapterType op( *A );

    //    typedef Dune::SeqILU0< MatrixContainer, VectorContainer, VectorContainer, 1 >
    //      SeqILU0Type;

    //    SeqILU0Type prec( *A, 1.0 );

    //    typedef Dune::CGSolver< VectorContainer >
    //      CG;

    //    CG cg( op, prec, 1e-4, 100, 2 );

    //    Dune::InverseOperatorResult res;

    //    // u_0 = A^(-1) ( F - G )
    //    cg.apply( *u0, *F, res );


    //    // postprocessing
    //    DiscreteFunctionType solution = Dune::FemTools::discreteFunctionFactory< DiscreteFunctionType >( discreteH1,
    //    *u0 );
    //    Dune::FemTools::writeDiscreteFunctionToVTK( solution, "solution" );

    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
