#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#define HAVE_DUNE_DETAILED_DISCRETIZATIONS 1

#include <iostream>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>
#include <dune/common/dynvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/discretefunction/projection/dirichlet.hh>

#include <dune/detailed/discretizations/space/continuouslagrange/fem.hh>
#include <dune/detailed/discretizations/la/containerfactory/eigen.hh>
#include <dune/detailed/discretizations/evaluation/elliptic.hh>
#include <dune/detailed/discretizations/evaluation/product.hh>
#include <dune/detailed/discretizations/localoperator/codim0.hh>
#include <dune/detailed/discretizations/localfunctional/codim0.hh>
#include <dune/detailed/discretizations/assembler/local/codim0.hh>
#include <dune/detailed/discretizations/space/constraints.hh>
#include <dune/detailed/discretizations/assembler/system.hh>
#include <dune/detailed/discretizations/discretefunction/default.hh>


const std::string id = "elliptic.continuousgalerkin";

#ifndef POLORDER
const int polOrder = 1;
#else
const int polOrder = POLORDER;
#endif
dune_static_assert((polOrder > 0), "ERROR: polOrder hast to be positive!");

using namespace Dune::Detailed::Discretizations;


/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(const std::string& filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << std::endl;
    file << "grid = "
         << "gridprovider.cube" << std::endl;
    file << "boundaryinfo = "
         << "boundaryinfo.alldirichlet" << std::endl;
    file << "[gridprovider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = [12; 12; 12]" << std::endl;
    file << "[boundaryinfo.idbased]" << std::endl;
    file << "dirichlet = [1; 2; 3]" << std::endl;
    file << "neumann = [4]" << std::endl;
    file << "[diffusion]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 1.0; 1.0]" << std::endl;
    file << "[force]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 1.0; 1.0]" << std::endl;
    file << "[dirichlet]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [0.1*x[0]; 0.0; 0.0]" << std::endl;
    file << "[neumann]" << std::endl;
    file << "order = 0" << std::endl;
    file << "variable = x" << std::endl;
    file << "expression = [1.0; 0.0; 0.0]" << std::endl;
    file << "[solver]" << std::endl;
    file << "type = bicgstab.ilut" << std::endl;
    file << "maxIter = 5000" << std::endl;
    file << "precision = 1e-12" << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()


int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::Fem::MPIManager::initialize(argc, argv);

    // parameter
    const std::string paramFilename = id + ".description";
    ensureParamFile(paramFilename);
    Dune::Stuff::Common::ExtendedParameterTree description(argc, argv, paramFilename);
    description.assertSub(id);

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO | Dune::Stuff::Common::LOG_CONSOLE);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();

    // timer
    Dune::Timer timer;

    info << "setting up grid:" << std::endl;
    typedef Dune::Stuff::GridProviderInterface<> GridProviderType;
    const GridProviderType* gridProvider =
        Dune::Stuff::GridProviders<>::create(description.get(id + ".grid", "gridprovider.cube"), description);
    typedef GridProviderType::GridType GridType;
    const std::shared_ptr<const GridType> grid = gridProvider->grid();
    typedef Dune::grid::Part::Leaf::Const<GridType> GridPartType;
    typedef typename GridPartType::GridViewType GridViewType;
    const GridPartType gridPart(*grid);
    typedef typename Dune::Stuff::GridboundaryInterface<typename GridPartType::GridViewType> BoundaryInfoType;
    const std::shared_ptr<const BoundaryInfoType> boundaryInfo(
        Dune::Stuff::Gridboundaries<typename GridPartType::GridViewType>::create(
            description.get(id + ".boundaryinfo", "boundaryinfo.alldirichlet"), description));
    info << "  took " << timer.elapsed() << " sec, has " << grid->size(0) << " entities" << std::endl;
    info << "visualizing grid... " << std::flush;
    timer.reset();
    gridProvider->visualize(description.get(id + ".boundaryinfo", "boundaryinfo.alldirichlet"),
                            description,
                            description.get(id + ".filename", id) + ".grid");
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

    const unsigned int DUNE_UNUSED(dimDomain) = GridProviderType::dim;
    const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef typename GridType::ctype DomainFieldType;
    typedef double RangeFieldType;

    typedef Dune::Stuff::FunctionExpression<DomainFieldType, dimDomain, RangeFieldType, dimRange>
        ExpressionFunctionType;
    const std::shared_ptr<const ExpressionFunctionType> diffusion(
        ExpressionFunctionType::create(description.sub("diffusion")));
    const std::shared_ptr<const ExpressionFunctionType> force(ExpressionFunctionType::create(description.sub("force")));
    const std::shared_ptr<const ExpressionFunctionType> dirichlet(
        ExpressionFunctionType::create(description.sub("dirichlet")));
    const std::shared_ptr<const ExpressionFunctionType> neumann(
        ExpressionFunctionType::create(description.sub("neumann")));

    // info << "initializing discrete function spaces... " << std::flush;
    // timer.reset();
    typedef ContinuousLagrangeSpace::FemWrapper<GridPartType, polOrder, RangeFieldType, dimRange> SpaceType;
    const SpaceType space(gridPart);

    // left hand side
    // * elliptic diffusion operator
    typedef LocalOperatorCodim0Integral<Evaluation::Elliptic, ExpressionFunctionType> EllipticOperatorType;
    const EllipticOperatorType diffusionOperator(*diffusion);
    // * right hand side
    //   * L2 force functional
    typedef LocalFunctionalCodim0Integral<Evaluation::Product, ExpressionFunctionType> L2FunctionalType;
    const L2FunctionalType forceFunctional(*force);
    //   * L2 neumann functional
    const L2FunctionalType neumannFunctional(*neumann);

    info << "initializing matrix (of size " << space.mapper().size() << "x" << space.mapper().size()
         << ") and vectors... " << std::flush;
    // timer.reset();
    typedef ContainerFactoryEigen<RangeFieldType> ContainerFactory;
    typedef ContainerFactory::RowMajorSparseMatrixType MatrixType;
    typedef ContainerFactory::DenseVectorType VectorType;
    std::shared_ptr<MatrixType> systemMatrix(ContainerFactory::createRowMajorSparseMatrix(space, space));
    std::shared_ptr<VectorType> forceVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> dirichletVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> neumannVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> rhsVector(ContainerFactory::createDenseVector(space));
    std::shared_ptr<VectorType> solutionVector(ContainerFactory::createDenseVector(space));
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "assembing system... " << std::flush;
    timer.reset();
    // * dirichlet boundary values
    typedef DiscreteFunctionDefault<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType dirichletProjection(space, dirichletVector, "dirichlet");
    Dune::Stuff::DiscreteFunction::project(*boundaryInfo, *dirichlet, dirichletProjection);

    // * local matrix assembler
    typedef LocalAssemblerCodim0Matrix<EllipticOperatorType> LocalEllipticOperatorMatrixAssemblerType;
    const LocalEllipticOperatorMatrixAssemblerType diffusionMatrixAssembler(diffusionOperator);
    // * local vector assemblers
    typedef LocalAssemblerCodim0Vector<L2FunctionalType> LocalL2FunctionalVectorAssemblerType;
    //   * force vector
    const LocalL2FunctionalVectorAssemblerType forceVectorAssembler(forceFunctional);
    //   * neumann vector
    const LocalL2FunctionalVectorAssemblerType neumannVectorAssembler(neumannFunctional);
    // * system assembler
    typedef SystemAssembler<SpaceType, SpaceType> SystemAssemblerType;
    SystemAssemblerType systemAssembler(space);
    systemAssembler.addLocalMatrixAssembler(diffusionMatrixAssembler, *systemMatrix);
    systemAssembler.addLocalVectorAssembler(forceVectorAssembler, *forceVector);
    systemAssembler.addLocalVectorAssembler(neumannVectorAssembler, *neumannVector);
    systemAssembler.assemble();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "applying constraints... " << std::flush;
    timer.reset();
    Constraints::Dirichlet<GridViewType, RangeFieldType> dirichletConstraints(
        *boundaryInfo, space.mapper().maxNumDofs(), space.mapper().maxNumDofs());
    rhsVector->backend() =
        forceVector->backend() + neumannVector->backend() - systemMatrix->backend() * dirichletVector->backend();
    systemAssembler.applyConstraints(dirichletConstraints, *systemMatrix, *rhsVector);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    info << "solving linear system (of size " << systemMatrix->rows() << "x" << systemMatrix->cols() << ")"
         << std::endl;
    const std::string solverType     = description.get("solver.type", "bicgstab");
    const unsigned int solverMaxIter = description.get("solver.maxIter", 5000);
    const double solverPrecision     = description.get("solver.precision", 1e-12);
    info << "  using '" << solverType << "'... " << std::flush;
    timer.reset();
    typedef typename Dune::Stuff::LA::Solver::Interface<MatrixType, VectorType> SolverType;
    std::shared_ptr<SolverType> solver(Dune::Stuff::LA::Solver::create<MatrixType, VectorType>(solverType));
    const unsigned int failure =
        solver->apply(*systemMatrix, *rhsVector, *solutionVector, solverMaxIter, solverPrecision);
    if (failure)
      DUNE_THROW(Dune::MathError, "\nERROR: linear solver '" << solverType << "' reported a problem!");
    if (solutionVector->size() != space.mapper().size())
      DUNE_THROW(Dune::MathError,
                 "\nERROR: linear solver '" << solverType << "' produced a solution of wrong size (is "
                                            << solutionVector->size()
                                            << ", should be "
                                            << space.mapper().size()
                                            << ")!");
    solutionVector->backend() += dirichletVector->backend();
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    const std::string solutionFilename = description.get(id + ".filename", id) + ".solution";
    const std::string solutionName     = id + ".solution";
    info << "writing solution to '" << solutionFilename;
    if (dimDomain == 1)
      info << ".vtp";
    else
      info << ".vtu";
    info << "'... " << std::flush;
    timer.reset();
    typedef DiscreteFunctionDefaultConst<SpaceType, VectorType> ConstDiscreteFunctionType;
    const std::shared_ptr<const ConstDiscreteFunctionType> solution(
        new ConstDiscreteFunctionType(space, solutionVector, solutionName));
    typedef Dune::VTKWriter<GridViewType> VTKWriterType;
    VTKWriterType vtkWriter(gridPart.gridView());
    vtkWriter.addVertexData(solution);
    vtkWriter.write(solutionFilename);
    info << "done (took " << timer.elapsed() << " sec)" << std::endl;

    // done
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
