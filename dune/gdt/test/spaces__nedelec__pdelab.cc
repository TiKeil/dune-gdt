// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_TIMED_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_DEBUG_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#endif

#include "config.h"

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <iostream>
#include <vector>
#include <string>


#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/alugrid.hh>

#include <dune/stuff/functions/global.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/periodicview.hh>

#include <dune/gdt/spaces/nedelec/dune-pdelab-wrapper.hh>

#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/local/integrands/curlcurl.hh>
#include <dune/gdt/local/integrands/elliptic.hh>


using namespace Dune;

TEST(nedelec, alugridsimplex2d)
{

   // instantiate alugrid
    const int dim = 3;
    typedef Dune::ALUGrid< dim, dim, simplex, conforming > GridType;
    unsigned int num_macro_cubes = 4;
    Stuff::Grid::Providers::Cube< GridType > grid_provider(0.0, 1.0, num_macro_cubes);
    auto& grid = grid_provider.grid();
    typedef GridType::LeafGridView LeafGridView;
    auto leafView = grid.leafGridView();
    typedef LeafGridView::Codim< 0 >::Entity EntityType;
    typedef LeafGridView::template Codim< 0 >::Iterator ElementLeafIterator;
    ElementLeafIterator it = leafView.template begin< 0 >();
    std::cout << "hi" << std::endl;


    //Fkt addieren und multipliziere
    typedef Dune::Stuff::Functions::Expression< EntityType, double, 3, double, 1 > ExpressionFct;
    typedef Dune::Stuff::Functions::Constant< EntityType, double, 3, double , 1 > ConstFct;
    ExpressionFct myfunction("x", "x[0]+x[1]", 1);
    // ExpressionFctType a("x", "3.0", 0);

    ConstFct a(1.0);
    auto newfct = myfunction+a;
    Dune::Stuff::Functions::Product< ConstFct, ExpressionFct > newfct2(a, myfunction);
    //    typedef Dune::Stuff::Functions::Expression< EntityType, double, 3, std::complex< double >, 1> ComplexExpressionFct;
    //    ComplexExpressionFct mycfct("x", "x[0]+i*x[1]", 1);
        //Fazit: geht genau so, wie es hier steht. RangeFiled muss bei beiden Fkt natuerlich gleich sein. Geht aber nicht fuer komplexe Fkt, da Fktsauswertungen
        // bzw. Ueberpruefungen gegen unedlich etc. (ohne Betrag) ausgefuehrt werden
        // es sollten die Parameterfkt fertig an curlcurldiskretisierung ubergeben werden!
        // Problem: kann noch nicht einmal komplex-wertige Expressionfkt anlegen, da deren evaluate-Methode std::isnan und std::isinf benutzt!!!


        //local evaluation curl curl testen
    const ConstFct permeab(1.0);\
    typedef Dune::Stuff::Functions::Expression< EntityType, double, 3, double, 3 > ExpressionFctVector;
    const std::vector< std::string > expressions{"0", "0", "x[0]"};
    const std::vector< std::vector< std::string > > gradients{{"0", "0", "0"}, {"0", "0", "0"}, {"1", "0", "0"}};
    const ExpressionFctVector testfct("x", expressions, 1, "stuff.globalfunction.expression", gradients);
    const std::vector< std::string > expressions1{"0", "0", "-x[0]"};
    const std::vector< std::vector< std::string > > gradients1{{"0", "0", "0"}, {"0", "0", "0"}, {"-1", "0", "0"}};
    const ExpressionFctVector testfct1("x", expressions1, 1, "stuff.globalfunction.expression", gradients1);
    const Dune::FieldVector< double, 3 > localpoint(0.5);
    Dune::DynamicMatrix< double > ret(1, 1, 0.0);


   // newfct.visualize(leafView, "id", false);
    Dune::GDT::LocalCurlCurlIntegrand< ConstFct > localeval(permeab);
    auto permeabtuple = localeval.localFunctions(*it);
    localeval.evaluate(permeabtuple, *(testfct.local_function(*it)), *(testfct1.local_function(*it)), localpoint, ret);
    std::cout << ret <<std::endl;
    //funktioniert!!! :-)

    std::cout<< expressions[2] << std::endl;

}

TEST(nedelec, alugridsimplex3d)
{

}
