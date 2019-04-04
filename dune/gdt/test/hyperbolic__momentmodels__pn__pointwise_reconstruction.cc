// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

#define USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR 0
#include <dune/gdt/test/hyperbolic_momentmodels_pn.hh>
#include <dune/gdt/test/pn-discretization.hh>

TYPED_TEST_CASE(HyperbolicPnTest, YaspGridTestCasesWithReconstruction);
TYPED_TEST(HyperbolicPnTest, check_with_pointwise_linear_reconstruction)
{
  this->run();
}
