#ifndef DUNE_FUNCTIONALS_CONSTRAINTS_DIRICHLET_HH
#define DUNE_FUNCTIONALS_CONSTRAINTS_DIRICHLET_HH

#include "localdefault.hh"

namespace Dune {

namespace Functionals {

namespace Constraints {

/** @brief Constraints for zero Dirichlet values on the entire boundary domain
 *
 * This class implements constraints on the degrees of freedom on a @link
 * Subspace::Linear "linear subspace" @endlink.
 *
 * Constraints efficiently implement functionals @f$ f_i:{\cal X}_H \to
 * \mathbb{R}@f$ for each degree of freedom @f$i\in {\cal C} \subset
 * \{1,\dots,H\}@f$, where @f$ H @f$ is the number of degree of freedoms in the
 * underlying discrete function space @f$ {\cal X}_H @f$. In case of this
 * Dirichlet constraints class, the set @f$ {\cal C} @f$ includes all degrees of
 * freedom "lying" on boundary edges and @f$ f_i(u) = 0 @f$ for all @f$ i=\cal
 * C @f$.
 *
 * @note The Dirichlet constraints only make sense on a finite element space,
 * not on a discontinuous discrete function space.
 *
 * @ingroup Constraints
 *
 * @tparam DiscFuncSpace discrete function space on which constraints shall
 * be applied.
 */
template <class DiscreteFunctionSpaceImp>
class DirichletZero
{
public:
  //! Discrete function space on which the Dirichlet constraints are applied
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  //! Underlying grid part
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  //! @brief dimension of the grid part
  static const int griddim = GridPartType::GridType::dimension;

  //! @brief Return type of local() method, implementing the LocalConstraints
  //! interface
  typedef Dune::Functionals::Constraints::LocalDefault<double, 2 * griddim, 6 * griddim> LocalConstraintsType;

public:
  /** @brief Constructor for the Dirichlet constraints
   *
   *  @param space    discrete function space object on which the Dirichlet
   *                  constraints are applied
   */
  DirichletZero(const DiscreteFunctionSpaceType& space)
    : space_(space)
    , gridPart_(space.gridPart())
  {
  }


  /** @brief returns a local constraint object for the entity \a en
   *
   * @param en Entity for which the local constraints shall be compted
   *
   * @returns local constraints object (copyable). c.f.
   * Constraints::LocalDefault for more details on the return type
   * implementation.
   */
  template <class Entity>
  const LocalConstraintsType local(const Entity& entity) const
  {
    typedef typename DiscreteFunctionSpaceType::HostSpaceType DiscFuncSpace;
    typedef typename DiscFuncSpace::BaseFunctionSetType BFS;
    typedef typename DiscFuncSpace::LagrangePointSetType LPSType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;

    const int faceCodim = 1;
    typedef typename LPSType::template Codim<faceCodim>::SubEntityIteratorType FaceDofIteratorType;

    const BFS& bfs             = space_.hostSpace().baseFunctionSet(entity);
    const unsigned int numCols = bfs.numBaseFunctions();
    LocalConstraintsType lc(numCols);

    const LPSType& lps = space_.hostSpace().lagrangePointSet(entity);

    /*    // get slave dof structure (for parallel runs)
     *    SlaveDofsType &slaveDofs = this->slaveDofs();
     *    const int numSlaveDofs = slaveDofs.size();*/

    // set of local boundary dofs
    std::set<unsigned int> localBoundaryDofs;

    // loop over all intersections and find dirichlet
    // boundaries
    const IntersectionIterator endIt = gridPart_.iend(entity);
    for (IntersectionIterator it = gridPart_.ibegin(entity); it != endIt; ++it) {
      // get intersection
      const Intersection& ii = *it;

      // skip non-boundary elements
      if (!ii.boundary())
        continue;

      // get local face number of boundary intersection
      const int face = ii.indexInInside();

      // get iterator over all local dofs on this face
      FaceDofIteratorType faceIt          = lps.template beginSubEntity<faceCodim>(face);
      const FaceDofIteratorType faceEndIt = lps.template endSubEntity<faceCodim>(face);

      // iterate over face dofs and set unit row
      for (; faceIt != faceEndIt; ++faceIt) {
        const int localDof = *faceIt;

        localBoundaryDofs.insert(localDof);

        /*        // clear all other columns
         *        const int globalDof = dfSpace.mapToGlobal( entity, localDof );

         *        // cancel all slave dofs (for parallel runs)
         *        for( int i = 0; i < numSlaveDofs; ++i )
         *        {
         *          // slave dofs are canceled
         *          if( globalDof == slaveDofs[ i ] )
         *            localMatrix.set( localDof, localDof, 0 );
         *        }*/
      }
    }


    /************************************************************************
     * iterate over local boundary dof set and fill local constraint matrix *
     ************************************************************************/
    typedef std::set<unsigned int>::const_iterator LBIterator;

    LBIterator lbend     = localBoundaryDofs.end();
    unsigned int numRows = 0;
    for (LBIterator lbit = localBoundaryDofs.begin(); lbit != lbend; ++lbit, ++numRows) {
      const unsigned int localDof = *lbit;
      for (unsigned int i = 0; i < numCols; ++i) {
        if (numRows == 0)
          lc.setColumnDofs(i, space_.hostSpace().mapToGlobal(entity, i));
        lc.setLocalMatrix(numRows, i, 0.0);
      }
      lc.setLocalMatrix(numRows, localDof, 1.0);
      lc.setRowDofs(numRows, space_.hostSpace().mapToGlobal(entity, localDof));
    }
    lc.setRowDofsSize(numRows);

    return lc;
  }

private:
  const DiscreteFunctionSpaceType& space_;
  const GridPartType& gridPart_;
}; // end class Dirichlet

} // end of namespace Constraints

} // end of namespace Functionals

} // end of namespace Dune

#endif /* end of include guard: DUNE_FUNCTIONALS_CONSTRAINTS_DIRICHLET_HH */
