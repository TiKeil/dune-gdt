#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_CODIM0_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_CODIM0_HH

#include <vector>
#include <utility>

#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/vector.hh>

#include "../basefunctionset/interface.hh"
#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace LocalFunctional {


// forward, to be used in the traits
template <class UnaryEvaluationImp>
class Codim0Integral;


template <class UnaryEvaluationImp>
class Codim0IntegralTraits
{
public:
  typedef Codim0Integral<UnaryEvaluationImp> derived_type;
  typedef LocalEvaluation::Codim0Interface<typename UnaryEvaluationImp::Traits, 1> UnaryEvaluationType;
};


template <class UnaryEvaluationImp>
class Codim0Integral : public LocalFunctional::Codim0Interface<Codim0IntegralTraits<UnaryEvaluationImp>>
{
public:
  typedef Codim0IntegralTraits<UnaryEvaluationImp> Traits;
  typedef typename Traits::UnaryEvaluationType UnaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  Codim0Integral(const UnaryEvaluationImp evaluation)
    : evaluation_(evaluation)
  {
  }

  template <class... Args>
  Codim0Integral(Args&&... args)
    : evaluation_(std::forward<Args>(args)...)
  {
  }

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  /**
   *  \brief      Applies the local functional.
   *  \tparam T   Traits of the BaseFunctionSetInterface implementation, representing the type of the testBase
   *  \attention  ret is assumed to be zero!
   */
  template <class T, class D, int d, class R, int r, int rC>
  void apply(const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase, Dune::DynamicVector<R>& ret,
             std::vector<Dune::DynamicVector<R>>& tmpLocalVectors) const
  {
    // local inducing function
    const auto& entity        = testBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    typedef Dune::QuadratureRules<D, d> VolumeQuadratureRules;
    typedef Dune::QuadratureRule<D, d> VolumeQuadratureType;
    const int quadratureOrder = evaluation().order(localFunctions, testBase);
    assert(quadratureOrder >= 0 && "Not implemented for negative integration orders!");
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), 2 * quadratureOrder + 1);
    // check vector and tmp storage
    const size_t size = testBase.size();
    assert(ret.size() >= size);
    assert(tmpLocalVectors.size() >= numTmpObjectsRequired_);
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector<D, d> x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight  = quadPointIt->weight();
      // clear tmp vector
      Dune::DynamicVector<R>& localVector = tmpLocalVectors[0];
      Dune::Stuff::Common::clear(localVector);
      // evaluate the local operation
      evaluation().evaluate(localFunctions, testBase, x, localVector);
      // compute integral
      for (size_t ii = 0; ii < size; ++ii) {
        ret[ii] += localVector[ii] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const UnaryEvaluationType& evaluation() const
  {
    return static_cast<const UnaryEvaluationType&>(evaluation_);
  }

  const UnaryEvaluationImp evaluation_;
}; // class Codim0Integral


} // namespace LocalFunctional
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_CODIM0_HH
