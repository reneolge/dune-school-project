#ifndef ELLIPT_FEMSCHEME_HH
#define ELLIPT_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/solver/cginverseoperator.hh>

// lagrange interpolation 
#include <dune/fem/operator/lagrangeinterpolation.hh>

/*********************************************************/

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// local includes
#include "probleminterface.hh" 

#include "model.hh"

#include "rhs.hh"
#include "elliptic.hh"

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "poisson-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};

// FemScheme 
//----------

/*******************************************************************************
 * template arguments are:
 * - GridPart: the part of the grid used to tesselate the 
 *             computational domain
 * - Model: description of the data functions and methods required for the
 *          elliptic operator (massFlux, diffusionFlux)
 *     Model::ProblemType boundary data, exact solution, 
 *                        and the type of the function space
 *******************************************************************************/
template < class Model > 
class FemScheme 
{
public:   
  //! type of the mathematical model 
  typedef Model ModelType ;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, f: \Omega -> R) 
  typedef typename ModelType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation 
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;

  /*********************************************************/

  //! define Laplace operator
  typedef EllipticOperator< DiscreteFunctionType, ModelType > EllipticOperatorType;

  FemScheme( GridPartType &gridPart, 
             const ModelType& implicitModel )
    : implicitModel_( implicitModel ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      rhs_( "rhs", discreteSpace_ ),
      // the elliptic operator (implicit)
      implicitOperator_( implicitModel_, discreteSpace_ ), 
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) )
  {
    // set all DoF to zero 
    solution_.clear();
  }

  DiscreteFunctionType &solution() 
  {
    return solution_;
  }
  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  //! setup the right hand side
  void prepare() 
  { 

    // assemble rhs
    assembleRHS ( implicitModel_.rightHandSide(), rhs_ );

  }

  //! solve the system 
  void solve ( bool assemble )
  { 
    // on-the-fly version (does not work with ISTL solvers)
    LinearInverseOperatorType solver( implicitOperator_, solverEps_, solverEps_ );
    // solve system 
    solver( rhs_, solution_ );
  }

protected:  
  const ModelType& implicitModel_;   // the mathematical model 

  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with 

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space 
  DiscreteFunctionType solution_;   // the unknown 
  DiscreteFunctionType rhs_;        // the right hand side 
  
  EllipticOperatorType implicitOperator_; // the implicit operator 

  const double solverEps_ ; // eps for linear solver 
};

#endif // end #if ELLIPT_FEMSCHEME_HH
