#include <math.h>   //cos (in myfilter)

#include <config.h>

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include header of elliptic solver
#include "femscheme.hh"
#include "poisson.hh"
#include "myfilter.hh"


// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
double algorithm ( HGridType &grid, int step )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > GridPartType;
  GridPartType gridPart(grid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, 
              HGridType::dimensionworld, 1 > FunctionSpaceType;

  // type of the mathematical model used 
  typedef DiffusionModel< FunctionSpaceType, GridPartType > ModelType;
  

  assert( problemPtr );
  ProblemType& problem = *problemPtr ;

  // implicit model for left hand side 
  ModelType implicitModel( problem, gridPart );

  // poisson solver
  typedef FemScheme< ModelType > SchemeType;
  SchemeType scheme( gridPart, implicitModel );

  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution) ; // tuple with pointers 
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  // setup the right hand side
  scheme.prepare();
  // solve once 
  scheme.solve( true );

  // write initial solve 
  dataOutput.write();

  // calculate error 
  double error = 0 ;

  {
    // calculate standard error 
    // select norm for error computation
    typedef Dune::Fem::L2Norm< GridPartType > NormType;
//    typedef Dune::Fem::H1Norm< GridPartType > NormType;									 FÜR H1-Norm diese Zeile benutzen
																							// TODO: CASE ZUM HIN UND HERSCHALTEN ZWISCHEN L2 UND H1
    NormType norm( gridPart );
    return norm.distance( gridExactSolution, scheme.solution() );							// L2-FEHLER BERECHNUNG MIT gridExactSolution
  }

  return error ;
}

template <class HGridType>
DiscreteFunctionType GetExactSolution ( HGridType &grid, int wdh_exact )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > GridPartType;
  GridPartType gridPart(grid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, 
              HGridType::dimensionworld, 1 > FunctionSpaceType;

  // type of the mathematical model used 
  typedef DiffusionModel< FunctionSpaceType, GridPartType > ModelType;
  

  assert( problemPtr );
  ProblemType& problem = *problemPtr ;

  // implicit model for left hand side 
  ModelType implicitModel( problem, gridPart );

  // poisson solver
  typedef FemScheme< ModelType > SchemeType;
  SchemeType scheme( gridPart, implicitModel );

  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution) ; // tuple with pointers 
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

for ( int step = 1; step <= wdh_exact; ++step )
  {
	  Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );								// Verfeinern um ExactSolution "genau" zu haben
  }
  // setup the right hand side
  scheme.prepare();
  // solve once 
  scheme.solve( true );

  DiscreteFunctionSpaceType gridExactSolution ( scheme.solution() );

  return gridExactSolution;
}
// mainfunktion
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line 
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files 
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file  
  Dune::Fem::Parameter::append( "../data/parameter" );

  // type of hierarchical grid 
  typedef Dune::GridSelector::GridType  HGridType ;

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static 
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser 
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance 
  grid.loadBalance();

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "poisson.level" );

  // number of global refinements to bisect grid width 
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid 
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "poisson.repeats", 0 );

  GetExactSolution ( grid , 10 )
  
  // calculate first step  
  double oldError = algorithm( grid, (repeats > 0) ? 0 : -1 );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected 
    // and all memory is adjusted correctly 
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    const double newError = algorithm( grid, step );
    const double eoc = log( oldError / newError ) / M_LN2;
    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      std::cout << "Error: " << newError << std::endl;
      std::cout << "EOC( " << step << " ) = " << eoc << std::endl;
    }
    oldError = newError;
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
