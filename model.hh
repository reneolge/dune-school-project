#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

#include "probleminterface.hh"

// DiffusionModel
// --------------

template< class FunctionSpace, class GridPart >
struct DiffusionModel
{
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef ProblemInterface< FunctionSpaceType > ProblemType ;

protected:
  enum FunctionId { rhs, bnd };
  template <FunctionId id>
  class FunctionWrapper;
public:
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhs>, GridPartType > RightHandSideType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bnd>, GridPartType > DirichletBoundaryType;

  //! constructor taking problem reference 
  DiffusionModel( const ProblemType& problem, const GridPart &gridPart )
    : problem_( problem ),
      gridPart_(gridPart),
      rhs_(problem_),
      bnd_(problem_)
  {
  }

  template< class Entity, class Point >
  void source ( const Entity &entity, 
                const Point &x,
                const RangeType &value, 
                RangeType &src ) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    RangeType m;
    problem_.m(xGlobal,m);
    for (unsigned int i=0;i<src.size();++i)
      src[i] = m[i]*value[i];
  }

  //! return the diffusive flux 
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity, 
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    // the flux is simply the identity 
    flux = gradient;
  }

  // return Fem :: Function for right hand side 
  RightHandSideType rightHandSide(  ) const 
  {
    return RightHandSideType( "right hand side", rhs_, gridPart_, 5 );  
  }
  
protected:
  template <FunctionId id>
  class FunctionWrapper : public Dune::Fem::Function< FunctionSpaceType, FunctionWrapper< id > >
  {
    const ProblemInterface<FunctionSpaceType>& impl_;
    public:   
    FunctionWrapper( const ProblemInterface<FunctionSpaceType>& impl )
    : impl_( impl ) {}
 
    //! evaluate function 
    void evaluate( const DomainType& x, RangeType& ret ) const 
    {
      if( id == rhs ) 
      {
        // call right hand side of implementation 
        impl_.f( x, ret );
      }
      else 
      {
        DUNE_THROW(Dune::NotImplemented,"FunctionId not implemented"); 
      }
    }
  };
   
  const ProblemType& problem_;
  const GridPart &gridPart_;
  FunctionWrapper<rhs> rhs_;
  FunctionWrapper<bnd> bnd_;
};

#endif // #ifndef ELLIPTC_MODEL_HH
