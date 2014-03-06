#ifndef ELLIPTIC_HH
#define ELLIPTIC_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <myfilter.hh>

// EllipticOperator
// ----------------

template< class DiscreteFunction, class Model >
struct EllipticOperator
: public virtual Dune::Fem::Operator< DiscreteFunction >
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;		// FunctionSpace geändert in DiscreteFunctionSpaceType

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType; 

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

public:
  //! constructor 
  EllipticOperator ( const ModelType &model, const DiscreteFunctionSpaceType &space )
  : model_( model )
  {}
      
  // prepare the solution vector 

  //! application operator 
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

protected:
  const ModelType &model () const { return model_; }

private:
  ModelType model_;
};

// Implementation of EllipticOperator
// ----------------------------------

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{

  w.clear();
   const DiscreteFunctionSpaceType &dfSpace = u.space();


  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    LocalFunctionType u_local = u.localFunction( entity);
    LocalFunctionType w_local = w.localFunction( entity );
    
    // obtain quadrature order 
    const int quadOrder = u_local.order() + w_local.order();
    
    QuadratureType quadrature( entity, quadOrder);
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
	  typedef School::Filter< HostGridPartType > FilterType;
      FilterType filter( hostGridPart );
  
	  //switch between omega
	  const bool enthalten = filter.contains(entity);
	  switch ( enthalten )
	  {
	    case 0: 																					// ERSTER FALL (innerhalb)
	      // obtain quadrature point
		  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

		  // evaluate f
		  RangeType u_x, newu_x ;
		  u_local.evaluate( quadrature[ pt ], u_x );
		  

		  //evaluate source and diffusive flux
		  model_.source(entity, x , u_x , newu_x );
		  
		  
			DomainType y = x;
			y /= x.two_norm();

			JacobianRangeType Ju;
			uJacobian( y, Ju );
			Ju.mv( y, ret );
			ret *= dimDomain - 1;
    
		   	RangeType uVal;
			u( y, uVal );
			ret += uVal;
		  }     // TODO ENDE IM FOLGENDEN ANGENOMMEN Hu IST HESSEMATRIX
		
		
		  RangeType IN(Ju[0,0]+100*Ju[1,1]);
		    
		  // BRAUCHEN WIR NICHT, DA DIREKT IN CODE GESCHRIEBENmodel_.diffusiveFlux(entity, x , u_x , Du_x , newDu_x );
		  
		  // multiply by quadrature weight
		  //newu_x *= quadrature.weight( pt ) * geometry.integrationElement( x );
		  //newDu_x *= quadrature.weight( pt ) * geometry.integrationElement( x );
		  newIN *= quadrature.weight( pt ) * geometry.integrationElement( x );
		  
		  // add f * phi_i to rhsLocal[ i ]
		  w_local.axpy( quadrature[ pt ], newIN );
	      break ;
	      
	    case 1: 																					// ZWEITER FALL (außerhalb)
	    
	      // obtain quadrature point
		  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

		  // evaluate f
		  RangeType u_x, newu_x ;
		  u_local.evaluate( quadrature[ pt ], u_x );
		  

		  //evaluate source and diffusive flux
		  model_.source(entity, x , u_x , newu_x );
		  
		  
		  //  TODO: HESSEMATRIX MUSS NOCH BERECHNET WERDEN! WAS PASSIERT HIER
			DomainType y = x;
			y /= x.two_norm();

			JacobianRangeType Ju;
			uJacobian( y, Ju );
			Ju.mv( y, ret );
			ret *= dimDomain - 1;
    
		    HessianRangeType Hu;
			hessian( y, Hu );
			for( int i = 0; i < dimDomain; ++i )
			{
			  DomainType ds = y;
			  ds *= -y[ i ];
			  ds[ i ] += 1;

			  for( int j = 0; j < dimRange; ++j )
			  {
				DomainType Hds;
				Hu[ j ].mv( ds, Hds );
				ret[ j ] -= ds * Hds;
			  }
			}
			RangeType uVal;
			u( y, uVal );
			ret += uVal;
		  }     // TODO ENDE IM FOLGENDEN ANGENOMMEN Hu IST HESSEMATRIX
		
		
		  RangeType IN( 39.047*Hu[0,0]-20.448*Hu[0,1]-20.448*H[1,0]+11.345*Hu[1,1]);
		    
		  // BRAUCHEN WIR NICHT, DA DIREKT IN CODE GESCHRIEBENmodel_.diffusiveFlux(entity, x , u_x , Du_x , newDu_x );
		  
		  // multiply by quadrature weight
		  //newu_x *= quadrature.weight( pt ) * geometry.integrationElement( x );
		  //newDu_x *= quadrature.weight( pt ) * geometry.integrationElement( x );
		  newIN *= quadrature.weight( pt ) * geometry.integrationElement( x );
		  
		  // add f * phi_i to rhsLocal[ i ]
		  w_local.axpy( quadrature[ pt ], newIN );
		  
	      break ;  
	      }
			
		
		
		
		
		
      
    }
  }
  
  
  
  
  
  
  
/********************************************************************************/
  // communicate data (in parallel runs)
  w.communicate();

}

#endif // #ifndef ELLIPTIC_HH
