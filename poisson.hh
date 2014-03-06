#ifndef POISSON_PROBLEMS_HH
#define POISSON_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "probleminterface.hh"

// -laplace u + u = f with Neumann boundary conditions on domain [0,1]^d
// Exact solution is u(x_1,...,x_d) = cos(2*pi*x_1)...cos(2*pi*x_d)
template <class FunctionSpace> 
class Enthalten : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:  
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& rechteSeite) const
  {
    rechteSeite = 0;
  }

  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(1);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x, 
                 RangeType& value) const 
  {
    value = RangeType( 0 );
  }
};
// -laplace u = f with zero Dirichlet boundary conditions on domain [0,1]^d
// Exsct solution is u(x_1,...,x_d) = sin(2*pi*x_1)...sin(2*pi*x_d)
template <class FunctionSpace> 
class SinusProduct : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:  
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi = 4*dimDomain*(M_PI*M_PI);
    for( int i = 0; i < dimDomain; ++i )
      phi *= std::sin( 2*M_PI*x[ i ] );
  }

 
  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const 
  {
    return true;
  }
  virtual bool hasDirichletBoundary () const 
  {
    return true ;
  }
};

// A problem on a unit sphere
template <class FunctionSpace> 
class SphereProblem : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:  
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& ret) const
  {
    DomainType y = x;
    y /= x.two_norm();

    JacobianRangeType Ju;
    uJacobian( y, Ju );
    Ju.mv( y, ret );
    ret *= dimDomain - 1;

    RangeType uVal;
    u( y, uVal );
    ret += uVal;
  }

 
  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(1);
  }


};

#endif // #ifndef POISSON_HH
