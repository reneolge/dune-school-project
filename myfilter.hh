#ifndef DUNE_FEM_SCHOOL_INTRO_FILTER_HH
#define DUNE_FEM_SCHOOL_INTRO_FILTER_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-fem includes
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>

namespace School
{

  // Filter
  // ------------

  /*
   * Given a center x and a radius r, a point y is in the filter, 
   * iff it lies within the ball of radius r around x, i. e.,
   * if y in B_r( x ).
   *
   * \note This class implements the filter interface is documented 
   *        in <dune/fem/gridpart/filter/filter.hh>
   */
  template< class HostGridPart >
  class Filter
  {
    typedef HostGridPart HostGridPartType;

    class BasicFilter;

    // this class is implemented through a helper class,
    // see file <dune/fem/gridpart/filter/basicfilterwrapper.hh>
    typedef Dune::Fem::BasicFilterWrapper< HostGridPartType, BasicFilter > ImplementationType;

  public:
    // single coordinate type
    typedef typename BasicFilter::ctype ctype;
    // global coordinate type
    typedef typename BasicFilter::GlobalCoordinateType GlobalCoordinateType;
 
    template< int cd >
    struct Codim
    {
      // entity of arbitrary codimension in host grid part
      typedef typename HostGridPartType::template Codim< cd >::EntityType EntityType;
    };

    // element (codim 0 entity) type
    typedef typename Codim< 0 >::EntityType EntityType;
   
    // constructor
    Filter ( const HostGridPartType &hostGridPart)
    : implementation_( hostGridPart ) 
    {}

    // returns true, if host grid part entity shall be contained in filtered grid part
    template< int cd >
    bool contains ( const typename Codim< cd >::EntityType &entity ) const
    {
      return implementation_.template contains< cd >( entity );
    }

    // returns true, if host grid part entity shall be contained in filtered grid part
    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      return contains< Entity::codimension >( entity );
    }

    // for more information, see <dune/fem/gridpart/filter/filter.hh>
    template< class Intersection >
    bool interiorIntersection ( const Intersection &intersection ) const
    {
      typedef typename Intersection::EntityPointer EntityPointerType;
      const EntityPointerType outside = intersection.outside();
      return contains( *outside );
    }

    // for more information, see <dune/fem/gridpart/filter/filter.hh>
    template< class Intersection >
    bool intersectionBoundary( const Intersection &intersection ) const
    {
      return true;
    }

    // for more information, see <dune/fem/gridpart/filter/filter.hh>
    template< class Intersection >
    bool intersectionNeighbor ( const Intersection &intersection ) const
    {
      return false;
    }

    // for more information, see <dune/fem/gridpart/filter/filter.hh>
    template< class Intersection >
    int intersectionBoundaryId ( const Intersection &intersection ) const
    {
      return 1;
    }

  private:
    ImplementationType implementation_;
  };

  // Implementation of FilterType< HostGridPart >::BasicFilter
  // ---------------------------------------------------------------

  template< class HostGridPart >
  class Filter< HostGridPart >::BasicFilter
  {
  public:
    typedef typename HostGridPart::ctype ctype;
    typedef Dune::FieldVector< ctype, HostGridPart::dimensionworld > GlobalCoordinateType;

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      const int codim = Entity::codimension;
      if( codim != 0 )
        DUNE_THROW( Dune::InvalidStateException, 
                    "BasicFilter::contains() only available for codim 0 entities." );
      return contains( entity.geometry().center() );
    }

  private:
  const double drehwinkel = 0.2976893343;
  const double streckung = 2.;
  double drehmatrix[2][2]=
	{
		{streckung*std::cos(drehwinkel) , -streckung*std::sin(drehwinkel)},
		{streckung*std::sin(drehwinkel) , streckung*std::cos(drehwinkel)}
	};
    bool contains ( const GlobalCoordinateType &x ) const
    {
      ctype punkt_gedreht (drehmatrix*x);
      
      return (-1 < punkt_gedreht[0] < 1 && -1 < punkt_gedreht[1] < 1);
    }
  };

} // namespace School

#endif // #ifndef DUNE_FEM_SCHOOL_INTRO_FILTER_HH
