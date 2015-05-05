/** 
 * \file geodesic_point_edge.hpp
 *  
 * \brief Definition of  class Geodesic_point_edge, which represents a
 * vertex  of a discrete  geodesic (an  open, simple  polygonal line).
 * The vertex is in the interior of an edge of the mesh over which the
 * geodesic is defined.
 *
 * \author
 * Mario Augusto de Souza Lizier \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * lizier at gmail (dot) com
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * mfsiqueira at gmail (dot) com
 *
 * \version 1.0
 * \date January 2010
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef GEODESIC_POINT_EDGE_HPP
#define GEODESIC_POINT_EDGE_HPP

#include "mesh_interface.hpp"        // common::Mesh_interface
#include "geodesic_point.hpp"        // Geodesic_point

#include <cassert>                   // assert
#include <cmath>                     // fabs


/** 
 * \defgroup DGLNameSpace Namespace dg.
 * @{
 */

/**
 * \namespace dg
 *
 * \brief The namespace dgl contains the definition and implementation
 * of       classes       Adg_mesh_interface,      Edg_mesh_interface,
 * Approx_geodesics,          Exact_geodesics,         Geodesic_point,
 * Geodesic_point_vertex,   Geodesic_point_edge,  Geodesic_point_face,
 * and Curve_point,  which allow us  to compute approximate  and exact
 * discrete  geodesics   over  triangle  surface   meshes  with  empty
 * boundary.
 */


namespace dg {

  /**
   * \class Geodesic_point_edge
   *
   * \brief This class represents a  vertex of a discrete geodesic (an
   * open, simple polygonal line). The  vertex lies in the interior of
   * an edge  of the  mesh over  which the geodesic  it belongs  to is
   * defined.
   */
  template< typename Mesh >
  class Geodesic_point_edge : public Geodesic_point< Mesh > {
  public:
    // ---------------------------------------------------------------
    //
    // Type name definitions.
    //
    // ---------------------------------------------------------------

    /**
     * \typedef MI
     *
     * \brief Definition of a type name for the mesh interface.
     */
    typedef common::Mesh_interface< Mesh > MI ;

    /**
     * \typedef Vertex
     *
     * \brief Definition of a type name for the mesh vertices.
     */
    typedef typename MI::Vertex Vertex ;


    /**
     * \typedef Halfedge
     *
     * \brief Definition of a type name for the mesh half-edges.
     */
    typedef typename MI::Halfedge Halfedge ;


    /**
     * \typedef Edge
     *
     * \brief Definition of a type name for the mesh edges.
     */
    typedef typename MI::Edge Edge ;


    /**
     * \typedef Face
     *
     * \brief Definition of a type name for the mesh faces.
     */
    typedef typename MI::Face Face ;


    /**
     * \typedef GP
     *
     * \brief Definition of a type name for the parent class.
     */
    typedef Geodesic_point< Mesh > GP ;


    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Geodesic_point_edge( Edge* e , double t )
     *
     * \brief Creates an  instance of this class.
     *
     * \param e Pointer  to an edge of the mesh  on which the geodesic
     * is defined.
     * \param t  The ratio that  indicates how far away  this geodesic
     * vertex is from the origin vertex of the first half-edge of edge
     * \a e.
     */
    Geodesic_point_edge( 
			Edge* e ,
			double t 
		       ) 
    {
#ifdef DEBUGMODE
      assert( 
	     ( e != 0 ) && 
	     ( fabs( t ) > 0 ) &&  ( fabs( t ) < 1 ) 
	    );
#endif

      _edge = e ;
      _t = t ;

      return ;
    }


    /**
     * \fn Geodesic_point_edge( const Geodesic_point_edge& ge )
     *
     * \brief Creates an instance  of this class form another instance
     * of this class.
     *
     * \param ge A reference to an instance of this class.
     */
    Geodesic_point_edge( const Geodesic_point_edge& ge )
    {
      _edge = ge.get_edge() ;
      _t = ge.get_t() ;

      return ;
    }


    /**
     * \fn ~Geodesic_point_edge()
     *
     * \brief Destroys an instance of this class.
     */    
    ~Geodesic_point_edge()
    {}


    /**
     * \fn inline Edge* get_edge() const
     *
     * \brief Returns  a pointer to  the mesh edge that  contains this
     * geodesic vertex.
     *
     * \return A pointer to the  mesh edge that contains this geodesic
     * vertex.
     */
    inline Edge* get_edge() const 
    {
      return _edge ;
    }


    /**
     * \fn inline void set_edge( Edge* e , double t )
     *
     * \brief Makes this geodesic vertex a vertex in the interior of a
     * mesh edge.
     *
     * \param e Pointer to an edge of the mesh over which the geodesic
     * is defined.
     *
     * \param t  The ratio  that indicates how  far away  the geodesic
     * vertex is from the origin mesh vertex of the first half-edge of
     * the edge \a e.
     */
    inline void set_edge( 
			 Edge* e ,
			 double t 
			) 
    {
#ifdef DEBUGMODE
      assert( 
	     ( e != 0 ) && 
	     ( fabs( t ) > 0 ) &&  ( fabs( t ) < 1 ) 
	    );
#endif

      _edge = e ;
      _t = t ;

      return ;
    }


    /**
     * \fn inline double get_t() const
     *
     * \brief  Returns  the ratio  that  indicates  how  far away  the
     * geodesic vertex  is from  the origin mesh  vertex of  the first
     * half-edge of the edge this vertex lies on.
     *
     * \return  The ratio  that indicates  how far  away  the geodesic
     * vertex is from the origin mesh vertex of the first half-edge of
     * the edge this vertex lies on.
     */
    inline double get_t() const 
    {
      return _t ;
    }


    /**
     * \fn inline bool on_vertex() const
     *
     * \brief Returns the  logic value false to indicate  that this is
     * not a geodesic point that coincides with a mesh vertex.
     *
     * \return The logic value false.
     */    
    inline bool on_vertex() const 
    {
      return false ;
    }


    /**
     * \fn inline bool on_edge() const
     *
     * \brief Returns the logic value  true to indicate that this is a
     * geodesic point that belongs to the interior of a mesh edge.
     *
     * \return The logic value true.
     */  
    inline bool on_edge() const
    {
      return true ;
    }


    /**
     * \fn inline bool on_face() const
     *
     * \brief  Returns the logic value false to indicate that this is a
     * geodesic point that  does not belong to the  interior of a mesh
     * face.
     *
     * \return The logic value false.
     */  
    inline bool on_face() const
    {
      return false ;
    }


    /**
     * \fn void get_coords( MI* m, double& x , double& y , double& z ) const
     *
     * \brief  Gets the  three  coordinates of  this geodesic  vertex,
     * which lies in the interior of an edge of the mesh.
     *
     * \param  m  Pointer to  the  mesh  over  which the  geodesic  is
     * defined.
     * \param x A reference to the X coordinate.
     * \param y A reference to the Y coordinate.
     * \param z A reference to the Z coordinate.
     */
    void get_coords(
		    MI* m , 
		    double& x , 
		    double& y , 
		    double& z 
		   ) 
      const
    {
#ifdef DEBUGMODE
      assert( m != 0 ) ;
#endif

      Halfedge* h1 ;
      Halfedge* h2 ;
      m->get_halfedges( _edge , h1 , h2 ) ;

      double x1 , y1 , z1 ;
      m->get_coords( m->get_org( h1 ) , x1 , y1 , z1 ) ;

      double x2 , y2 , z2 ;
      m->get_coords( m->get_org( h2 ) , x2 , y2 , z2 ) ;
      
      x = x1 + _t * ( x2 - x1 ) ;
      y = y1 + _t * ( y2 - y1 ) ;
      z = z1 + _t * ( z2 - z1 ) ;
      
      return ;
    }


    /**
     * \fn bool in_the_same_mesh_element( MI* m , GP* v , std::list< GP* >& lv ) const
     *
     * \brief  Tests whether  this vertex  and another  given geodesic
     * vertex lie  in the interior  of the same  mesh edge, or  in the
     * interior of distinct edges tha share the same mesh face. If so,
     * place both vertices  in the given list of  geodesic vertices if
     * thei  positions are  distinct; otherwise,  inserts only  one of
     * them in the  list.  In either case, returns  true.  However, if
     * the  vertices do  not lie  in the  interior of  the  same edge,
     * returns the logic value false.
     *
     * \param  m  Pointer to  the  mesh  that  implements the  generic
     * interface.
     * \param v Pointer to a geodesic vertex.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return The  logic value true  if this geodesic vertex  and the
     * given  geodesic vertex  lie in  the interior  of the  same mesh
     * edge, and the logic value false otherwise.
     */
    bool in_the_same_mesh_element(
				  MI* m ,
				  GP* v ,
				  std::list< GP* >& lv 
				 )
      const
    {
#ifdef DEBUGMODE
      assert( m != 0 ) ;
      assert( v != 0 ) ;
      assert( lv.empty() ) ;
#endif

      if ( v->on_edge() ) {
	Geodesic_point_edge< Mesh >* gpv = 
	  dynamic_cast< Geodesic_point_edge< Mesh >* >( v ) ;

	Edge* e = gpv->get_edge() ;

	if ( _edge == e ) {
	  double t = gpv->get_t() ;
	  if ( _t == t ) {
	    lv.push_back( v ) ;
	  }
	  else {
	    lv.push_back( const_cast< Geodesic_point_edge< Mesh >* >( this ) ) ;
	    lv.push_back( v ) ;
	  }

	  return true ;
	}
	else {
	  //
	  // Are the edges in the same face?
	  //
	  Halfedge* hi1 ;
	  Halfedge* hi2 ;
	  m->get_halfedges( _edge , hi1 , hi2 ) ;

	  Halfedge* he1 ;
	  Halfedge* he2 ;
	  m->get_halfedges( e , he1 , he2 ) ;	
	  if ( 
	      ( m->get_next( hi1 ) == he1 ) || 
	      ( m->get_prev( hi1 ) == he1 ) ||
	      ( m->get_next( hi1 ) == he2 ) || 
	      ( m->get_prev( hi1 ) == he2 ) ||
	      ( m->get_next( hi2 ) == he1 ) || 
	      ( m->get_prev( hi2 ) == he1 ) ||
	      ( m->get_next( hi2 ) == he2 ) || 
	      ( m->get_prev( hi2 ) == he2 )
	     ) {
	    lv.push_back( const_cast< Geodesic_point_edge< Mesh >* >( this ) ) ;
	    lv.push_back( v ) ;
	    
	    return true ;
	  }
	}
      }
      
      return false ;
    }


    /**
     * \fn void get_adj_vertices( MI* m , std::vector< Vertex* >& a ) const
     *
     * \brief Finds  all mesh  vertices in the  star of the  mesh edge
     * that contains  this geodesic vertex in its  interior (there are
     * four vertices).
     *
     * \param m Pointer to the mesh over which the geodesic is defined.
     * \param a A reference to a list of mesh vertices.
     */
    void get_adj_vertices( 
			  MI* m , 
			  std::vector< Vertex* >& a 
			 ) 
      const
    {
#ifdef DEBUGMODE
      assert( m != 0 ) ;
#endif

      Halfedge* h1 ;
      Halfedge* h2 ;
      m->get_halfedges( _edge , h1 , h2 ) ;

      a.push_back( m->get_org( h1 ) ) ;
      a.push_back( m->get_org( m->get_prev( h2 ) ) ) ;
      a.push_back( m->get_org( h2 ) ) ;
      a.push_back( m->get_org( m->get_prev( h1 ) ) ) ;
    
      return ;
    }

  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Geodesic_point_edge()
     *
     * \brief Creates an instance of this class.
     */
    Geodesic_point_edge() 
    {
      _edge = 0 ;
      _t = -1 ;

      return ;
    }


    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    Edge* _edge ;  ///< Pointer to the mesh edge containing this vertex in its interior (if any).

    double _t ;  ///< Ratio between the distance of this vertex to one end point of a mesh edge and the length of the edge.

  } ;

}

/** @} */ //end of group class.

#endif	/* GEODESIC_POINT_EDGE_HPP */
