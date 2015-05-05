/** 
 * \file geodesic_point_vertex.hpp
 *  
 * \brief Definition of  class Geodesic_point_vertex, which represents
 * a vertex of  a discrete geodesic (an open,  simple polygonal line).
 * The  vertex coincides  with a  vertex of  the mesh  over  which the
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

#ifndef ADG_GEODESIC_POINT_VERTEX_HPP
#define ADG_GEODESIC_POINT_VERTEX_HPP

#include "mesh_interface.hpp"         // common::Mesh_interface
#include "geodesic_point.hpp"         // Geodesic_point_vertex

#include <cassert>                    // assert


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
   * \class Geodesic_point_vertex
   *
   * \brief This class represents a  vertex of a discrete geodesic (an
   * open, simple polygonal line) that  coincides with a vertex of the
   * mesh over which the geodesic is defined.
   */
  template< typename Mesh >
  class Geodesic_point_vertex : public Geodesic_point< Mesh > {
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
     * \fn Geodesic_point_vertex( Vertex* v )
     *
     * \brief Creates an  instance of this class.
     *
     * \param v Pointer to a vertex  of the mesh on which the geodesic
     * is defined.
     */
    Geodesic_point_vertex( Vertex* v )
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
#endif

      _vertex = v ;

      return ;
    }


    /**
     * \fn Geodesic_point_vertex( const Geodesic_point_vertex& gv )
     *
     * \brief Creates an instance  of this class form another instance
     * of this class.
     *
     * \param gv A reference to an instance of this class.
     */
    Geodesic_point_vertex( const Geodesic_point_vertex& gv )
    {
      _vertex = gv.get_vertex() ;
    }


    /**
     * \fn ~Geodesic_point_vertex()
     *
     * \brief Destroys an instance of this class.
     */    
    ~Geodesic_point_vertex()
    {}


    /**
     * \fn inline Vertex* get_vertex() const
     *
     * \brief Returns a pointer to the mesh vertex that corresponds to
     * this geodesic vertex. 
     *
     * \return A pointer  to the mesh vertex that  corresponds to this
     * geodesic vertex.
     */
    inline Vertex* get_vertex() const 
    {
      return _vertex ;
    }


    /**
     * \fn inline void set_vertex( Vertex* v )
     *
     * \brief Assigns a mesh vertex to this geodesic vertex.
     *
     * \param  v  Pointer to  a  vertex of  the  mesh  over which  the
     * geodesic is defined.
     */
    inline void set_vertex( Vertex* v ) 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
#endif

      _vertex = v ;

      return ;
    }


    /**
     * \fn inline bool on_vertex() const
     *
     * \brief Returns the logic value  true to indicate that this is a
     * geodesic point that coincides with a mesh vertex.
     *
     * \return  The logic  value true.
     */    
    inline bool on_vertex() const 
    {
      return true ;
    }


    /**
     * \fn inline bool on_edge() const
     *
     * \brief Returns the logic value false to indicate that this is a
     * geodesic point that  does not belong to the  interior of a mesh
     * edge.
     *
     * \return The logic value false.
     */  
    inline bool on_edge() const
    {
      return false ;
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
     * \brief  Gets the  three  coordinates of  the  mesh vertex  that
     * coincides with this geodesic vertex.
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

      m->get_coords( _vertex , x , y , z ) ;
      
      return ;
    }


    /**
     * \fn bool in_the_same_mesh_element( MI* m , GP* v , std::list< GP* >& lv ) const
     *
     * \brief  Tests whether  this vertex  and another  given geodesic
     * vertex correspond to  the same mesh vertex.  If  so, place this
     * geodesic  vertex inside a  given vertex  list, and  returns the
     * logic value  true. Otherwise, it  does nothing and  returns the
     * logic value false.
     *
     * \param  m  Pointer to  the  mesh  that  implements the  generic
     * interface.
     * \param v Pointer to a geodesic vertex.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return  The logic  value  true if  the  given geodesic  vertex
     * coincides with the same mesh  vertex as this one, and the logic
     * value false otherwise.
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

      if ( v->on_vertex() ) {
	Geodesic_point_vertex< Mesh >* gpv = 
	  dynamic_cast< Geodesic_point_vertex< Mesh >* >( v ) ;

	if ( _vertex == gpv->get_vertex() ) {
	  lv.push_back( v ) ;
	  return true ;
	}
      }
      
      return false ;
    }


    /**
     * \fn void get_adj_vertices( MI* m , std::vector< Vertex* >& a ) const
     *
     * \brief  Finds  all mesh  vertices  adjacent  to  the mesh  vertex
     * correponding to this geodesic vertex, and place them in the given
     * list.
     *
     * \param m Pointer to the mesh over which the geodesic is defined.
     * \param a A reference to a list of adjacent mesh vertices.
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

      Halfedge* h1 = m->get_halfedge( _vertex ) ;
      Halfedge* h2 = h1 ;
      
      do {
	a.push_back( m->get_org( m->get_next( h2 ) ) ) ;
	h2 = m->get_next( m->get_mate( h2 ) ) ;
      } 
      while ( h2 != h1 ) ;
    
      return ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Geodesic_point_vertex()
     *
     * \brief Creates an instance of this class.
     */
    Geodesic_point_vertex() 
    {
      _vertex = 0 ;
    }


    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    Vertex* _vertex ;  ///< Pointer to the mesh vertex corresponding to this vertex (if any).

  } ;

}

/** @} */ //end of group class.

#endif	/* GEODESIC_POINT_VERTEX_HPP */

