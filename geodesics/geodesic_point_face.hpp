/** 
 * \file geodesic_point_face.hpp
 *  
 * \brief Definition of  class Geodesic_point_face, which represents a
 * vertex  of a discrete  geodesic (an  open, simple  polygonal line).
 * The vertex belongs to the interior of a face of the mesh over which
 * the geodesic is defined.
 *
 * \author
 * Mario Augusto de Souza Lizier
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

#ifndef GEODESIC_POINT_FACE_HPP
#define GEODESIC_POINT_FACE_HPP

#include "mesh_interface.hpp"         // common::Mesh_interface
#include "geodesic_point.hpp"         // Geodesic_point

#include <cassert>                    // assert
#include <cmath>                      // fabs


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
   * \class Geodesic_point_face
   *
   * \brief This class represents a  vertex of a discrete geodesic (an
   * open, simple polygonal line). The  vertex lies in the interior of
   * a  face of  the mesh  over which  the geodesic  it belongs  to is
   * defined.
   */
  template< typename Mesh >
  class Geodesic_point_face : public Geodesic_point< Mesh > {
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
     * \fn Geodesic_point_face( Face* f , double u , double v )
     *
     * \brief Creates an  instance of this class.
     *
     * \param f Pointer to a face of the mesh on which the geodesic is
     * defined.
     * \param u  First barycentric coordinate of a  point with respect
     * to the affine  frame defined by the vertices  of the given mesh
     * face.
     * \param v Second barycentric  coordinate of a point with respect
     * to the affine  frame defined by the vertices  of the given mesh
     * face.
     */
    Geodesic_point_face( 
			Face* f ,
			double u ,
			double v
		       ) 
    {
#ifdef DEBUGMODE
      assert( 
	     ( f != 0 ) && 
	     ( fabs( u ) > 0 ) &&  ( fabs( v ) > 0 ) &&
	     ( fabs( u + v ) < 1.0 )
	    ) ;
#endif

      _face = f ;
      _u = u ;
      _v = v ;

      return ;
    }


    /**
     * \fn Geodesic_point_face( const Geodesic_point_face& gf )
     *
     * \brief Creates an instance  of this class form another instance
     * of this class.
     *
     * \param gf A reference to an instance of this class.
     */
    Geodesic_point_face( const Geodesic_point_face& gf )
    {
      _face = gf.get_face() ;
      _u = gf.get_u() ;
      _v = gf.get_v() ;

      return ;
    }


    /**
     * \fn ~Geodesic_point_face()
     *
     * \brief Destroys an instance of this class.
     */    
    ~Geodesic_point_face()
    {}


    /**
     * \fn inline Face* get_face() const
     *
     * \brief Returns  a pointer to  the mesh face that  contains this
     * geodesic vertex.
     *
     * \return A pointer to the  mesh face that contains this geodesic
     * vertex.
     */
    inline Face* get_face() const 
    {
      return _face ;
    }


    /**
     * \fn inline void set_face( Face* f , double u , double v )
     *
     * \brief Makes this geodesic vertex a vertex in the interior of a
     * mesh face.
     *
     * \param f Pointer to a mesh face.
     * \param u  First barycentric coordinate of a  point with respect
     * to the affine  frame defined by the vertices  of the given mesh
     * face.
     * \param v Second barycentric  coordinate of a point with respect
     * to the affine  frame defined by the vertices  of the given mesh
     * face.
     */
    inline void set_face( 
			 Face* f ,
			 double u ,
			 double v
			) 
    {
#ifdef DEBUGMODE
      assert( 
	     ( f != 0 ) && 
	     ( fabs( u ) > 0 ) &&  ( fabs( v ) > 0 ) &&
	     ( fabs( u + v ) < 1.0 )
	    ) ;
#endif

      _face = f ;
      _u = u ;
      _v = v ;

      return ;
    }


    /**
     * \fn inline double get_u() const
     *
     * \brief  Returns  the   first  barycentric  coordinate  of  this
     * geodesic vertex with respect to the affine frame defined by the
     * vertices of the mesh face that contains this vertex.
     *
     * \return  The  first  barycentric  coordinate of  this  geodesic
     * vertex with respect to the affine frame defined by the vertices
     * of the mesh face that contains this vertex.
     */
    inline double get_u() const 
    {
      return _u ;
    }


    /**
     * \fn inline double get_v() const
     *
     * \brief  Returns  the  second  barycentric  coordinate  of  this
     * geodesic vertex with respect to the affine frame defined by the
     * vertices of the mesh face that contains this vertex.
     *
     * \return  The  second barycentric  coordinate  of this  geodesic
     * vertex with respect to the affine frame defined by the vertices
     * of the mesh face that contains this vertex.
     */
    inline double get_v() const 
    {
      return _v ;
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
     * \brief Returns the logic value  true to indicate that this is a
     * geodesic point that belongs to the interior of a mesh face.
     *
     * \return The logic value true.
     */  
    inline bool on_face() const
    {
      return true ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Geodesic_point_face()
     *
     * \brief Creates an  instance of this class.
     */
    Geodesic_point_face() 
    {
      _face = 0 ;

      _u = -1 ;
      _v = -1 ;

      return ;
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
      
      Halfedge* h = m->get_halfedge( _face ) ;

      double x1 , y1 , z1 ;
      m->get_coords( m->get_org( h ) , x1 , y1 , z1 ) ;

      double x2 , y2 , z2 ;
      m->get_coords( m->get_org( m->get_next( h ) ) , x2 , y2 , z2 ) ;

      double x3 , y3 , z3 ;
      m->get_coords( m->get_org( m->get_prev( h ) ) , x3 , y3 , z3 ) ;

      double w = 1 - _u - _v ;
    
      /*
       * The following is to reduce the impact of cancelation.
       */
      if ( fabs( w ) < 1e-14 ) {
	w = 0 ;
      }
      else if ( fabs( 1 - w ) < 1e-14 ) {
	w = 1 ;
      }
    
      x = _u * x1 + _v * x2 + w * x3 ;
      y = _u * y1 + _v * y2 + w * y3 ;
      z = _u * z1 + _v * z2 + w * z3 ;

      return ;
    }


    /**
     * \fn bool in_the_same_mesh_element( MI* m , GP* v , std::list< GP* >& lv ) const
     *
     * \brief  Tests whether  this vertex  and another  given geodesic
     * vertex correspond to  the same mesh vertex.  If  so, place this
     * geodesic  vertex inside  a given  vertex list,  and  return the
     * logic value  true. Otherwise, do  nothing and return  the logic
     * value false.
     *
     * \param  m  Pointer to  the  mesh  over  which the  geodesic  is
     * defined.
     * \param v Pointer to a geodesic vertex.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return  The logic  value  true if  the  given geodesic  vertex
     * coincides with  the same mesh  vertex as this  geodesic vertex,
     * and the logic value false otherwise.
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

      if ( v->on_face() ) {
	Geodesic_point_face< Mesh >* gpv = 
	  dynamic_cast< Geodesic_point_face< Mesh >* >( v ) ;

	if ( _face == gpv->get_face() ) {
	  if ( ( _u == gpv->get_u() ) && ( _v == gpv->get_v() ) ) {
	    lv.push_back( const_cast< Geodesic_point_face< Mesh >* >( this ) ) ;
	  }
	  else {
	    lv.push_back( const_cast< Geodesic_point_face< Mesh >* >( this ) ) ;
	    lv.push_back( v ) ;
	  }
	  return true ;
	}
      }
      
      return false ;
    }


    /**
     * \fn void get_adj_vertices( MI* m , std::vector< Vertex* >& a ) const
     *
     * \brief Finds all mesh  vertices incident to the face containing
     * this geodesic vertex, and place them in the given list.
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
 
      Halfedge* h1 = m->get_halfedge( _face ) ;

      a.push_back( m->get_org( h1 ) ) ;
      a.push_back( m->get_org( m->get_next( h1 ) ) ) ;
      a.push_back( m->get_org( m->get_prev( h1 ) ) ) ;
    
      return ;
    }

  private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    Face* _face ;  ///< Pointer to the mesh face containing this vertex in its interior (if any).
    
    double _u ; ///< First barycentric coordinate of this vertex with respect to the affine frame defined by the vertices of a mesh face.

    double _v ; ///< Second barycentric coordinate of this vertex with respect to the affine frame defined by the vertices of a mesh face.

  } ;


}

/** @} */ //end of group class.

#endif	/* GEODESIC_POINT_FACE_HPP */
