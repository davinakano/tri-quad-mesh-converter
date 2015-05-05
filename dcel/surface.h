/** 
 * \file surface.h
 *  
 * \brief Definition  and implementation  of the class  Surface, which
 * represents  a triangle  surface mesh  (i.e., a  simplicial surface)
 * with empty boundary.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
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

#ifndef SURFACE_H
#define SURFACE_H

#include "vertex.h"    // Vertex
#include "halfedge.h"  // Halfedge
#include "edge.h"      // Edge
#include "face.h"      // Face

#include <cassert>     // assert
#include <list>        // std::list
#include <map>         // std::map


/** 
* \defgroup DCELNameSpace Namespace dcel.
* @{
*/

/**
* \namespace dcel
*
* \brief The namespace dcel contains the definition and implementation
* of  all classes  of the  data structure  Doubly Connected  Edge List
* (DCEL), which  can be  used to represent  surface meshes  with empty
* boundary.
*
*/

namespace dcel {


  /**
   * \class Surface
   *
   * \brief This  class represents a triangle surface  mesh with empty
   * boundary.
   */
  template < 
             typename VAttrib = int ,
             typename FAttrib = int ,
             typename EAttrib = int ,
             typename HAttrib = int 
           >
  class Surface {
  public:
    // ---------------------------------------------------------------
    //
    // Type definitions
    //
    // ---------------------------------------------------------------

    /**
     * \typedef Vertex
     *
     * \brief Defines  Vertex as  an alias for  dcel::Vertex< VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Vertex< VAttrib, FAttrib , EAttrib , HAttrib > Vertex ;

    /**
     * \typedef Face
     *
     * \brief  Defines  Face  as  an alias  for  dcel::Face<  VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Face< VAttrib, FAttrib , EAttrib , HAttrib > Face ;

    /**
     * \typedef Edge
     *
     * \brief  Defines  Edge  as  an alias  for  dcel::Edge<  VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Edge< VAttrib, FAttrib , EAttrib , HAttrib > Edge ;

    /**
     * \typedef Halfedge
     *
     * \brief  Defines  Halfedge   as  an  alias  for  dcel::Halfedge<
     * VAttrib, FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Halfedge< VAttrib, FAttrib , EAttrib , HAttrib > Halfedge ;

    /**
     * \typedef VertexIterator
     *
     * \brief Defines  an iterator  for the list  of vertices  of this
     * mesh.
     */    
    typedef typename std::list< Vertex* >::const_iterator VertexIterator ;


    /** 
     * \typedef EdgeIterator
     *
     * \brief Defines an iterator for the list of edge of this mesh.
     */    
    typedef typename std::list< Edge* >::const_iterator EdgeIterator ;


    /** 
     * \typedef FaceIterator
     *
     * \brief Defines an iterator for the list of faces of this mesh.
     */    
    typedef typename std::list< Face* >::const_iterator FaceIterator ;


    /** 
     * \typedef VMAP
     *
     * \brief  Defines a hash  table type  for vertex  identifiers and
     * pointers.
     */
    typedef typename std::map< unsigned , Vertex* > VMAP ;


    /** 
     * \typedef HMAP
     *
     * \brief  Defines  a  hash  table of  half-edge  identifiers  and
     * pointers.
     */
    typedef typename std::map< std::pair< unsigned , unsigned > ,
                                               Halfedge* > HMAP ;
 

    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Surface( unsigned nverts , double* lverts , unsigned nfaces , unsigned* lfaces )
     *
     * \brief Creates an instance of this class.
     *
     * \param nverts The number of vertices of this surface.
     * \param lverts An array with the vertex coordinates.
     * \param nfaces The number of faces of this surface.
     * \param lfaces An array with the vertex identifiers of each face
     * of this surface.
     */
    Surface(
	    unsigned nverts ,
	    double* lverts ,
	    unsigned nfaces ,
	    unsigned* lfaces
	   ) ;


    /**
     * \fn ~Surface()
     *
     * \brief Releases the memory held fy an instance of this class.
     */
    virtual ~Surface() ;


    /** 
     * \fn inline unsigned get_number_of_vertices() const
     *
     * \brief Returns the number of vertices of this surface.
     *
     * \return The number of vertices of this surface.
     */
    inline unsigned get_number_of_vertices() const
    {
      return _list_of_vertices.size() ;
    }


    /** 
     * \fn inline unsigned get_number_of_edges() const
     *
     * \brief Returns the number of edges of this surface.
     *
     * \return The number of edges of this surface.
     *
     */
    inline unsigned get_number_of_edges() const
    { 
      return _list_of_edges.size() ; 
    }


    /**
     * \fn inline unsigned get_number_of_faces() const 
     *
     * \brief Returns the number of faces of this surface
     *
     * \return The number of faces of this surface
     */
    inline unsigned get_number_of_faces() const
    { 
      return _list_of_faces.size() ; 
    }


    /**
     * \fn inline VertexIterator vertices_begin() const
     *
     * \brief Returns a vertex iterator set to the first vertex of the
     * list of vertices of this surface.
     *
     * \return A vertex  list iterator set to the  first vertex of the
     * list.
     */
    inline VertexIterator vertices_begin() const
    { 
      return _list_of_vertices.begin() ;
    }


    /**
     * \fn inline EdgeIterator edges_begin() const 
     *
     * \brief Returns  an edge iterator set  to the first  edge of the
     * list of edges of this surface.
     *
     * \return  An edge list  iterator set  to the  first edge  of the
     * list.
     */
    inline EdgeIterator edges_begin() const
    {
      return _list_of_edges.begin() ;  
    }


    /**
     * \fn FaceIterator faces_begin() const
     *
     * \brief Returns  a face  iterator set to  the first face  of the
     * list of faces of this surface.
     *
     * @return A face list iterator set to the first face of the list.
     */
    FaceIterator faces_begin() const
    {
      return _list_of_faces.begin() ;  
    }


    /**
     * \fn inline VertexIterator vertices_end() const
     *
     * \brief Returns  a vertex iterator  set past the last  vertex of
     * the list of vertices of this surface.
     *
     * \return A vertex list iterator  set past the last vertex of the
     * list.
     */
    inline VertexIterator vertices_end() const
    { 
      return _list_of_vertices.end() ;
    }


    /**
     * \fn inline EdgeIterator edges_end() const 
     *
     * \brief Returns  an edge iterator set past the last  edge of the
     * list of edges of this surface.
     *
     * \return An  edge list  iterator set past  the last edge  of the
     * list.
     */
    inline EdgeIterator edges_end() const
    {
      return _list_of_edges.end() ;  
    }


    /**
     * \fn FaceIterator faces_end() const
     *
     * \brief Returns  a face iterator set  past the last  face of the
     * list of faces of this surface.
     *
     * @return A face list iterator set past the last face of the list.
     */
    FaceIterator faces_end() const
    {
      return _list_of_faces.end() ;  
    }

    /**
     * \fn void add_vertex( Vertex* v )
     *
     * \brief Performs an insertion of a vertex in the mesh.
     *
     * \param v A pointer of the vertex to be inserted.
     */
    void add_vertex( Vertex* v ) {
        _list_of_vertices.push_back( v ) ;
    }


    /**
     * \fn void add_edge( Edge* e )
     *
     * \brief Performs an insertion of an edge in the mesh.
     *
     * \param e A pointer of the edge to be inserted.
     */
    void add_edge( Edge* e ) {
        _list_of_edges.push_back( e ) ;
    }


    /**
     * \fn void add_face( Face* f )
     *
     * \brief Performs an insertion of a face in the mesh.
     *
     * \param f A pointer of the face to be inserted.
     */
    void add_face( Face* f ) {
        _list_of_faces.push_back( f ) ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /** 
     * \fn void create_vertices( unsigned nverts , double* lverts , VMAP& vmap )
     *
     * \brief Creates the vertices of this surface.
     *
     * \param nverts The number of vertices of this surface.
     * \param lverts An array with the vertex coordinates.
     * \param vmap A hash table with vertex identifiers and pointers.
     */
    void create_vertices(
			 unsigned nverts , 
			 double* lverts , 
			 VMAP& vmap 
			) ;


    /**
     * \fn void create_faces( unsigned nfaces , unsigned* lfaces , HMAP& hmap , VMAP& vmap ) 
     * 
     * \brief Creates the faces and half-edges of this surface.
     *
     * \param nfaces The number of faces of this surface.
     * \param lfaces An array with the vertex identifiers of each face
     * of this surface.
     * \param  hmap  A  hash  table  with  half-edge  identifiers  and
     * pointers.
     * \param vmap A hash table with vertex identifiers and pointers.
     */
    void create_faces(
		      unsigned nfaces , 
		      unsigned* lfaces , 
		      HMAP& hmap , 
		      VMAP& vmap 
		     ) ;

    /**
     * \fn void create_edges( HMAP& hmap ) 
     *
     * \brief Creates the edges of this surface.
     *
     * \param  hmap  A  hash  table  with  half-edge  identifiers  and
     * pointers.
     */
    void create_edges( HMAP& hmap ) ;


    /** 
     * \fn bool is_consistent() const
     *
     * \brief  Performs  a   topological  consistency  check  in  this
     * surface.
     *
     * \return  A logic value  true if  this surface  is topologically
     * consistent and false otherwise.
     */
    bool is_consistent() const ;


  protected:

    // ---------------------------------------------------------------
    //
    // Protected data members
    //
    // ---------------------------------------------------------------

    std::list< Vertex* > _list_of_vertices ;  ///< List of vertices of this surface.


    std::list< Edge* > _list_of_edges ;  ///< List of edges of this surface.


    std::list< Face* > _list_of_faces ;  ///< List of faces of this surface.

  } ;


  /**
   * \fn Surface< VAttrib , FAttrib , EAttrib , HAttrib >::Surface( unsigned nverts , double* lverts , unsigned nfaces , unsigned* lfaces )
   *
   * \brief Creates an instance of this class.
   *
   * \param nverts The number of vertices of this surface.
   * \param lverts An array with the vertex coordinates.
   * \param nfaces The number of faces of this surface.
   * \param lfaces An  array with the vertex identifiers  of each face
   * of this surface.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib 
           >
  Surface< VAttrib , FAttrib , EAttrib , HAttrib >::Surface(
     unsigned nverts ,
     double* lverts ,
     unsigned nfaces ,
     unsigned* lfaces 
    )
  {
    //
    // Create the vertices of the surface.
    //

    /** Define a hash table to keep track of vertex IDs. */
    VMAP vmap ;

    /** Create the vertices and initialize the table. */
    create_vertices( nverts , lverts , vmap ) ;

    //
    // Create the faces and halfedges of the surface.
    //

    /** Define a hash table to keep track of the halfedge IDs. */
    HMAP hmap ;

    /** Create the faces and halfedges and initialize the table. */
    create_faces( nfaces , lfaces , hmap , vmap ) ;

    /** Release memory associated with the vertex ID table. */
    vmap.clear();

    //
    // Create the edges of the surface.
    //
    create_edges( hmap ) ;

    /** Release memory associated with the halfedge ID table. */
    hmap.clear();

    /** Check if this surface is topologically consistent. */
    assert( is_consistent() ) ;

    return ;
  }


  /**
   * \fn Surface< VAttrib , FAttrib , EAttrib , HAttrib >::~Surface()
   *
   * \brief Releases the memory held by an instance of this class.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib 
           >
  Surface< VAttrib , FAttrib , EAttrib , HAttrib >::~Surface()
  {
    /**
     * Releases memory associated with the vertices.
     */
    for( VertexIterator vit = vertices_begin() ; vit != vertices_end() ; ++vit ) {
      Vertex* vertex = *vit ;

      if ( vertex != 0 ) delete vertex ;
    }

    _list_of_vertices.clear() ;

    /**
     * Releases memory associated with the edges and half-edges.
     */
    for( EdgeIterator eit = edges_begin() ; eit != edges_end() ; ++eit ) {
      Edge* edge = *eit ;

      Halfedge* h1 = edge->get_first_halfedge() ;
      Halfedge* h2 = edge->get_second_halfedge() ;

      if ( h1 != 0 ) delete h1 ;
      if ( h2 != 0 ) delete h2 ;

      if ( edge != 0 ) delete edge ;
    }

    _list_of_edges.clear() ;

    /**
     * Releases memory associated with the faces.
     */
    for( FaceIterator fit = faces_begin() ; fit != faces_end() ; ++fit ) {
      Face* face = *fit ;

      if ( face != 0 ) delete face ;
    }

    _list_of_faces.clear() ;

    return;
  }


  /** 
   * \fn void \fn Surface< VAttrib , FAttrib , EAttrib , HAttrib >::create_vertices( unsigned nverts , double* lverts , VMAP& vmap )
   *
   * \brief Creates the vertices of this surface.
   *
   * \param nverts The number of vertices of this surface.
   * \param lverts An array with the vertex coordinates.
   * \param vmap A hash table with vertex identifiers and pointers.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib 
           >
  void
  Surface< VAttrib , FAttrib , EAttrib , HAttrib >::create_vertices(
     unsigned nverts ,
     double* lverts ,
     VMAP& vmap
    )
  {
    /**
     * For each set of coordinates, ( x  , y , z ), in the given array
     * of  coordinates,  create  a   vertex  element  with  the  given
     * coordinates.
     */
    for ( unsigned i = 0 ; i < nverts ; i++ ) {
      /** Get the index of the first coordinate of the i-th vertex. */
      const unsigned j = 3 * i;

      /** Allocate memory for the vertex. */
      Vertex* vertex = new Vertex( 
				  lverts[ j     ] ,
				  lverts[ j + 1 ] ,
				  lverts[ j + 2 ] ,
				  0
				 ) ;

      /** Insert the vertex and its identifier into a hash table. */
      vmap.insert( std::make_pair( i , vertex ) ) ;
    }

    return ;
  }


  /**
   * \fn void \fn Surface< VAttrib , FAttrib , EAttrib , HAttrib >::create_faces( unsigned nfaces , unsigned* lfaces , HMAP& hmap , VMAP& vmap ) 
   * 
   * \brief Creates the faces and half-edges of this surface.
   *
   * \param nfaces The number of faces of this surface.
   * \param lfaces An  array with the vertex identifiers  of each face
   * of this surface.
   * \param hmap A hash table with half-edge identifiers and pointers.
   * \param vmap A hash table with vertex identifiers and pointers.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib 
           >
  void
  Surface< VAttrib , FAttrib , EAttrib , HAttrib >::create_faces(
     unsigned nfaces ,
     unsigned* lfaces ,
     HMAP& hmap ,
     VMAP& vmap
    )
  {
    /**
     * Loop over  the array  of face vertex  identifiers and  for each
     * iteration i, create  the i-th face ( i.e.,  a triangle) of this
     * surface.
     */
    for ( unsigned i = 0 ; i < nfaces ; i++ ) {
      /** Create a face vertex counter. */
      unsigned j = 3 * i ;

      /** Allocate memory for a new face. */
      Face* face = new Face( 0 ) ;

      //
      // Create the three half-edges of the i-th face.
      //

      Halfedge* h[ 3 ] ;

      for ( unsigned k = 0 ; k < 3 ; k++ ) {
        /** 
         * Find the ID of the k-th vertex of the i-th face. 
         */
        typename VMAP::const_iterator vit = vmap.find( lfaces[ j + k ] ) ;

        /** 
         * If k-th vertex cannot be found, throw an exception. 
         */
	assert( vit != vmap.end() ) ;

        Vertex* vertex = vit->second ;

        /** 
         * Allocate memory for the k-th half-edge. 
         */
        h[ k ] = new Halfedge(
                              vertex ,
			      0 ,
			      face ,
			      0 ,
			      0
                             ) ;

        /** 
         * Set the half-edge pointer of the vertex pointed by vertex.
         */
        vertex->set_halfedge( h[ k ] ) ;
      }

      /** 
       * Update the pointers to previous and next half-edges. 
       */
      h[ 0 ]->set_prev( h[ 2 ] ) ;
      h[ 2 ]->set_next( h[ 0 ] ) ;

      h[ 1 ]->set_prev( h[ 0 ] ) ;
      h[ 0 ]->set_next( h[ 1 ] ) ;

      h[ 2 ]->set_prev( h[ 1 ] ) ;
      h[ 1 ]->set_next( h[ 2 ] ) ;

      /** 
       * Set the pointer to the first half-edge of the i-th face.
       */
      face->set_halfedge( h[ 0 ] ) ;


      /** 
       * Add information to the table of halfedge IDs. 
       */
      for ( unsigned k = 0 ; k < 3 ; k++ ) {
        const unsigned v1 = lfaces[ j +     k             ] ;
        const unsigned v2 = lfaces[ j + ( ( k + 1 ) % 3 ) ] ;

        hmap.insert( std::make_pair( std::make_pair( v1 , v2 ) , h[ k ] ) ) ;
      }

      /** 
       * Insert face into the list of faces of this surface.
       */
      _list_of_faces.push_back( face ) ;
    }

    /** 
     * Place in  the list of vertices  only the vertices  that are the
     * origin  of some  half-edge. The  ones  that are  not should  be
     * ignored.
     */
    for ( typename VMAP::iterator vit = vmap.begin() ; 
	  vit != vmap.end() ; vit++ ) {
      if ( vit->second->get_halfedge() != 0 ) {
        _list_of_vertices.push_back( vit->second ) ;
      }
      else {
	delete vit->second ;
	vit->second = 0 ;
      }
    }

    return ;
  }


  /**
   * \fn void Surface< VAttrib , FAttrib , EAttrib , HAttrib >::create_edges( HMAP& hmap ) 
   *
   * \brief Creates the edges of this surface.
   *
   * \param hmap A hash table with half-edge identifiers and pointers.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib 
           >
  void
  Surface< VAttrib , FAttrib , EAttrib , HAttrib >::create_edges( 
     HMAP& hmap 
    )
  {
    /**
     * Loop over the  set of half-edges of this  surface and, for each
     * of  them, create  an edge  (if the  edge has  not  been created
     * already)  and   update  their  edge   and  half-edge  pointers,
     * respectively.
     */
    for ( typename HMAP::iterator hit = hmap.begin(); 
	  hit != hmap.end(); ++hit ) {
      /**
       * If the  edge pointer  of the current  half-edge has  not been
       * assigned an edge yet, do  it now. Otherwise, skip to the next
       * iteration.
       */
      if ( hit->second->get_edge() == 0 ) {
        /**
         * Create a new edge and  initializes its pointer to its first
         * half-edge.
         */
	Edge* edge = new Edge( hit->second , 0 ) ;

        /** 
         * Initializes the edge pointer of the current half-edge. 
         */
        hit->second->set_edge( edge ) ;

        /** 
         * Insert the  newly created  edge into the  list of  edges of
         * this surface.
         */
        _list_of_edges.push_back( edge ) ;

        /** 
         * Get the half-edge vertex identifiers.
         */
        unsigned v1 = hit->first.first ;
        unsigned v2 = hit->first.second ;

        /** 
	 * Initialize the edge pointer of the second half-edge.
	 */
        typename HMAP::iterator ait = 
	  hmap.find( std::make_pair( v2 , v1 ) ) ;

        if ( ait != hmap.end() ) {
          ait->second->set_edge( edge ) ;
          edge->set_second_halfedge( ait->second ) ;
        }

      }
    }

    return ;
  }


  /** 
   * \fn bool \fn Surface< VAttrib , FAttrib , EAttrib , HAttrib >::is_consistent() const
   *
   * \brief Performs a topological consistency check in this surface.
   *
   * \return  A logic  value  true if  this  surface is  topologically
   * consistent and false otherwise.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib
           >
  bool
  Surface< VAttrib , FAttrib , EAttrib , HAttrib >::is_consistent() const
  {
    //
    // Check half-edge and face pointers. 
    //

    for ( FaceIterator fit = faces_begin(); fit != faces_end() ; ++fit ) {
      /** 
       * Get  the pointer  to  the  first half-edge  of  the cycle  of
       * half-edges of the face.
       */
      Halfedge* h1 = ( *fit )->get_halfedge() ;

      /** 
       * Pointer to the first half-edge cannot be null. 
       */
      if ( h1 == 0 ) return false ;

      /**
       * Get the other two half-edges of the face. 
       */
      Halfedge* h2 = h1->get_next() ;
      Halfedge* h3 = h1->get_prev() ;

      /** 
       * Pointers  to  the  next  and previous  half-edges  cannot  be
       * null. 
       */
      if ( ( h2 == 0 ) || ( h3 == 0 ) ) return false ;

      /**
       * Pointers to face cannot be null. 
       */
      if ( ( h1->get_face() == 0 ) || ( h2->get_face() == 0 ) ||
	   ( h3->get_face() == 0 ) )
      {
	return false ;
      }

      /**
       * Pointer to face must be the same.
       */
      if ( ( h1->get_face() != h2->get_face() ) ||
	   ( h1->get_face() != h3->get_face() ) ||
	   ( h2->get_face() != h3->get_face() ) )
      {
	return false ;
      }

      /**
       * Pointers to the previous and next must be consistent.
       */
      if ( h2->get_prev() != h1 ) return false ;

      if ( h3->get_next() != h1 ) return false ;

      if ( h2->get_next() != h3 ) return false ;

      if ( h3->get_prev() != h2 ) return false ;

      /**
       * Pointers to origin vertices cannot be null.
       */
      if ( h1->get_origin() == 0 ) return false ;

      if ( h2->get_origin() == 0 ) return false ;

      if ( h3->get_origin() == 0 ) return false ;

      /**
       * Pointers to edges cannot be null.
       */
      if ( h1->get_edge() == 0 ) return false ;

      if ( h2->get_edge() == 0 ) return false ;

      if ( h3->get_edge() == 0 ) return false ;

      /**
       * Check edge pointers consistency.
       */
      if ( 
	  ( h1 != h1->get_edge()->get_first_halfedge()  ) &&
	  ( h1 != h1->get_edge()->get_second_halfedge() ) 
	 )
      {
	return false ;
      }

      if ( 
	  ( h2 != h2->get_edge()->get_first_halfedge()  ) &&
	  ( h2 != h2->get_edge()->get_second_halfedge() ) 
	 )
      {
	return false ;
      }

      if ( 
	  ( h3 != h3->get_edge()->get_first_halfedge()  ) &&
	  ( h3 != h3->get_edge()->get_second_halfedge() ) 
	 )
      {
	return false ;
      }

      /**
       * Check mate half-edge pointers.
       */
      if ( h1->get_mate() == 0 ) return false ;

      if ( h2->get_mate() == 0 ) return false ;

      if ( h3->get_mate() == 0 ) return false ;

      if ( h1->get_mate()->get_mate() != h1 ) return false ;

      if ( h2->get_mate()->get_mate() != h2 ) return false ;

      if ( h3->get_mate()->get_mate() != h3 ) return false ;

    } 
    // end of the for loop for checking faces.

    //
    // Check vertex pointers.
    //

    for ( VertexIterator vit = vertices_begin() ; vit != vertices_end() ; ++vit ) {
      /** 
       * Get  the pointer  to one  of  the half-edges  of the  current
       * vertex.
       */ 
      Halfedge* h1 = ( *vit )->get_halfedge() ;

      /**
       * The half-edge pointer cannot be null.
       */
      if ( h1 == 0 ) return false ;

      /** 
       * The half-edge  pointer to the  origin vertex must  agree with
       * the pointer to the current vertex.
       */
      if ( h1->get_origin() != *vit ) return false ;

      /** 
       * Get the  next half-edge in a counterclowise  traversal of the
       * edges incident to the current vertex.
       */

      Halfedge* h2 = h1->get_prev()->get_mate() ;
      do {
	/** 
	 * The half-edge pointer to  the origin vertex must agree with
	 * the pointer to the current vertex.
	 */
	if ( h2->get_origin() != *vit ) {
	  return false;
	}
	else {
	  h2 = h2->get_prev()->get_mate() ;
        }
      }
      while ( h2 != h1) ;
    }
    // end of the for loop for checking vertices

    /**
     * Check the edge pointers.
     */
    for( EdgeIterator eit = edges_begin(); eit != edges_end() ; ++eit ) {
      /** 
       * Get  the pointer  to the first half-edge of the edge.
       */ 
      Halfedge* h1 = (*eit)->get_first_halfedge() ;

      /** 
       * The pointer to the first half-edge cannot be null.
       */
      if ( h1 == 0 ) return false ;

      /**
       * Check the edge pointer of the first half-edge.
       */
      if ( h1->get_edge() != *eit ) return false ;


      /** 
       * Get the pointer to the second half-edge of the edge.
       */ 
      Halfedge* h2 = (*eit)->get_second_halfedge() ;

      /** 
       * The pointer to the second half-edge cannot be null.
       */
      if ( h2 == 0 ) return false ;

      /**
       * Check the edge pointer of the first half-edge.
       */
      if ( h2->get_edge() != *eit ) return false ;

      /**
       * The first and second half-edges cannot be the same.
       */
      if ( h1 == h2 ) return false ;

      /**
       * The first and second half-edges must be mates of each other.
       */
      if ( ( h1 != h2->get_mate() ) || ( h2 != h1->get_mate() ) )
      {
	return false ;
      }
    }
    // end of the for loop for checking edges

    return true;
  }

}

/** @} */ //end of group class.

#endif   // SURFACE_H
