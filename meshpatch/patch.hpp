/** 
 * \file patch.hpp
 *  
 * \brief  Definition and  implementation  of the  class Patch,  which
 * represents  a triangulation of  a surface  patch homeomorphic  to a
 * disk.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date October 2010
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef PATCH_HPP
#define PATCH_HPP

#include "patchvertex.hpp"   // meshpatch::PatchVertex
#include "patchedge.hpp"     // meshpatch::PatchEdge
#include "mypair.hpp"        // meshpatch::MyPair

#include <cassert>           // assert
#include <vector>            // std::vector
#include <list>              // std::list
#include <map>               // std::map

#include <iostream>

/** 
 * \defgroup MESHPATCHNameSpace Namespace meshpatch.
 * @{
 */

/**
 * \namespace meshpatch
 *
 * \brief   The  namespace  meshpatch   contains  the   definition  and
 * implementation of  all classes of the data  structure MeshPatch data
 * structure,  which can  be used  to  represent a  triangulation of  a
 * surface patch homeomorphic  to a disk.
 */


// --------------------------------------------------------------------
//
// The MeshPatch data  structure is based on the  "directed edge" data
// structure. See  paper "Directed Edges --  A Scalable Representation
// for  Triangular Meshes"  by  S.  Campagna,  L.  Kobbelt, and  H.-P.
// Seidel at the  Journal of graphics Tools, Vol. 3,  No. 4, pp. 1-12,
// 1998.
//
// --------------------------------------------------------------------


namespace meshpatch {


  /**
   * \class Patch
   *
   * \brief This  class represents a triangulation of  a surface patch
   * homeomorphic to a closed disk in \f$R^2\f$.
   */
  template < 
             typename VAttrib = int ,
             typename EAttrib = int ,
             typename PAttrib = int 
           >
  class Patch {
  public:
    // ---------------------------------------------------------------
    //
    // Type definitions
    //
    // ---------------------------------------------------------------

    /**
     * \typedef PatchVertex
     *
     * \brief     Defines    PatchVertex     as    an     alias    for
     * meshpatch::PatchVertex< VAttrib, EAttrib , PAttrib >.
     */    
    typedef meshpatch::PatchVertex< VAttrib, EAttrib , PAttrib > PatchVertex ;

    /**
     * \typedef PatchEdge
     *
     * \brief Defines PatchEdge  as an alias for meshpatch::PatchEdge<
     * VAttrib, FAttrib , EAttrib , HAttrib >.
     */    
    typedef meshpatch::PatchEdge< VAttrib, EAttrib , PAttrib > PatchEdge ;


    /**
     * \typedef Pair
     *
     * \brief Defines Pair as an alias for meshpatch::MyPair< unsigned
     * , unsigned >.
     */    
    typedef meshpatch::MyPair< unsigned , unsigned > Pair ;


    /** 
     * \typedef HMAP
     *
     * \brief  Defines  a  hash  table of  half-edge  identifiers  and
     * pointers.
     */
    typedef typename std::map< Pair , unsigned > HMAP ;
 

    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Patch( unsigned nverts , unsigned nfaces , unsigned* lfaces , double* vcoords )
     *
     * \brief Creates an instance of this class.
     *
     * \param nverts The number of vertices of this surface patch.
     * \param nfaces The number of faces of this surface patch.     
     * \param lfaces  An  array with  the vertex  identifiers  of each
     * face.
     * \param vcoords The vertex coordinates.
     */
    Patch(
	  unsigned nverts ,
	  unsigned nfaces ,
	  unsigned* lfaces ,
          double* vcoords
	 ) ;


    /**
     * \fn Patch( const Patch< VAttrib , EAttrib , PAttrib >& patch  )
     *
     * \brief Creates an instance of this class.
     *
     * \param patch A pointer to an object from this class.
     */
    Patch( const Patch< VAttrib , EAttrib , PAttrib >& patch ) ;


    /**
     * \fn ~Patch()
     *
     * \brief Releases the memory held by an instance of this class.
     */
    virtual ~Patch() 
    {}


    /** 
     * \fn inline unsigned get_number_of_vertices() const
     *
     * \brief Returns the number of vertices of this surface patch.
     *
     * \return The number of vertices of this surface patch.
     */
    inline unsigned get_number_of_vertices() const
    {
      return _patch_vts.size() ;
    }


    /** 
     * \fn inline unsigned get_number_of_edges() const
     *
     * \brief Returns the number of edges of this surface patch.
     *
     * \return The number of edges of this surface patch.
     *
     */
    inline unsigned get_number_of_edges() const
    { 
      return ( _patch_egs.size() + _patch_bds.size() ) >> 1 ; 
    }


    /** 
     * \fn inline unsigned get_number_of_halfedges() const
     *
     * \brief Returns the number of half-edges of this surface patch.
     *
     * \return The number of half-edges of this surface patch.
     *
     */
    inline unsigned get_number_of_halfedges() const
    { 
      return _patch_egs.size() ; 
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
      return _patch_egs.size() / 3 ; 
    }


    /** 
     * \fn inline unsigned get_number_of_boundary_edges() const
     *
     * \brief  Returns the number  of boundary  edges of  this surface
     * patch.
     *
     * \return The number of boundary edges of this surface patch.
     *
     */
    inline unsigned get_number_of_boundary_edges() const
    { 
      return _patch_bds.size() ; 
    }


    /**
     * \fn inline unsigned get_vertex_edge( unsigned i ) const
     *
     * \brief  Returns the  index of  an edge  incident with  the i-th
     * vertex.
     *
     * \param i The index of a surface patch vertex.
     *
     * \return The  index of  a surface patch  edge incident  with the
     * i-th vertex.
     */
    inline unsigned get_vertex_edge( unsigned i ) const
    {
      assert( i < _patch_vts.size() ) ;

      return _patch_vts[ i ].get_edge() ;
    }


    /**
     * \fn inline void get_vertex_coords( unsigned i , double& x , double& y , double& z ) const
     *
     * \brief  Returns the  index of  an edge  incident with  the i-th
     * vertex.
     *
     * \param i The index of a surface patch vertex.
     * \param x First Cartesian coordinate of the i-th vertex.
     * \param y Second Cartesian coordinate of the i-th vertex.
     * \param z Third Cartesian coordinate of the i-th vertex.
     *
     * \return The  index of  a surface patch  edge incident  with the
     * i-th vertex.
     */
    inline void get_vertex_coords( unsigned i , double& x , double& y , double& z  ) const
    {
      assert( i < _patch_vts.size() ) ;

      x = _patch_vts[ i ].x() ;
      y = _patch_vts[ i ].y() ;
      z = _patch_vts[ i ].z() ;

      return ;
    }


    /**
     * \fn inline VAttrib& get_vertex_attributes( unsigned i )
     *
     * \brief Returns the attributes associated with the i-th vertex.
     *
     * \param i The index of a surface patch vertex.
     *
     * \return The attributes associated with the i-th vertex.
     */
    inline VAttrib& get_vertex_attributes( unsigned i )
    {
      assert( i < _patch_vts.size() ) ;

      return _patch_vts[ i ].get_attributes() ;
    }


    /**
     * \fn inline unsigned get_origin_vertex( unsigned i ) const
     *
     * \brief Returns the index of  the origin vertex of the i-th edge
     * of this surface patch.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The index of the origin vertex of the i-th edge of this
     * surface patch.
     */
    inline unsigned get_origin_vertex( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;

      return _patch_egs[ i ].get_origin() ;
    }


    /**
     * \fn inline unsigned get_mate( unsigned i ) const
     *
     * \brief Returns the  index of the mate edge of  the i-th edge of
     * this surface patch.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The  index of the  mate edge of  the i-th edge  of this
     * surface patch.
     */
    inline unsigned get_mate( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;
      assert( !is_boundary_edge( i ) ) ;

      return _patch_egs[ i ].get_mate() ;
    }


    /**
     * \fn inline unsigned get_next( unsigned i ) const
     *
     * \brief Returns the index of the edge following the i-th edge of
     * this surface patch in the same face cycle.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The index  of the edge following the  i-th edge of this
     * surface patch in the same face cycle.
     */
    inline unsigned get_next( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;

      return _patch_egs[ _patch_egs[ i ].get_prev() ].get_prev() ;  
    }


    /**
     * \fn inline unsigned get_prev( unsigned i ) const
     *
     * \brief Returns the index of the edge preceding the i-th edge of
     * this surface patch in the same face cycle.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The index  of the edge preceding the  i-th edge of this
     * surface patch in the same face cycle.
     */
    inline unsigned get_prev( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;

      return _patch_egs[ i ].get_prev() ;  
    }


    /**
     * \fn inline EAttrib& get_edge_attributes( unsigned i ) const
     *
     * \brief Returns  the set of attributes associated  with the i-th
     * edge.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The set of attributes associated with the i-th edge.
     */
    inline EAttrib& get_edge_attributes( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;

      return _patch_egs[ i ].get_attributes() ;  
    }


    /**
     * \fn inline unsigned get_boundary_id( unsigned i ) const
     *
     * \brief Returns  the index of the  i-th edge in  the sequence of
     * boundary edges.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The index of the i-th edge of this surface patch in the
     * sequence of boundary edges.
     */
    inline unsigned get_boundary_id( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;
      assert( is_boundary_edge( i ) ) ;

      return _patch_egs[ i ].get_mate() ;
    }


    /**
     * \fn inline unsigned get_ith_boundary_edge( unsigned i ) const
     *
     * \brief Returns the index of the i-th boundary edge.
     *
     * \param i The index of a surface patch boundary edge.
     *
     * \return The  index of  the i-th boundary  edge of  this surface
     * patch.
     */
    inline unsigned get_ith_boundary_edge( unsigned i ) const
    {
      assert( i < _patch_bds.size() ) ;

      return _patch_bds[ i ] ;  
    }


    /**
     * \fn void get_vertex_star( unsigned i , std::vector< unsigned >& star ) const
     *
     * \brief Gets a list with the indices of all edges incident to
     * the i-th vertex.
     *
     * \param i The index of a surface patch vertex.
     * \param star A reference to a vector of edge indices.
     *
     */
    void get_vertex_star( unsigned i , std::vector< unsigned >& star ) const
    {
      assert( i < _patch_vts.size() ) ;

      unsigned eg1 = _patch_vts[ i ].get_edge() ;
      unsigned eg2 = eg1 ;
      bool endloop = false ;
      do {
	star.push_back( eg2 ) ;
	eg2 = get_prev( eg2 ) ;
	if ( is_boundary_edge( eg2 ) ) {
	  star.push_back( eg2 ) ;
	  endloop = true ;
	}
	else {
	  eg2 = get_mate( eg2 ) ;
	}
      }
      while ( ( eg2 != eg1 ) && !endloop ) ;

      return ;  
    }


    /**
     * \fn inline PAttrib& get_attributes()
     *
     * \brief  Returns  the set  of  attributes  associated with  this
     * surface patch.
     *
     * \return  The set  of  attributes associated  with this  surface
     * patch.
     */
    inline PAttrib& get_attributes()
    {
      return _attributes ;  
    }


    /**
     * \fn inline void set_attributes( PAttrib& attrib )
     *
     * \brief Assigns a set of attributes with this surface patch.
     *
     * \param attrib A set of attributes.
     */
    inline void set_attributes( PAttrib& attrib )
    {
      _attributes = attrib ;  
    }


    /**
     * \fn inline bool is_boundary_edge( unsigned i ) const
     *
     * \brief Decides whether the i-th edge of this surface patch is a
     * boundary edge.
     *
     * \param i The index of a surface patch edge.
     *
     * \return The  logic value  true if the  i-th edge is  a boundary
     * edge, and the logic value false otherwise.
     */
    inline bool is_boundary_edge( unsigned i ) const
    {
      assert( i < _patch_egs.size() ) ;

      return _patch_egs[ i ].is_boundary_edge() ;  
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn inline void set_vertex_edge( unsigned i , unsigned j )
     *
     * \brief Assigns an edge to the i-th vertex of this surface patch.
     *
     * \param i The index of a surface patch vertex.
     * \param j The index of a surface patch edge.
     *
     */
    inline void set_vertex_edge( unsigned i , unsigned j )
    {
      assert( i < _patch_vts.size() ) ;
      assert( j < _patch_egs.size() ) ;

      _patch_vts[ i ].set_edge( j ) ;

      return ;
    }


    /**
     * \fn inline void set_mate( unsigned i , unsigned j )
     *
     * \brief Assigns a mate to the i-th edge of this surface patch.
     *
     * \param i The index of a surface patch edge.
     * \param j The index of another surface patch edge (the mate).
     *
     */
    inline void set_mate( unsigned i , unsigned j )
    {
      assert( i < _patch_egs.size() ) ;
      assert( j < _patch_egs.size() ) ;

      _patch_egs[ i ].set_mate( j ) ;
      _patch_egs[ j ].set_mate( i ) ;

      _patch_egs[ i ].set_boundary_status( false ) ;
      _patch_egs[ j ].set_boundary_status( false ) ;

      return ;
    }


    /**
     * \fn void create_edges( unsigned nfaces , unsigned* lfaces , HMAP& hmap ) 
     * 
     * \brief Creates the edges of this surface patch.
     *
     * \param nfaces The number of faces of the surface patch.
     * \param lfaces  The set  of vertex indices  of each face  of the
     * surface patch.
     * \param hmap A hash table with vertex and edge identifiers.
     *
     */
    void create_edges( unsigned nfaces , unsigned* lfaces , 
      HMAP& hmap ) ;


    /**
     * \fn void identify_edges( const HMAP& hmap ) 
     *
     * \brief Identifies mate edges and boundary edges.
     *
     * \param hmap A hash table with vertex and edge identifiers.
     */
    void identify_edges( const HMAP& hmap ) ;


  protected:

    // ---------------------------------------------------------------
    //
    // Protected data members
    //
    // ---------------------------------------------------------------

    std::vector< PatchVertex > _patch_vts ;  ///< Sequence of vertices of this surface patch.

    std::vector< PatchEdge > _patch_egs ;  ///< Sequence of edges of this surface patch.

    std::vector< unsigned > _patch_bds ;  ///< Sequence of boundary edges of this surface patch.

    PAttrib _attributes ;  ///< Set of attributes associated with this patch.

  } ;


  /**
   * \fn Patch< VAttrib , EAttrib , PAttrib >::Patch( unsigned nverts , unsigned nfaces , unsigned* lfaces , double* vcoords )
   *
   * \brief Creates an instance of this class.
   *
   * \param nverts The number of vertices of this surface patch.
   * \param nfaces The number of faces of this surface patch.     
   * \param  lfaces An  array with  the vertex  identifiers  of each
   * face.
   * \param vcoords The vertex coordinates.
   */
  template < 
            typename VAttrib ,
            typename EAttrib ,
            typename PAttrib 
           >
  Patch< VAttrib , EAttrib , PAttrib >::Patch(
					      unsigned nverts ,
					      unsigned nfaces ,
					      unsigned* lfaces ,
					      double* vcoords
					     )
  {
    //
    // Create the vertices of the surface.
    //
    assert( nverts >= 3 ) ;
    assert( nfaces >= 1 ) ;

    _patch_vts.resize( nverts ) ;

    //
    // Set the values of the vertex coordinates.
    //
    for ( unsigned i = 0 ; i < nverts ; i++ ) {
      _patch_vts[ i ].set_x_coord( vcoords[ 3 * i     ] ) ;
      _patch_vts[ i ].set_y_coord( vcoords[ 3 * i + 1 ] ) ;
      _patch_vts[ i ].set_z_coord( vcoords[ 3 * i + 2 ] ) ;
    }
    
    //
    // Create the edges of this surface patch.
    //

    /* Define a hash table to keep track of the half-edge IDs. */
    HMAP hmap ;

    create_edges( nfaces , lfaces , hmap ) ;

    //
    // Identify mate edges.
    //
    identify_edges( hmap ) ;

    return ;
  }


  /**
   * \fn Patch< VAttrib , EAttrib , PAttrib >::Patch( const Patch< VAttrib , EAttrib , PAttrib >& patch  )
   *
   * \brief Creates an instance of this class.
   *
   * \param patch A reference to an object to be cloned.
   */
  template < 
            typename VAttrib ,
            typename EAttrib ,
            typename PAttrib 
           >
  Patch< VAttrib , EAttrib , PAttrib >::Patch( const Patch< VAttrib , EAttrib , PAttrib >& patch )
  {
    assert( patch.get_number_of_vertices() >= 3 ) ;
    assert( patch.get_number_of_faces() >= 1 ) ;

    _patch_vts = patch._patch_vts ;
    _patch_egs = patch._patch_egs ;
    _patch_bds = patch._patch_bds ;

    this->_attributes = patch._attributes ;

    return ;
  }


  /**
   * \fn void Patch< VAttrib , EAttrib , PAttrib >::create_edges( unsigned nfaces , unsigned* lfaces , HMAP& hmap ) 
   * 
   * \brief Creates the edges of this surface patch.
   *
   * \param nfaces The number of faces of the surface patch.
   * \param  lfaces The  set of  vertex indices  of each  face  of the
   * surface patch.
   * \param hmap A hash table with vertex and edge identifiers.
   *
   */
  template < 
             typename VAttrib ,
             typename EAttrib ,
             typename PAttrib 
           >
  void
  Patch< VAttrib , EAttrib , PAttrib >::create_edges(
						     unsigned nfaces ,
						     unsigned* lfaces ,
						     HMAP& hmap 
						    )
  {
    /* Define a counter for the number of edges. */
    unsigned j = 0 ;

    /*
     * For each face, create its three edges, and assign an ID to each
     * of them.
     */
    for ( unsigned i = 0 ; i < nfaces ; i++ ) {
      unsigned v1 = lfaces[ 3 * i     ] ;
      unsigned v2 = lfaces[ 3 * i + 1 ] ;
      unsigned v3 = lfaces[ 3 * i + 2 ] ;

      _patch_egs.push_back( PatchEdge( v1 , 0 , 3 * i + 2 ) ) ;
      _patch_egs.push_back( PatchEdge( v2 , 0 , 3 * i     ) ) ;
      _patch_egs.push_back( PatchEdge( v3 , 0 , 3 * i + 1 ) ) ;

      assert( hmap.find( Pair( v1 , v2 ) ) == hmap.end() ) ;
      assert( hmap.find( Pair( v2 , v3 ) ) == hmap.end() ) ;
      assert( hmap.find( Pair( v3 , v1 ) ) == hmap.end() ) ;
						     
      hmap.insert( std::make_pair( Pair( v1 , v2 ) , 3 * i     ) ) ;
      hmap.insert( std::make_pair( Pair( v2 , v3 ) , 3 * i + 1 ) ) ;
      hmap.insert( std::make_pair( Pair( v3 , v1 ) , 3 * i + 2 ) ) ;

      set_vertex_edge( v1 , 3 * i     ) ;
      set_vertex_edge( v2 , 3 * i + 1 ) ;
      set_vertex_edge( v3 , 3 * i + 2 ) ;

      j += 3 ;
    }

    return ;
  }


  /**
   * \fn void Patch< VAttrib , EAttrib , PAttrib >::identify_edges( const HMAP& hmap ) 
   *
   * \brief Identifies mate edges and boundary edges.
   *
   * \param hmap A hash table with vertex and edge identifiers.
   */
  template < 
             typename VAttrib ,
             typename EAttrib ,
             typename PAttrib 
           >
  void
  Patch< VAttrib , EAttrib , PAttrib >::identify_edges( const HMAP& hmap ) 
  {
    /*
     * Find the mate of each edge, and determine the boundary edges.
     */

    /* 
     * Define a variable for storing the index of one boundary edge.
     */
    unsigned last_id = 0 ;

    /* Define a counter for boundary edges. */
    unsigned j = 0 ;

    typename HMAP::const_iterator hit ;
    for ( hit = hmap.begin() ; hit != hmap.end() ; ++hit ) {
      unsigned v1 = ( *hit ).first.first() ;
      unsigned v2 = ( *hit ).first.second() ;
      unsigned id = ( *hit ).second ;

      /* Does the edge have a mate? */
      typename HMAP::const_iterator hitaux = 
	hmap.find( Pair( v2 , v1 ) ) ;

      if ( hitaux != hmap.end() ) {
	/* Indentify mate edges. */
	set_mate( id , ( *hitaux ).second ) ;
      }
      else {
	/* Set boundary edge status. */
	_patch_egs[ id ].set_boundary_status( true ) ; 

	/* 
	 * Make  sure the  edge  pointer  of the  origin  vertex of  a
	 * boundary edge points to the  boundary edge. If so, the code
	 * for computing the star of  a vertex can be more efficiently
	 * written.
	 */
	set_vertex_edge( v1 , id ) ;

	/* Keep the index of the boundary edge. */
	last_id = id ;

	/* Update the boundary edge counter. */
	++j ;
      }      
    }

    /*
     * If there  is any boundary edge,  we place them  in the boundary
     * edge array in  the same order they occur  in a counterclockwise
     * traversal of  the boundary.  The boundary  of the triangulation
     * is  expected  to be  a  simple,  closed  curve. Otherwise,  the
     * following loop will detect the problem and an exception will be
     * thrown.
     */
    if ( j > 0 ) {
      /* 
       * Get all edges that follow  the first one in the same boundary
       * curve. If  the triangulation  is consistent, there  should be
       * only one boundary curve.
       */

      /* 
       * Define a marker for the first boundary edge to be processed. 
       */
      unsigned first_id = last_id ;

      /* Define a boundary edge counter. */
      unsigned k = 0 ;

      do {
	/* Insert the edge index in the boundary edge array. */
	_patch_bds.push_back( last_id ) ;
	_patch_egs[ last_id ].set_mate( k ) ;

	/* Increment the boundary edge counter. */
	++k ;

	/* Get the next boundary edge. */
	unsigned next_id = get_next( last_id )  ;
	while ( !is_boundary_edge( next_id ) && ( next_id != last_id ) ) {
	  next_id = get_next( get_mate( next_id ) ) ;
	}

	/* 
	 * Make sure we did find a boundary edge other than "last_id".
	 */
	assert( next_id != last_id ) ;

	/* Get the index of the next boundary edge. */
	last_id = next_id ;
      }
      while ( last_id != first_id ) ;

      /* 
       * To  make sure  we only  have a  simple, closed  curve  as the
       * boundary of  the triangulation, we compare the  values of "j"
       * and "k".
       */
      assert( j == k ) ;
    }

    return ;
  }


}

/** @} */ //end of group class.

#endif   // PATCH_HPP
