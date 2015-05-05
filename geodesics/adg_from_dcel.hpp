/** 
 * \file Dm_from_dcel.hpp
 *  
 * \brief  Definition  and implementation  of  the class  Dm_from_dcel,
 * which implements an abstract interface for a dense mesh represented
 * by the DCEL data structure and that can be parametrized over a base
 * mesh.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * mfsiqueira at gmail (dot) com
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

#ifndef ADG_FROM_DCEL_HPP
#define	ADG_FROM_DCEL_HPP

#include "dcel/surface.h"                       // dcel::Surface
#include "adg_mesh_interface.hpp"              // dm2bm::Dense_mesh_interface
#include "types.h"                            // Vertex and Face attributes

#include <cassert>                           // assert


/** 
 * \defgroup DM2BMPFROMDCELNameSpace Namespace dm2bmfromdcel.
 * @{
 */

/**
 * \namespace dm2bmfromdcel
 *
 * \brief  The  namespace dm2bmfromdcel  contains  the definition  and
 * implementation of two classes  representing triangles meshes in the
 * DCEL data  structure. One  class is supposed  to represent  a dense
 * mesh, while the other is  supposed to represent a coarse mesh. Both
 * classes have  attributes and methods that enable  us to parametrize
 * the dense mesh over the coarse mesh using a discrete geodesic based
 * approach. 
 */

namespace adgfromdcel {


  /**
   * \typedef DMesh
   *
   * \brief Definition of a type name for the surface mesh.
   */
  typedef dcel::Surface<myVertexAttributes,myFaceAttributes,myEdgeAttributes,myHalfedgeAttributes> TMesh ;

 
  /**
   * \class Dm_from_dcel
   *
   * \brief This class represents a  triangle mesh using the DCEL data
   * structure, implements  the abstract interface Adg_mesh_interface,
   * and  contains  some  special   attributes,  which  enable  us  to
   * manipulate triangle meshes that  can be parametrized over coarser
   * meshes.
   */
  class adg_from_dcel : public dg::Adg_mesh_interface< TMesh > {
  public:
    // ---------------------------------------------------------------
    //
    // Type name definitions for the generic mesh components.
    //
    // ---------------------------------------------------------------


    /**
     * \typedef Vertex
     *
     * \brief Definition of a type name for the mesh vertices.
     */
    typedef TMesh::Vertex Vertex ;


    /**
     * \typedef Halfedge
     *
     * \brief Definition of a type name for the mesh half-edges.
     */
    typedef TMesh::Halfedge Halfedge ;


    /**
     * \typedef Edge
     *
     * \brief Definition of a type name for the mesh edges.
     */
    typedef TMesh::Edge Edge ;


    /**
     * \typedef Face
     *
     * \brief Definition of a type name for the mesh faces.
     */
    typedef TMesh::Face Face ;


    /**
     * \typedef VertexIterator
     *
     * \brief Definition of a type name for the vertex iterators.
     */
    typedef TMesh::VertexIterator VertexIterator ;


    /**
     * \typedef EdgeIterator
     *
     * \brief Definition of a type name for the edge iterators.
     */
    typedef TMesh::EdgeIterator EdgeIterator ;


    /**
     * \typedef FaceIterator.
     *
     * \brief Definition of a type name for the face iterators.
     */
    typedef TMesh::FaceIterator FaceIterator ;


    /**
     * \fn Dm_from_dcel( DMesh* m )
     *
     * \brief Creates an instance of this class.
     *
     * \param m Pointer to a mesh stored in a DCEL.
     */    
    adg_from_dcel( TMesh* m ) : _mesh( m )
    {
      assert( _mesh != 0 ) ;
    }


    /**
     * \fn ~Dm_from_dcel()
     *
     * \brief Destroys an instance of this class.
     */
    virtual ~adg_from_dcel()
    {}


    /**
     * \fn inline Vertex* get_org( Halfedge* h ) const
     *
     * \brief Returns  the origin vertex  of a given half-edge  of the
     * mesh.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the origin vertex of a given half-edge.
     */
    inline Vertex* get_org( Halfedge* h )  const 
    {
#ifdef DEBUGMODE
      assert( h != 0 ) ;
//      assert( h->get_face()->get_surface() == _mesh ) ;
#endif

      return h->get_origin() ;
    }


    /**
     * \fn inline Edge* get_edge( Halfedge* h ) const
     *
     * \brief Returns the  edge a given half-edge of  the mesh belongs
     * to.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the edge a given half-edge belongs to.
     */
    inline Edge* get_edge( Halfedge* h ) const 
    {
#ifdef DEBUGMODE
      assert( h != 0 ) ;
   //   assert( h->get_face()->get_surface() == _mesh ) ;
#endif

      return h->get_edge() ;
    }


    /**
     * \fn inline Face* get_face( Halfedge* h ) const
     *
     * \brief Returns the  face a given half-edge of  the mesh belongs
     * to.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the face a given half-edge belongs to.
     */
    inline Face* get_face( Halfedge* h ) const
    {
#ifdef DEBUGMODE
      assert( h != 0 ) ;
      //assert( h->get_face()->get_surface() == _mesh ) ;
#endif

      return h->get_face() ;
    }


    /**
     * \fn inline Halfedge* get_prev( Halfedge* h ) const
     *
     * \brief Returns the half-edge that precedes a given half-edge of
     * the mesh  in the  face cycle of  half-edges that  contains both
     * half-edges.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return  A pointer  to  the half-edge  that  precedes a  given
     * half-edge in their common face half-edge cycle.
     *
     */
    inline Halfedge* get_prev( Halfedge* h ) const 
    {
#ifdef DEBUGMODE
      assert( h != 0 ) ;
      //assert( h->get_face()->get_surface() == _mesh ) ;
#endif

      return h->get_prev() ;
    }


    /**
     * \fn inline Halfedge* get_next( Halfedge* h ) const
     *
     * \brief Returns the half-edge that succeeds a given half-edge of
     * the mesh  in the  face cycle of  half-edges that  contains both
     * half-edges.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return  A pointer  to  the half-edge  that  succeeds a  given
     * half-edge in their common face half-edge cycle.
     */
    inline Halfedge* get_next( Halfedge* h ) const 
    {
#ifdef DEBUGMODE
      assert( h != 0 ) ;
      //assert( h->get_face()->get_surface() == _mesh ) ;
#endif

      return h->get_next() ;
    }


    /**
     * \fn Halfedge* get_mate( Halfedge* h ) const
     *
     * \brief  Returns the  mate  of  a given  half-edge.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the mate half-edge of a given half-edge.
     *
     */
    inline Halfedge* get_mate( Halfedge* h ) const 
    {
#ifdef DEBUGMODE
      assert( h != 0 ) ;
      //assert( h->get_face()->get_surface() == _mesh ) ;
#endif
 
      return h->get_mate() ;
    }


    /**
     * \fn inline Halfedge* get_halfedge( Face* f ) const
     *
     * \brief Returns  the first half-edge of the  cycle of half-edges
     * of a given face of the mesh.
     *
     * \param f Pointer to a face of the mesh.
     *
     * \return A pointer to the first half-edge of a given face.
     */
    inline Halfedge* get_halfedge( Face* f ) const
    {
#ifdef DEBUGMODE
      assert( f != 0 ) ;
      //assert( f->get_surface() == _mesh ) ;
#endif

      return f->get_halfedge() ;
    }


    /**
     * \fn inline Halfedge* get_halfedge( Vertex* v ) const
     *
     * \brief Returns one half-edge with origin at a given vertex.  It
     * is  assumed  that  this  method  will always  return  the  same
     * half-edge  (as  many  half-edges  may  share  the  same  origin
     * vertex).
     *
     * \param v Pointer to a vertex of the mesh.
     *
     * \return  A pointer  to one  half-edge  with origin  at a  given
     * vertex.
     */
    inline Halfedge* get_halfedge( Vertex* v ) const
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      return v->get_halfedge() ;
    }


    /**
     * \fn inline VertexIterator vertices_begin() const
     *
     * \brief Returns a vertex iterator set to the initial vertex of a
     * vertex sequence of the mesh.
     *
     * \return A vertex iterator set to the initial vertex of a vertex
     * sequence of the mesh.
     */
    inline VertexIterator vertices_begin() const 
    {
      return _mesh->vertices_begin() ;
    }


    /**
     * \fn inline bool is_done( const VertexIterator& iterator ) const
     *
     * \brief Returns  a logic value  true if a given  vertex iterator
     * has  reached  the  end  of  a  vertex  sequence  of  the  mesh;
     * otherwise, it returns the logic value false.
     *
     * \param iterator A vertex iterator.
     *
     * \return  A logic  value true  if  a given  vertex iterator  has
     * reached the end of a vertex sequence of the mesh; otherwise, it
     * returns the logic value false.
     */
    inline bool is_done( const VertexIterator& iterator ) const
    {
      return ( iterator == _mesh->vertices_end() ) ;
    }


    /**
     * \fn inline void move_forward( VertexIterator& iterator ) const
     *
     * \brief Makes  the iterator point  to the vertex  succeeding its
     * current vertex in a vertex sequence of the mesh.
     *
     * \param iterator A reference to a vertex iterator.
     */
    inline void move_forward( VertexIterator& iterator ) const
    {
      ++iterator ;
    }


    /**
     * \fn inline Vertex* get_vertex( const VertexIterator& iterator ) const
     *
     * \brief Returns  the current vertex  of a given  vertex iterator
     * for a vertex sequence of the mesh.
     *
     * \param iterator A reference to a vertex iterator.
     *
     * \return  A pointer  to the  current  vertex of  a given  vertex
     * iterator for a vertex sequence of the mesh.
     */
    inline Vertex* get_vertex( const VertexIterator& iterator ) const
    {
      return *iterator ;
    }


    /**
     * \fn virtual EdgeIterator edges_begin() const
     *
     * \brief Returns an  edge iterator set to the  initial edge of an
     * edge sequence of the mesh.
     *
     * \return An  edge iterator  set to the  initial edge of  an edge
     * sequence of the mesh.
     *
     */
    inline EdgeIterator edges_begin() const 
    {
      return _mesh->edges_begin() ;
    }


    /**
     * \fn inline bool is_done( const EdgeIterator& iterator ) const
     *
     * \brief Returns A logic value  true if a given edge iterator has
     * reached the end of an  edge sequence of the mesh; otherwise, it
     * returns the logic value false.
     *
     * \param iterator An edge iterator.
     *
     * \return A logic value true if a given edge iterator has reached
     * the end of an edge  sequence of the mesh; otherwise, it returns
     * the logic value false.
     */
    inline bool is_done( const EdgeIterator& iterator ) const
    {
      return ( iterator == _mesh->edges_end() ) ;
    }


    /**
     * \fn inline void move_forward( EdgeIterator& iterator ) const
     *
     * \brief  Makes the  iterator point  to the  edge  succeeding its
     * current edge in an edge sequence of the mesh.
     *
     * \param iterator A reference to an edge iterator.
     */
    inline void move_forward( EdgeIterator& iterator ) const 
    {
      ++iterator ;
    }


    /**
     * \fn inline Edge* get_edge( const EdgeIterator& iterator ) const
     *
     * \brief Returns the current edge of a given edge iterator for an
     * edge sequence of the mesh.
     *
     * \param iterator A reference to an edge iterator.
     *
     * \return A pointer to the  current edge of a given edge iterator
     * for an edge sequence of the mesh.
     */
    inline Edge* get_edge( const EdgeIterator& iterator ) const
    {  
      return *iterator ;
    }


    /**
     * \fn inline FaceIterator faces_begin() const
     *
     * \brief Returns  a face  iterator set to  the initial face  of a
     * face sequence of the mesh.
     *
     * \return  A face  iterator set  to the  initial face  of  a face
     * sequence of the mesh.
     */
    inline FaceIterator faces_begin() const 
    {
      return _mesh->faces_begin() ;
    }


    /**
     * \fn virtual bool is_done( const FaceIterator& iterator ) const
     *
     * \brief Returns A logic value  true if a given face iterator has
     * reached the end  of a face sequence of  the mesh; otherwise, it
     * returns the logic value false.
     *
     * \param iterator A face iterator.
     *
     * \return A logic value true if a given face iterator has reached
     * the end of  a face sequence of the  mesh; otherwise, it returns
     * the logic value false.
     */
    inline bool is_done( const FaceIterator& iterator ) const 
    {
      return ( iterator == _mesh->faces_end() ) ;
    }


    /**
     * \fn inline void move_forward( FaceIterator& iterator ) const
     *
     * \brief  Makes the  iterator point  to the  face  succeeding its
     * current face in a face sequence of the mesh.
     *
     * \param iterator A reference to a face iterator.
     */
    inline void move_forward( FaceIterator& iterator ) const 
    {
      ++iterator ;
    }


    /**
     * \fn inline Face* get_face( const FaceIterator& iterator ) const
     *
     * \brief Returns the current face  of a given face iterator for a
     * face sequence of the mesh.
     *
     * \param iterator A reference to a face iterator.
     *
     * \return A pointer to the  current face of a given face iterator
     * for a face sequence of the mesh.
     */
    inline Face* get_face( const FaceIterator& iterator ) const 
    {
      return *iterator ;
    }


    /**
     * \fn inline void get_halfedges( Edge* e, Halfedge*& h1, Halfedge*& h2 ) const
     *
     * \brief Returns both halfedges that define a given edge.
     *
     * \param e A reference to an edge iterator.
     * \param h1 A pointer to receive the first halfedge.
     * \param h2 A pointer to receive the second halfedge.
     *
     * \return Two pointers to the halfedges of a given edge.
     */
    inline void get_halfedges( 
			      Edge* e , 
			      Halfedge*& h1 ,
			      Halfedge*& h2 
			     ) 
      const 
    {
#ifdef DEBUGMODE
      assert( e != 0 ) ;
      //assert( e->get_first_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      h1 = e->get_first_halfedge() ;
      h2 = e->get_second_halfedge() ;

      return ;
    }


    /**
     * \fn inline void get_coords( Vertex* v, double& x , double& y , double& z ) const
     *
     * \brief Returns the three coordinates of a given vertex.
     *
     * \param v A reference to the given vertex.
     * \param x Returns X coordinate.
     * \param y Returns Y coordinate.
     * \param z Returns Z coordinate.
     */
    inline void get_coords( 
			   Vertex* v ,
			   double& x ,
			   double& y ,
			   double& z 
			  ) 
      const 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      x = v->x() ;
      y = v->y() ;
      z = v->z() ;

      return ;
    }


    /**
     * \fn inline void set_coords( Vertex* v, double x, double y, double z )
     *
     * \brief Define the three coordinates of a given vertex.
     *
     * \param v A reference to the given vertex;
     * \param x X coordinate;
     * \param y Y coordinate;
     * \param z Z coordinate;
     *
     */
    inline void set_coords( 
			   Vertex* v ,
			   double x ,
			   double y ,
			   double z 
			  ) 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      v->set_x_coord( x ) ;
      v->set_y_coord( y ) ;
      v->set_z_coord( z ) ;

      return ;
    }


    /**
     * \fn inline void flip_edge( Edge* e )
     *
     * \brief Flips a given edge of the base mesh.
     *
     * \param e Pointer to the mesh edge to be flipped.
     */
    inline void flip_edge( Edge* e )
    {
#ifdef DEBUGMODE
      assert( e != 0 ) ;
      //assert( e->get_first_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      assert( false );
      //_mesh->flip_edge( e ) ;

    }


    /**
     * \fn inline void remove_vertex( Vertex* v )
     *
     * \brief Removes a vertex and all of its incident faces. 
     *
     * \param v A pointer to the vertex to be removed.
     */
    inline void remove_vertex( Vertex* v )
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

        assert ( false );
     //_mesh->remove_vertex( v ) ;

    }


    /**
     * \fn inline bool is_alive( Vertex* v ) const
     *
     * \brief Checks if the vertex is in the alive set.
     *
     * \param v A pointer to a given vertex.
     *
     * \return  The logic value  true if  the given  vertex is  in the
     * alive set, and the logic value false otherwise.
     */
    inline bool is_alive( Vertex* v ) const 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      return v->get_attributes().alive() ;
    }


   /**
    * \fn inline void set_alive( Vertex* v, bool a )
    *
    * \brief  Assigns  a  given   logic  value  to  the  vertex  flag
    * (supposedly a vertex attribute) that indicates if the vertex is
    * in the alive set.
    *
    * \param v A pointer to a given vertex.
    * \param a A logic value (i.e., true or false).
    */
    inline void set_alive( Vertex* v , bool a ) 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      v->get_attributes().set_alive( a ) ;
    }


    /**
     * \fn inline bool is_close( Vertex* v ) const
     *
     * \brief Checks if the vertex is in the alive close.
     *
     * \param v A pointer to a given vertex.
     *
     * \return  The logic value  true is  the given  vertex is  in the
     * close set, and the logic value false otherwise.
     */
    inline bool is_close( Vertex* v ) const 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      return v->get_attributes().close() ;
    }


    /**
     * \fn inline void set_close( Vertex* v, bool a )
     *
     * \brief  Assigns  a  given   logic  value  to  the  vertex  flag
     * (supposedly a vertex attribute) that indicates if the vertex is
     * in the close set.
     *
     * \param v A pointer to a given vertex.
     * \param a A logic value (i.e., true or false).
     */
    inline void set_close( 
			  Vertex* v ,
			  bool a 
			 ) 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      v->get_attributes().set_close( a ) ;

      return ;
    }


    /**
     * \fn inline double get_dist( Vertex* v ) const
     *
     * \brief Returns  the distance  function value associated  with a
     * given vertex.
     *
     * \param v A pointer to a given vertex.
     *
     * \return The distance funtion value at the given vertex.
     */
    inline double get_dist( Vertex* v ) const 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      return v->get_attributes().dist() ;
    }


    /**
     * \fn inline void set_dist( Vertex* v, double d )
     *
     * \brief Assigns a given  distance value to the distance function
     * attribute of a given vertex.
     *
     * \param v A pointer to a given vertex.
     * \param d A distance value (i.e., a real number).
     */
    inline void set_dist( 
			 Vertex* v ,
			 double d 
			 ) 
    {
#ifdef DEBUGMODE
      assert( v != 0 ) ;
      //assert( v->get_halfedge()->get_face()->get_surface() == _mesh ) ;
#endif

      v->get_attributes().set_dist( d ) ;
    }


    /**
     * \fn inline void activate( Face* f )
     *
     * \brief Make a mesh face active.
     *
     * \param f A pointer to a face.
     */
    inline void activate( Face* f )
    {
#ifdef DEBUGMODE
      assert( f != 0 ) ;
      //assert( f->get_surface() == _mesh ) ;
#endif

      f->get_attributes().activate() ;
    }


    /**
     * \fn inline void deactivate( Face* f )
     *
     * \brief Make a mesh face active.
     *
     * \param f A pointer to a face.
     */
    inline void deactivate( Face* f )
    {
#ifdef DEBUGMODE
      assert( f != 0 ) ;
      //assert( f->get_surface() == _mesh ) ;
#endif

      f->get_attributes().deactivate() ;
    }
      

    /**
     * \fn inline bool is_active( Face* f ) const
     *
     * \brief Decides whether a given mesh face is active.
     *
     * \param f A pointer to a face.
     *
     * \return The  logic value  true if the  face is active,  and the
     * logic value false otherwise.
     */
    inline bool is_active( Face* f ) const
    {
#ifdef DEBUGMODE
      assert( f != 0 ) ;
      //assert( f->get_surface() == _mesh ) ;
#endif

      return f->get_attributes().is_active() ;
    }

  private:
    // ---------------------------------------------------------------
    //
    // Private data members.
    //
    // ---------------------------------------------------------------

    TMesh* _mesh ;   ///< Pointer to the actual mesh.

  } ;

}

/** @} */ //end of group class.

#endif	/* DM_FROM_DCEL_HPP */

