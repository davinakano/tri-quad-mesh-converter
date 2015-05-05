/** 
 * \file mesh_interface.hpp
 *  
 * \brief  Definition   and  implementation  of   the  template  class
 * Mesh_interface, which  is an abstract interface for  a generic data
 * structure for representing triangle surface meshes in the Euclidean
 * space.
 *
 * \author
 * Marcelo Ferreira Siqueira
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

#ifndef MESH_INTERFACE_HPP
#define	MESH_INTERFACE_HPP


/** 
 * \defgroup COMMONNameSpace Namespace common.
 * @{
 */

/**
 * \namespace common
 *
 * \brief   The   namespace  common   contains   the  definition   and
 * implementation  of several  general classes,  which  contains basic
 * algorithms and/or  basic data structures  that can be used  by more
 * specialized classes.
 */

namespace common {

  /**
   * \class Mesh_interface
   *
   * \brief This  class represents the  template class Mesh_interface,
   * which  is an  abstract interface  for a  half-edge  based generic
   * triangle mesh.
   */
  template< typename Mesh >
  class Mesh_interface {
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
    typedef typename Mesh::Vertex Vertex ;


    /**
     * \typedef Halfedge
     *
     * \brief Definition of a type name for the mesh half-edges.
     */
    typedef typename Mesh::Halfedge Halfedge ;


    /**
     * \typedef Edge
     *
     * \brief Definition of a type name for the mesh edges.
     */
    typedef typename Mesh::Edge Edge ;


    /**
     * \typedef Face
     *
     * \brief Definition of a type name for the mesh faces.
     */
    typedef typename Mesh::Face Face ;


    /**
     * \typedef VertexIterator
     *
     * \brief Definition of a type name for the vertex iterators.
     */
    typedef typename Mesh::VertexIterator VertexIterator ;


    /**
     * \typedef EdgeIterator
     *
     * \brief Definition of a type name for the edge iterators.
     */
    typedef typename Mesh::EdgeIterator EdgeIterator ;


    /**
     * \typedef FaceIterator.
     *
     * \brief Definition of a type name for the face iterators.
     */
    typedef typename Mesh::FaceIterator FaceIterator ;


    // ---------------------------------------------------------------
    //
    // Public pure virtual methods 
    //
    // ---------------------------------------------------------------


    /**
     * \fn virtual Vertex* get_org( Halfedge* h ) const
     *
     * \brief Returns the  origin vertex of a given  half-edge of the
     * mesh.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the origin vertex of a given half-edge.
     */
    virtual Vertex* get_org( Halfedge* h ) const = 0 ;


    /**
     * \fn virtual Edge* get_edge( Halfedge* h ) const
     *
     * \brief Returns the  edge a given half-edge of  the mesh belongs
     * to.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the edge a half-edge belongs to.
     */
    virtual Edge* get_edge( Halfedge* h ) const = 0 ;


    /**
     * \fn virtual Face* get_face( Halfedge* h ) const
     *
     * \brief Returns the  face a given half-edge of  the mesh belongs
     * to.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the face a half-edge belongs to.
     */
    virtual Face* get_face( Halfedge* h ) const = 0 ;


    /**
     * \fn virtual Halfedge* get_prev( Halfedge* h ) const
     *
     * \brief Returns the half-edge that precedes a given half-edge of
     * the mesh  in the  face cycle of  half-edges that  contains both
     * half-edges.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return  A  pointer to  the  half-edge  that  precedes a  given
     * half-edge in their common face half-edge cycle.
     *
     */
    virtual Halfedge* get_prev( Halfedge* h ) const = 0 ;


    /**
     * \fn virtual Halfedge* get_next( Halfedge* h ) const
     *
     * \brief Returns the half-edge that succeeds a given half-edge of
     * the mesh  in the  face cycle of  half-edges that  contains both
     * half-edges.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return  A  pointer to  the  half-edge  that  succeeds a  given
     * half-edge in their common face half-edge cycle.
     */
    virtual Halfedge* get_next( Halfedge* h ) const = 0 ;


    /**
     * \fn virtual Halfedge* get_mate( Halfedge* h ) const
     *
     * \brief  Returns the  mate  of  a given  half-edge.
     *
     * \param h Pointer to a half-edge of the mesh.
     *
     * \return A pointer to the mate half-edge of a half-edge.
     */
    virtual Halfedge* get_mate( Halfedge* h ) const = 0 ;


    /**
     * \fn virtual Halfedge* get_halfedge( Face* face ) const
     *
     * \brief Returns  the first half-edge of the  cycle of half-edges
     * of a given face of the mesh.
     *
     * \param face Pointer to a face of the mesh.
     *
     * \return A pointer to the first half-edge of a face.
     */
    virtual Halfedge* get_halfedge( Face* face ) const = 0 ;


    /**
     * \fn virtual Halfedge* get_halfedge( Vertex* vertex ) const
     *
     * \brief Returns one half-edge with origin at a given vertex.  It
     * is  assumed  that  this  method  will always  return  the  same
     * half-edge  (as  many  half-edges  may  share  the  same  origin
     * vertex).
     *
     * \param vertex Pointer to a vertex of the mesh.
     *
     * \return  A pointer  to one  half-edge with  origin at  a given
     * vertex.
     */
    virtual Halfedge* get_halfedge( Vertex* vertex ) const = 0 ;


    /**
     * \fn virtual VertexIterator vertices_begin() const
     *
     * \brief Returns a vertex iterator set to the initial vertex of a
     * vertex sequence of the mesh.
     *
     * \return A vertex iterator set to the initial vertex of a vertex
     * sequence of the mesh.
     */
    virtual VertexIterator vertices_begin( ) const = 0 ;


    /**
     * \fn virtual bool is_done( const VertexIterator& iterator ) const
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
    virtual bool is_done( const VertexIterator& iterator ) const = 0 ;


    /**
     * \fn virtual void move_forward( VertexIterator& iterator ) const
     *
     * \brief Makes  the iterator point  to the vertex  succeeding its
     * current vertex in a vertex sequence of the mesh.
     *
     * \param iterator A reference to a vertex iterator.
     */
    virtual void move_forward( VertexIterator& iterator ) const = 0 ;


    /**
     * \fn virtual Vertex* get_vertex( const VertexIterator& iterator ) const
     *
     * \brief Returns  the current vertex  of a given  vertex iterator
     * for a vertex sequence of the mesh.
     *
     * \param iterator A reference to a vertex iterator.
     *
     * \return A  pointer to the  current vertex of a  vertex iterator
     * for a vertex sequence of the mesh.
     */
    virtual Vertex* get_vertex( const VertexIterator& iterator ) const = 0 ;


    /**
     * \fn virtual EdgeIterator edges_begin() const
     *
     * \brief Returns an  edge iterator set to the  initial edge of an
     * edge sequence of the mesh.
     *
     * \return An  edge iterator  set to the  initial edge of  an edge
     * sequence of the mesh.
     */
    virtual EdgeIterator edges_begin( ) const = 0 ;


    /**
     * \fn virtual bool is_done( const EdgeIterator& iterator ) const
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
    virtual bool is_done( const EdgeIterator& iterator ) const = 0 ;


    /**
     * \fn virtual void move_forward( EdgeIterator& iterator ) const
     *
     * \brief  Makes the  iterator point  to the  edge  succeeding its
     * current edge in an edge sequence of the mesh.
     *
     * \param iterator A reference to an edge iterator.
     */
    virtual void move_forward( EdgeIterator& iterator ) const = 0 ;


    /**
     * \fn virtual Edge* get_edge( const EdgeIterator& iterator ) const
     *
     * \brief Returns the current edge of a given edge iterator for an
     * edge sequence of the mesh.
     *
     * \param iterator A reference to an edge iterator.
     *
     * \return A pointer  to the current edge of  an edge iterator for
     * an edge sequence of the mesh.
     */
    virtual Edge* get_edge( const EdgeIterator& iterator ) const = 0 ;


    /**
     * \fn virtual FaceIterator faces_begin() const
     *
     * \brief Returns  a face  iterator set to  the initial face  of a
     * face sequence of the mesh.
     *
     * \return  A face  iterator set  to the  initial face  of  a face
     * sequence of the mesh.
     */
    virtual FaceIterator faces_begin( ) const = 0 ;


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
    virtual bool is_done( const FaceIterator& iterator ) const = 0 ;


    /**
     * \fn virtual void move_forward( FaceIterator& iterator ) const
     *
     * \brief  Makes the  iterator point  to the  face  succeeding its
     * current face in a face sequence of the mesh.
     *
     * \param iterator A reference to a face iterator.
     */
    virtual void move_forward( FaceIterator& iterator ) const = 0 ;


    /**
     * \fn virtual Face* get_face( const FaceIterator& iterator ) const
     *
     * \brief Returns the current face  of a given face iterator for a
     * face sequence of the mesh.
     *
     * \param iterator A reference to a face iterator.
     *
     * \return A pointer to the current  face of a face iterator for a
     * face sequence of the mesh.
     *
     */
    virtual Face* get_face( const FaceIterator& iterator ) const = 0 ;


    /**
     * \fn virtual void get_halfedges( Edge* e , Halfedge*& h1 , Halfedge*& h2 ) const = 0
     *
     * \brief Returns both halfedges that define a given edge.
     *
     * \param e Pointer to a mesh edge.
     * \param h1 A pointer to receive the first halfedge
     * \param h2 A pointer to receive the second halfedge
     *
     * \return Two pointers to the halfedges of a given edge.
     */
    virtual void get_halfedges( Edge* e , Halfedge*& h1 , Halfedge*& h2 ) const = 0 ;

    /**
     * \fn virtual void get_coords( Vertex* v, double& x , double& y , double& z ) const = 0
     *
     * \brief Returns the three coordinates of a given vertex.
     *
     * \param v Pointer to a vertex.
     * \param x A reference to the X coordinate.
     * \param y A reference to the Y coordinate.
     * \param z A reference to the Z coordinate.
     */
    virtual void get_coords( Vertex* v , double& x , double& y , double& z ) const = 0 ;


    /**
     * \fn virtual void set_coords( Vertex* v , double x , double y , double z ) = 0
     *
     * \brief  Assigns values  to  the three  coordinates  of a  given
     * vertex.
     *
     * \param v Pointer to a vertex.
     * \param x The X coordinate.
     * \param y The Y coordinate.
     * \param z The Z coordinate.
     */
    virtual void set_coords( Vertex* v , double x , double y , double z ) = 0 ;


    /**
     * \fn virtual void flip_edge( Edge* e ) = 0
     *
     * \brief Flips a given edge of the mesh.
     *
     * \param e Pointer to the mesh edge to be flipped.
     */
    virtual void flip_edge( Edge* e ) = 0 ;


    /**
     * \fn virtual void remove_vertex( Vertex* v ) = 0
     *
     * \brief Removes a vertex and all of its incident faces. 
     *
     * \param v A pointer to the vertex to be removed.
     */
    virtual void remove_vertex( Vertex* v ) = 0 ;


    /**
     * \fn Mesh_interface()
     *
     * \brief Creates an  instance of this class. This  is an abstract
     * class, so this constructor  will never be called to instantiate
     * an object  of this class,  but by the constructor  of inherited
     * classes.
     */
    Mesh_interface() 
    {}


    /**
     * \fn ~Mesh_interface()
     *
     * \brief Destroys an instance of this class.
     */
    virtual ~Mesh_interface() 
    {}

  } ;

}

/** @} *///end of group class.

#endif	/* MESH_INTERFACE_HPP */

