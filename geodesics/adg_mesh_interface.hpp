/** 
 * \file adg_mesh_interface.hpp
 *  
 * \brief  Definition   and  implementation  of   the  template  class
 * Adg_mesh_interface,  which is  an abstract  interface  for triangle
 * surface  meshes   where  approximate  discrete   geodesics  can  be
 * computed.
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

#ifndef ADG_MESH_INTERFACE_HPP
#define	ADG_MESH_INTERFACE_HPP

#include "mesh_interface.hpp"  // common::Mesh_interface


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


  using common::Mesh_interface ;


  /**
   * \class Adg_mesh_interface
   *
   * \brief  This  class  represents   an  abstract  interface  for  a
   * half-edge  based  generic triangle  mesh  over which  approximate
   * discrete geodesics can be computed.
   */
  template< typename Mesh >
  class Adg_mesh_interface : public Mesh_interface< Mesh > {
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
     * \fn virtual bool is_alive( Vertex* v ) const = 0
     *
     * \brief Checks if the vertex is in the alive set.
     *
     * \param v A pointer to a vertex.
     *
     * \return  The logic value  true if  the given  vertex is  in the
     * alive set, and the logic value false otherwise.
     */
    virtual bool is_alive( Vertex* v ) const = 0 ;


    /**
     * \fn virtual void set_alive( Vertex* v, bool a ) = 0
     *
     * \brief  Assigns  a  given   logic  value  to  the  vertex  flag
     * (supposedly a vertex attribute) that indicates if the vertex is
     * in the alive set.
     *
     * \param v A pointer to a vertex.
     * \param a A logic value (i.e., true or false).
     */
    virtual void set_alive( Vertex* v , bool a ) = 0 ;


    /**
     * \fn virtual bool is_close( Vertex* v ) const = 0
     *
     * \brief Checks if the vertex is in the alive close.
     *
     * \param v A pointer to a vertex.
     *
     * \return  The logic value  true is  the given  vertex is  in the
     * close set, and the logic value false otherwise.
     */
    virtual bool is_close( Vertex* v ) const = 0 ;


    /**
     * \fn virtual void set_close( Vertex* v , bool a ) = 0
     *
     * \brief  Assigns  a  given   logic  value  to  the  vertex  flag
     * (supposedly a vertex attribute) that indicates if the vertex is
     * in the close set.
     *
     * \param v A pointer to a vertex.
     * \param a A logic value (i.e., true or false).
     */
    virtual void set_close( Vertex* v , bool a ) = 0 ;


    /**
     * \fn virtual double get_dist( Vertex* v ) const = 0
     *
     * \brief Returns  the distance  function value associated  with a
     * given vertex.
     *
     * \param v A pointer to a vertex.
     *
     * \return The distance funtion value at the given vertex.
     */
    virtual double get_dist( Vertex* v ) const = 0 ;


    /**
     * \fn virtual void set_dist( Vertex* v, double d ) = 0
     *
     * \brief Assigns a given  distance value to the distance function
     * attribute of a given vertex.
     *
     * \param v A pointer to a given vertex.
     * \param d A distance value (i.e., a real number).
     */
    virtual void set_dist( Vertex* v , double d ) = 0 ;


    /**
     * \fn virtual void activate( Face* f ) = 0
     *
     * \brief Make a mesh face active.
     *
     * \param f A pointer to a face.
     */
    virtual void activate( Face* f ) = 0 ;


    /**
     * \fn virtual void deactivate( Face* f ) = 0
     *
     * \brief Make a mesh face active.
     *
     * \param f A pointer to a face.
     */
    virtual void deactivate( Face* f ) = 0 ;


    /**
     * \fn virtual bool is_active( Face* f ) const = 0
     *
     * \brief Decides whether a given mesh face is active.
     *
     * \param f A pointer to a face.
     *
     * \return The  logic value  true if the  face is active,  and the
     * logic value false otherwise.
     */
    virtual bool is_active( Face* f ) const = 0 ;


    /**
     * \fn Adg_mesh_interface()
     *
     * \brief Creates an  instance of this class. This  is an abstract
     * class, so this constructor  will never be called to instantiate
     * an object  of this class,  but by the constructor  of inherited
     * classes.
     */
    Adg_mesh_interface() 
    {}


    /**
     * \fn ~Adg_mesh_interface()
     *
     * \brief Destroys an instance of this class.
     */
    virtual ~Adg_mesh_interface() 
    {}

  } ;

}

/** @} *///end of group class.

#endif	/* ADG_MESH_INTERFACE_HPP */

