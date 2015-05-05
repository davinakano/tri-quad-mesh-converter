/** 
 * \file geodesic_point.hpp
 *  
 * \brief   Definition   and    implementation   of   abstract   class
 * Geodesic_point, which  represents a  vertex of a  discrete geodesic
 * (an open, simple polygonal line). The vertex can be a vertex of the
 * mesh over  which the geodesic is  defined, or it can  belong to the
 * interior of a mesh edge or face.
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

#ifndef GEODESIC_POINT_HPP
#define GEODESIC_POINT_HPP

#include "mesh_interface.hpp"      // common::Mesh_interface

#include <vector>                  // std::vector


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
   * \class Geodesic_point
   *
   * \brief This  class represents represents  a vertex of  a discrete
   * geodesic (an  open, simple polygonal  line). The vertex can  be a
   * vertex of the mesh over which  the geodesic is defined, or it can
   * belong to the interior of a mesh edge or face.
   */
  template< typename Mesh >
  class Geodesic_point {
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

 
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Geodesic_point()
     *
     * \brief Creates an  instance of this class.
     */
    Geodesic_point() 
    {}


    /**
     * \fn virtual ~Geodesic_point()
     *
     * \brief Destroys an instance of this class.
     */    
    virtual ~Geodesic_point()
    {}


    /**
     * \fn virtual bool on_vertex() const = 0
     *
     * \brief Determines  whether this geodesic vertex is  a vertex of
     * the mesh over which the geodesic is defined.
     *
     * \return  The logic  value true  if  this geodesic  vertex is  a
     * vertex of the mesh over  which the geodesic is defined, and the
     * logic value false otherwise.
     */    
    virtual bool on_vertex() const = 0 ;


    /**
     * \fn virtual bool on_edge() const = 0 ;
     *
     * \brief Determines  whether this geodesic vertex  belongs to the
     * interior  of an edge  of the  mesh over  which the  geodesic is
     * defined.
     *
     * \return The logic value true  if this geodesic vertex is in the
     * interior  of an edge  of the  mesh over  which the  geodesic is
     * defined, and the logic value false otherwise.
     */  
    virtual bool on_edge() const = 0 ;


    /**
     * \fn virtual bool on_face() const = 0 ;
     *
     * \brief Determines  whether this geodesic vertex  belongs to the
     * interior  of a  face of  the mesh  over which  the  geodesic is
     * defined.
     *
     * \return The logic value true  if this geodesic vertex is in the
     * interior  of a  face of  the mesh  over which  the  geodesic is
     * defined, and the logic value false otherwise.
     */  
    virtual bool on_face() const = 0 ;


    /**
     * \fn virtual void get_coords( MI* m, double& x , double& y , double& z ) const = 0
     *
     * \brief Gets the three coordinates of this geodesic point.
     *
     * \param  m  Pointer to  the  mesh  over  which the  geodesic  is
     * defined.
     * \param x A reference to the X coordinate.
     * \param y A reference to the Y coordinate.
     * \param z A reference to the Z coordinate.
     */
    virtual void get_coords( MI* m, double& x , double& y , 
      double& z ) const = 0 ;


    /**
     * \fn virtual bool in_the_same_mesh_element( MI* m , Geodesic_point< Mesh >* v , std::list< Geodesic_point< Mesh >* >& lv ) const = 0
     *
     * \brief Tests whether this  geodesic vertex and a given geodesic
     * vertex are  the same, or lie  in the interior of  the same mesh
     * edge, or lie in the  interior of distinct mesh edges that share
     * a face,  or lie in  the interior of  the same face. If  so, the
     * method  inserts  this geodesic  vertex  in  the  given list  of
     * vertices (and  possibly the given vertex too),  and returns the
     * logic value true. Otherwise, the method returns the logic value
     * false.
     *
     * \param  m  Pointer  to   the  mesh  that  implements  the  mesh
     * interface.
     * \param v Pointer to a geodesic vertex.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return  The logic  value true  if this  geodesic vertex  and a
     * given geodesic vertex  are the same, or lie  in the interior of
     * the same  mesh edge,  or lie in  the interior of  distinct mesh
     * edges that  share a face,  or lie in  the interior of  the same
     * face. Otherwise, it returns the logic value false.
     */
    virtual bool in_the_same_mesh_element( MI* m , 
      Geodesic_point< Mesh >* v , 
      std::list< Geodesic_point< Mesh >* >& lv ) const = 0 ;


    /**
     * \fn virtual void get_adj_vertices( MI* m , std::vector< Vertex* >& a ) const = 0
     *
     * \brief  Finds all  mesh vertices  adjacent to  the  mesh vertex
     * correponding to  this geodesic vertex. If  this geodesic vertex
     * lies  in the  interior of  an edge  of the  mesh, it  finds all
     * vertices  in the edge  star. Finally,  if this  geodesic vertex
     * lies in the interior of a  mesh face, finds the vertices of the
     * face.
     *
     * \param m Pointer to the mesh over which the geodesic is defined.
     * \param a A reference to a list of mesh vertices.
     */
    virtual void get_adj_vertices( MI* m , std::vector< Vertex* >& a ) 
      const = 0 ;

  } ;

}

/** @} */ //end of group class.

#endif	/* GEODESIC_POINT_HPP */

