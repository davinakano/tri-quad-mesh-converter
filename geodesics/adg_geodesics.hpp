/** 
 * \file adg_geodesics.hpp
 *  
 * \brief  Definition  and implementation  of  class Geodesics,  which
 * computes  approximate  discrete  geodesics  defined  over  triangle
 * surface meshes.
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

#ifndef ADG_GEODESICS_HPP
#define	ADG_GEODESICS_HPP

#include <algorithm>                     // std::min, std::max, find
#include <list>                          // std::list
#include <vector>                        // std::vector
#include <cmath>                         // sqrt, acos, sin, fabs

#include "adg_mesh_interface.hpp"        // Adg_mesh_interface

#include "geodesic_point.hpp"            // Geodesic_point
#include "geodesic_point_vertex.hpp"     // Geodesic_point_vertex
#include "geodesic_point_edge.hpp"       // Geodesic_point_edge
#include "geodesic_point_face.hpp"       // Geodesic_point_face

#include "adg_curve_point.hpp"           // Curve_point
#include "priority_queue.hpp"            // dynamic_priority_queue
#include "geometric.hpp"                 // Geometric

#include <iostream>

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

  using common::Geometric ;
  using common::dynamic_priority_queue ;

  /**
   * \class Approx_geodesics
   *
   * \brief  This  class represents  a  finder  of discrete  geodesics
   * defined on a  triangle surface mesh. The class  is independent of
   * the surface mesh data structure, which is a template parameter of
   * the class.
   */
  template< typename Mesh >
  class Approx_geodesics {
  private:
    // ---------------------------------------------------------------
    //
    // Constant definitions
    //
    // ---------------------------------------------------------------

    const double _MYPI ; ///< The constant PI.
    const double _TOL_T_EDGE ; ///< Minimum allowed point-in-edge movement threshold.
    const double _TOL_E_NEAR_V ; ///< Threshold for snapping point-in-edge to a mesh vertex.
    const double _SHORTENING_RATIO ;  ///< Threshold for considering two curves the same.

  public:
    // ---------------------------------------------------------------
    //
    // Type name definitions for the generic mesh components.
    //
    // ---------------------------------------------------------------

    /**
     * \typedef MI
     *
     * \brief Definition of a type name for the mesh interface.
     */
    typedef Adg_mesh_interface< Mesh > MI ;

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
     * \typedef VertexIterator
     *
     * \brief Definition of a type name for the vertex iterators.
     */
    typedef typename MI::VertexIterator VertexIterator ;


    /**
     * \typedef EdgeIterator
     *
     * \brief Definition of a type name for the edge iterators.
     */
    typedef typename MI::EdgeIterator EdgeIterator ;


    /**
     * \typedef FaceIterator.
     *
     * \brief Definition of a type name for the face iterators.
     */
    typedef typename MI::FaceIterator FaceIterator ;


    /**
     * \typedef GP
     *
     * \brief Definition of a type name for a Geodesic_point, which is
     * a geodesic vertex that can be  either a mesh vertex, a point in
     * the interior  of a mesh edge, or  a point in the  interior of a
     * mesh face.
     */
    typedef Geodesic_point< Mesh > GP ;


    /**
     * \typedef GPV
     *
     * \brief Definition  of a type name  for a Geodesic_point_vertex,
     * which is a geodesic vertex. This vertex coincides with a vertex
     * of the mesh over which the geodesic it belongs to is defined.
     */
    typedef Geodesic_point_vertex< Mesh > GPV ;


    /**
     * \typedef GPE
     *
     * \brief  Definition of  a type  name for  a Geodesic_point_edge,
     * which is a geodesic vertex that lies in the interior of an edge
     * of the mesh over which the geodesic it belongs to is defined.
     */
    typedef Geodesic_point_edge< Mesh > GPE ;


    /**
     * \typedef GPF
     *
     * \brief  Definition of  a type  name for  a Geodesic_point_face,
     * which is a geodesic vertex that  lies in the interior of a face
     * of the mesh over which the geodesic it belongs to is defined.
     */
    typedef Geodesic_point_face< Mesh > GPF ;


    /**
     * \typedef CP
     *
     * \brief  Definition of  a  type name  for  a Curve_point,  which
     * represents  a  geodesic  vertex  during the  execution  of  the
     * algorithm  for computing  a discrete  geodesics  connecting two
     * mesh points.
     */
    typedef Curve_point< Mesh > CP ;


    /**
     * \typedef INTERSECTION_TYPE
     *
     * \brief Definition of a type name for all possible return values
     * of  the algorithm that  computes line  segment to  line segment
     * intersection.
     */
    typedef typename Geometric::INTERSECTION_TYPE INTERSECTION_TYPE ;


    /**
     * \typedef SIDE
     *
     * \brief  Definition of  a type  name for  the three  cases  of a
     * geodesic passing through the star of a mesh vertex: it can pass
     * either on  the right side  of the star vertex  (RIGHT_SIDE), or
     * the left side (LEFT_SIDE), or through the vertex (NO_SIDE).
     */
    typedef typename CP::SIDE SIDE ;

    /**
     * \typedef Path
     *
     * \brief  Definition  of a  type  name  for  the polygonal  curve
     * defining a discrete geodesic.
     */
    typedef typename std::list< CP > Path ;

    /**
     * \typedef Path_iterator
     *
     * \brief Definition of a type name for an iterator of a polygonal
     * curve defining a discrete geodesic.
     */
    typedef typename std::list< CP >::iterator Path_iterator ;

    /**
     * \typedef SortedSet
     *
     * \brief Definition of a type name for the priority queue used by
     * the Fast  Marching Method  (FMM) during the  first step  of the
     * geodesic computation.
     */
    typedef dynamic_priority_queue< Vertex* , double > SortedSet ;


    // ---------------------------------------------------------------
    //
    // Public methods.
    //
    // ---------------------------------------------------------------

    /**
     * \fn Approx_geodesics( MI* m , double error , int niter )
     *
     * \brief Creates an  instance of this class.
     *
     * \param  m  Pointer to  the  mesh  over  which the  geodesic  is
     * defined.
     * \param error Upper bound  for the maximum error associated with
     * any  geodesic  vertex.   The  error  value  associated  with  a
     * geodesic vertex  is a measure  of how straight the  geodesic is
     * around the vertex.
     * \param niter  Upper bound for  the number of iterations  of the
     * path correction algorithm.
     */
    Approx_geodesics( MI* m , double error , int niter ) :
      _MYPI( acos( -1 ) ) ,
      _TOL_T_EDGE( 1e-3 ) ,         // The value of this constant should not be changed!
      _TOL_E_NEAR_V( 1e-1 ) ,        // The value of this constant should not be changed!
      _SHORTENING_RATIO( 1.0 - 1e-10 ) ,
      _error( error ) , 
      _max_iter( niter ) , 
      _mesh( m ) ,
      _num_iter( 0 ) ,
      _h_prev( 0 ) ,
      _h_next( 0 ) ,
      _prev_x( 0 ) ,
      _prev_y( 0 ) ,
      _curr_x( 0 ) ,
      _curr_y( 0 ) ,
      _next_x( 0 ) ,
      _next_y( 0 )
    {}


    /**
     * \fn ~Approx_geodesics()
     *
     * \brief Destroys an instance of this class.
     */
    virtual ~Approx_geodesics( ) {
    }


    /**
     * \fn void compute( GP* v1 , GP* v2 , std::list< GP* >& lv )
     *
     * \brief Computes  a discrete geodesic over a  triangle mesh. The
     * initial and final  vertices of the geodesic can  be vertices of
     * the mesh,  or points in the  interior of edges or  faces of the
     * mesh.
     *
     * \param v1 The initial point of the discrete geodesic.
     * \param v2 The final point of the discrete geodesic.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return A polygonal curve  consisting of the geodesic vertices,
     * which are vertices of the mesh and/or points in the interior of
     * a mesh edge or face.
     */
    void compute( GP* v1 , GP* v2 , std::list< GP* >& lv ) ;


    /**
     * \fn void compute_first_approximation( GP* vi, GP* ve , std::list< GP* >& lv )
     *
     * \brief Computes a fast,  but rough approximation for a discrete
     * geodesic defined  over a triangle mesh.  The  initial and final
     * vertices  of the  geodesic  can  be vertices  of  the mesh  the
     * geodesic is defined  on, or points in the  interior of edges or
     * faces of the mesh.
     *
     * \param vi The initial vertex of the discrete geodesic.
     * \param ve The final vertex of the discrete geodesic.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return  A  polygonal curve  consisting  of geodesic  vertices,
     * which are vertices of the mesh (there is vertex inside edges or
     * faces).
     *
     * \sa compute_first_approximation.
     */
    void compute_first_approximation( GP* v1 , GP* v2 , 
      std::list< GP* >& lv ) ;


    /**
     * \fn void compute_restricted_first_approximation( GP* vi, GP* ve , std::list< GP* >& lv )
     *
     * \brief Computes a fast,  but rough approximation for a discrete
     * geodesic defined over a restricted set of faces of the triangle
     * mesh. The restricted faces are the active one.  Face activation
     * must be  performed before this  method is invoked.  The initial
     * and final vertices of the  geodesic can be vertices of the mesh
     * the geodesic is defined on,  or points in the interior of edges
     * or faces of the mesh.
     *
     * \param vi The initial vertex of the discrete geodesic.
     * \param ve The final vertex of the discrete geodesic.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return  A  polygonal curve  consisting  of geodesic  vertices,
     * which are vertices of the mesh (there is vertex inside edges or
     * faces).
     */
    void compute_restricted_first_approximation( GP* vi , GP* v2 , 
      std::list< GP* >& lv ) ;


    /**
     * \fn void shorten_geodesic( std::list< GP* >& lv )
     *
     * \brief Shortens a  geodesic as much as possible  using the path
     * correction algorithm.
     *
     * \param lv A list of vertices of the curve to be shortened. This
     * list  will contain the  vertices of  the shortened  curve after
     * this method is finished.
     *
     * \sa compute_first_approximation.
     */
    void shorten_geodesic( std::list< GP* >& lv ) ;


    /**
     * \fn void shorten_geodesic( const std::list< GP* >& lv_in , std::list< GP* >& lv_out )
     *
     * \brief Shortens a  geodesic as much as possible  using the path
     * correction algorithm.
     *
     * \param lv_in A list of vertices of the curve to be shortended.
     * \param lv_out A list of verttices of the shortended curve.
     *
     * \sa compute_first_approximation.
     */
    void shorten_geodesic( const std::list< GP* >& lv_in , std::list< GP* >& lv_out ) ;


    /**
     * \fn double compute_length() const
     *
     * \brief  Calculates the length of this geodesic.
     *
     * \return The length of this geodesic.
     */
    double compute_length() const ;


    /**
     * \fn unsigned get_total_number_of_iterations() const
     *
     * \brief Returns the number of iterations of the path correction
     * algorithm.
     *
     * \return  The  number  of  iterations  of  the  path  correction
     * algorithm.
     */
    unsigned get_total_number_of_iterations() const
    {
      return _num_iter ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private methods.
    //
    // ---------------------------------------------------------------

    /**
     * \fn inline void get_limits( Halfedge*& prev , Halfedge*& next ) const 
     *
     * \brief Gets the limit  half-edges associated with this geodesic
     * vertex.  These half-edges identify  the first and the last mesh
     * triangle of the sequence  of mesh triangles associated with the
     * current  geodesic vertex  and unfolded  by the  geodesic finder
     * path correction algorithm.
     *
     * \param prev A reference to the first limiting half-edge.
     * \param next A reference to the second limiting half-edge.
     */
    inline void get_limits( 
			   Halfedge*& prev , 
			   Halfedge*& next 
			  )
    const 
    {
      prev = _h_prev ;
      next = _h_next ;
    }
  

    /**
     * \fn inline void set_limits( Halfedge* prev , Halfedge* next )
     *
     * \brief Assigns  values to the limit  half-edges associated with
     * this geodesic  vertex. These half-edges identify  the first and
     * the  last  mesh triangle  of  the  sequence  of mesh  triangles
     * associated  with  this  geodesic  vertex and  unfolded  by  the
     * geodesic finder path correction algorithm.
     *
     * \param prev Pointer to the first limiting half-edge.
     * \param next Pointer to the second limiting half-edge.
     */
    inline void set_limits( 
			   Halfedge* prev ,
			   Halfedge* next 
			  ) 
    {
      _h_prev = prev ;
      _h_next = next ;
      
      return ;
    }


    /**
     * \fn inline void get_2d_coords( double& prev_x , double& prev_y , double& curr_x , double& curr_y , double& next_x , double& next_y ) const 
     *
     * \brief Gets the  planar coordinates of this vertex,  and of the
     * previous and next vertices.  These coordinates are defined with
     * respect to  the unfolded sequence of  triangles associated with
     * this geodesic vertex.
     *
     * \param prev_x A  reference to the X coordinate  of the previous
     * geodesic  vertex  with  respect  to the  unfolded  sequence  of
     * triangles associated with this geodesic vertex.
     * \param prev_y A  reference to the Y coordinate  of the previous
     * geodesic  vertex  with  respect  to the  unfolded  sequence  of
     * triangles associated with this geodesic vertex.
     * \param curr_x A reference ro  the X coordinate of this geodesic
     * vertex  with respect  to  its associated  unfolded sequence  of
     * triangles.
     * \param curr_y A reference to  the Y coordinate of this geodesic
     * vertex  with respect  to  its associated  unfolded sequence  of
     * triangles.
     * \param  next_x A  reference to  the  X coordinate  of the  next
     * geodesic  vertex  with  respect  to the  unfolded  sequence  of
     * triangles associated with this geodesic vertex.
     * \param  next_y A  reference to  the  Y coordinate  of the  next
     * geodesic  vertex  with  respect  to the  unfolded  sequence  of
     * triangles associated with this geodesic vertex.
     */
    inline void get_2d_coords( 
			      double& prev_x ,
			      double& prev_y ,
			      double& curr_x ,
			      double& curr_y ,
			      double& next_x ,
			      double& next_y 
			     ) 
      const 
    {
      prev_x = _prev_x ; 
      prev_y = _prev_y ;
      curr_x = _curr_x ;
      curr_y = _curr_y ;
      next_x = _next_x ;
      next_y = _next_y ;
      
      return ;
    }


    /**
     * \fn inline void set_2d_coords( double prev_x , double prev_y , double curr_x , double curr_y , double next_x , double next_y )
     *
     * \brief Assigns  planar coordinates to  this vertex, and  to the
     * previous and  next vertices.  These coordinates  are given with
     * respect to  the unfolded sequence of  triangles associated with
     * this geodesic vertex.
     *
     * \param prev_x The X  coordinate of the previous geodesic vertex
     * with respect  to the unfolded sequence  of triangles associated
     * with this geodesic vertex.
     * \param prev_y The Y  coordinate of the previous geodesic vertex
     * with respect  to the unfolded sequence  of triangles associated
     * with this geodesic vertex.
     * \param  curr_x The X  coordinate of  this geodesic  vertex with
     * respect to its associated unfolded sequence of triangles.
     * \param  curr_y The Y  coordinate of  this geodesic  vertex with
     * respect to its associated unfolded sequence of triangles.
     * \param next_x The X coordinate of the next geodesic vertex with
     * respect to  the unfolded sequence of  triangles associated with
     * this geodesic vertex.
     * \param next_y The Y coordinate of the next geodesic vertex with
     * respect to  the unfolded sequence of  triangles associated with
     * this geodesic vertex.
     */
    inline void set_2d_coords(
			      double prev_x ,
			      double prev_y ,
			      double curr_x ,
			      double curr_y ,
			      double next_x ,
			      double next_y
			     ) 
    {
      _prev_x = prev_x ;
      _prev_y = prev_y ;
      _curr_x = curr_x ;
      _curr_y = curr_y ;
      _next_x = next_x ;
      _next_y = next_y ;
      
      return ;
    }

    // ---------------------------------------------------------------
    //
    // Methods related to the Fast Marching Method.
    //
    // ---------------------------------------------------------------

    /**
     * \fn void dupdate( Vertex* v )
     *
     * \brief Update  the value  of the distance  function at  a given
     * vertex.
     *
     * \param v A vertex of the mesh.
     */
    bool dupdate( Vertex* v ) ;


    /**
     * \fn void initialize_all_attributes()
     *
     * \brief Sets the set tag of all mesh vertices.
     */
    void initialize_all_attributes() ;


    /**
     * \fn void compute_initial_front( GP* vi, std::vector< Vertex* >& li )
     *
     * \brief Sets up the initial propagation front, which consists of
     * the vertices  connected to the  initial vertex of  the geodesic
     * being computed.
     *
     * \param  vi  Pointer  to  the  initial vertex  of  the  discrete
     * geodesic.
     * \param li A reference to a list of vertices.
     */
    void compute_initial_front( GP* vi , std::vector< Vertex* >& li ) ;


    /**
     * \fn void compute_final_front( GP* ve, std::vector< Vertex* >& le , bool& alive ) const
     *
     * \brief Computes the front  around the final vertex and verifies
     * if any of the front vertices is alive.
     *
     * \param  ve  Pointer  to  the  initial vertex  of  the  discrete
     * geodesic.
     * \param le  A reference  to a list  of vertices adjacent  to the
     * final vertex.
     * param  alive Flag  to indicate  that an  alive vertex  has been
     * found.
     */
    void compute_final_front( GP* ve , std::vector< Vertex* >& le ,
      bool& alive ) const ;


    /**
     * \fn bool compute_close( std::vector< Vertex* >& li, std::vector< Vertex* >& le )
     *
     * \brief Sweeps  the initial front  over the mesh to  compute the
     * distance  function  value  in  the vertices  belonging  to  the
     * current far set.
     *
     * \param li List of vertices of the initial front.
     * \param le List of vertices adjacent ot the final vertex.
     *
     * \return Return true if this computation was done successfully.
     */
    bool compute_close( std::vector< Vertex* >& li , 
      std::vector< Vertex* >& le ) ;


    /**
     * \fn bool in_the_same_mesh_element( GP* vi, GP* ve , std::list< GP* >& lv ) const
     *
     * \brief If the initial and  final vertices of a geodesic are the
     * same  mesh element  or belong  to  the same  mesh element,  the
     * method computes  the geodesic  and returns true.  Otherwise, it
     * returns false.
     *
     * \param vi The initial vertex of the discrete geodesic.
     * \param ve The final vertex of the discrete geodesic.
     * \param lv A reference to a list of geodesic vertices.
     *
     * \return The logic value true  if the initial and final vertices
     * are  in  the same  mesh  element,  and  the logic  value  false
     * otherwise.
     */
    bool in_the_same_mesh_element( GP* vi , GP* ve , 
      std::list< GP* >& lv ) const ;


    /**
     * \fn void compute_path( GP* vi , GP* ve , std::vector< Vertex* >& li , std::vector< Vertex* >& le , std::list< GP* >& lv )
     *
     * \brief Computes the  set of mesh vertices that  forms the first
     * approximation to the discrete geodesic being computed.
     *
     * \param vi The initial vertex of the geodesic.
     * \param ve The final vertex of the geodesic.
     * \param li The set of vertices of the initial front.
     * \param le The set of vertices adjacent to the final vertex.
     * \param lv A reference to the list of geodesic vertices.
     */
    void compute_path( GP* vi , GP* ve , std::vector< Vertex* >& li ,
      std::vector< Vertex* >& le , std::list< GP* >& lv ) ;


    /**
     * \fn void get_coords( GP* v, double& x , double& y , double& z ) const
     *
     * \brief Returns  the three coordinates of a  given geodesic vertex.
     *
     * \param v Pointer to a geodesic vertex.
     * \param x A reference to the X coordinate.
     * \param y A reference to the Y coordinate.
     * \param z A reference to the Z coordinate.
     */
    void get_coords( GP* v, double& x , double& y , double& z ) const ;


    /**
     * \fn void get_restricted_adj_vertices( Vertex* v , std::vector< Vertex* >& a ) const
     *
     * \brief Finds all mesh vertices connected to a given mesh vertex
     * by an edge that is incident to a marked face.
     *
     * \param v A pointer to a geodesic vertex.
     * \param a A reference to a list of adjacent mesh vertices.
     */
    virtual void get_restricted_adj_vertices( Vertex* v , std::vector< Vertex* >& a ) const ;


    /**
     * \fn void get_restricted_adj_vertices( GP* p , std::vector< Vertex* >& a ) const
     *
     * \brief  Finds all mesh  vertices adjacent  to the  mesh element
     * that contains  a given geodesic vertex.  In  particular, if the
     * geodesic vertex is a mesh vertex, the method finds the adjacent
     * mesh  vertices.  If  the geodesic  vertex is  a mesh  edge, the
     * method finds  the vertices in  the edge star.  Finally,  if the
     * geodesic  vertex is  a mesh  face,  the method  finds the  face
     * vertices.
     *
     * \param p A pointer to a geodesic vertex.
     * \param a A reference to a list of vertices.
     */
    void get_restricted_adj_vertices(  GP* p ,
      std::vector< Vertex* >& a ) const ;


    /**
     * \fn bool is_incident_to_active_face( Vertex* v ) const
     *
     * \brief Decides whether a given vertex is incident to an active face.
     *
     * \param v A pointer to a vertex.
     *
     * \return The logic  value true if the given  vertex is incident to
     * an active face, and the logic value false otherwise.
     */
    bool is_incident_to_active_face( Vertex* v ) const ;


    // ---------------------------------------------------------------
    //
    // Methods related to the path correction algorithm.
    //
    // ---------------------------------------------------------------

    /**
     * \fn void start_shortening()
     *
     * \brief Carry on the actual geodesic shortening process.
     *
     * \sa compute_first_approximation, shorten_geodesic
     */
    void start_shortening() ;


    /**
     * \fn bool eliminate_vertices( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief  Eliminates all  geodesic vertices  that follow  a given
     * vertex  and  share  one  mesh  face with  its  predecessor  and
     * successor vertices.
     *
     * \param  prev  The  predecessor  of  the  given  vertex  in  the
     * geodesic.
     * \param node The given vertex.
     * \param next The successor of the given vertex in the geodesic.
     *
     * \return The logic value true if the end of the geodesic has not
     * been  reached during  the  elimination process,  and the  logic
     * value false otherwise.
     */
    bool eliminate_vertices( Path_iterator& prev , 
      Path_iterator& node , Path_iterator& next ) ;


    /**
     * \fn void pre_processing( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief Calculates  the error  associated with a  given geodesic
     * vertex.
     *
     * \param  prev  The  predecessor  of  the  given  vertex  in  the
     * geodesic.
     * \param node The given vertex.
     * \param next The successor of the given vertex in the geodesic.
     *
     * \return The logic value true if the current vertex has not been
     * removed from the geodesic, and the logic value false otherwise.
     */
    void pre_processing( Path_iterator& prev , Path_iterator& node ,
        Path_iterator& next ) ;

    /**
     * \fn bool pos_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief Shortens the geodesic by straightening it in the segment
     * defined  by three  given  vertices, which  are named  previous,
     * current, and  next. The current vertex belongs  to the interior
     * of an edge.
     *
     * \param prev The previous vertex.
     * \param node The current vertex.
     * \param next The next vertex.
     *
     * \return The logic value true if the geodesic has been shortened
     * by  moving  the  current   vertex  along  the  interior  of  an
     * edge. Otherwise, the logic value false is returned.
     *
     * \sa pos_on_v
     */
    bool pos_on_e( Path_iterator& prev , Path_iterator& node ,
        Path_iterator& next ) ;


    /**
     * \fn void pos_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief Shortens the geodesic by straightening it in the segment
     * defined  by three  given  vertices, which  are named  previous,
     * current,  and next. The  current vertex  coincides with  a mesh
     * vertex.
     *
     * \param prev The previous vertex.
     * \param node The current vertex.
     * \param next The next vertex.
     *
     * \sa pos_on_e
     */
    void pos_on_v( Path_iterator& prev , Path_iterator& node ,
        Path_iterator& next ) ;


    /**
     * \fn bool snap_to_vertex( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief  Snaps  a  given  geodesic  vertex, which  lies  in  the
     * interior  of a  mesh edge,  to a  mesh vertex  if  the geodesic
     * vertex and  its immediate neighbors  are too close to  the mesh
     * vertex  and if they  belong to  the concave  local side  of the
     * geodesic.
     *
     * \param prev The previous vertex.
     * \param node The current vertex.
     * \param next The next vertex.
     *
     * \return The logic  value true if the given  geodesic vertex has
     * been  snapped to  a  mesh  vertex, and  the  logic value  false
     * otherwise.
     */
    bool snap_to_vertex( Path_iterator& prev , Path_iterator& node ,
      Path_iterator& next ) ;


    /**
     * \fn bool can_be_removed( const CP& cp_prev , const CP& cp_node , const CP& cp_next ) const
     *
     * \brief Decides  whether a given geodesic vertex  can be removed
     * from the current geodesic path. The decision takes into account
     * the  predecessor and  successor of  the vertex  in  the current
     * geodesic.
     *
     * \param cp_prev Reference to the  previous vertex.
     * \param cp_node Reference to a vertex.
     * \param cp_next Reference to the next vertex.
     *
     * \return The logic value true  if a given geodesic vertex can be
     * removed, and the logic value false otherwise.
     */
    bool can_be_removed( const CP& cp_prev , const CP& cp_node , 
      const CP& cp_next ) const ;


    /**
     * \fn void update_error_info_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief Updates  the information needed for  computing the error
     * associated with  a given vertex, which belongs  to the interior
     * of a mesh edge.
     *
     * \param prev The predecessor of the given vertex in the geodesic.
     * \param node The given vertex.
     * \param next The successor of the given vertex in the geodesic.
     *
     * \sa pre_processing
     */
    void update_error_info_on_e( Path_iterator& prev , 
      Path_iterator& node , Path_iterator& next ) ;


    /**
     * \fn void compute_error_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief  Computes the  error  associated with  a given  geodesic
     * vertex that belongs to the interior of a mesh edge.
     *
     * \param prev The predecessor of the given vertex in the geodesic.
     * \param node The given vertex.
     * \param next The successor of the given vertex in the geodesic.
     *
     * \sa pre_processing, compute_error_on_v
     */
    void compute_error_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next ) ;


    /**
     * \fn void update_error_info_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next , double& langle , double& rangle )
     *
     * \brief Updates  the information needed for  computing the error
     * associated with a geodesic given vertex, which coincides with a
     * mesh vertex.
     *
     * \param  prev  The  predecessor  of  the  given  vertex  in  the
     * geodesic.
     * \param node The given vertex.
     * \param next The successor of the given vertex in the geodesic.
     * \param  langle  The angle  defined  by  the  left side  of  the
     * geodesic at the given vertex.
     * \param  rangle The  angle  defined  by the  right  side of  the
     * geodesic at the given vertex.
     *
     * \sa pre_processing
     */
    void update_error_info_on_v( Path_iterator& prev , 
      Path_iterator& node , Path_iterator& next , double& langle , 
      double& rangle ) ;


    /**
     * \fn void compute_error_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
     *
     * \brief  Computes the  error  associated with  a given  geodesic
     * vertex, which coincides with a mesh vertex.
     *
     * \param  prev  The  predecessor  of  the  given  vertex  in  the
     * geodesic.
     * \param node The given vertex.
     * \param next The successor of the given vertex in the geodesic.
     *
     * \sa pre_processing, update_error_info_on_v
     */
    void compute_error_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next ) ;


    /**
     * \fn void search_halfedges_on_e( const CP& cp_prev, const CP& cp_node, const CP& cp_next, Halfedge*& h_prev, Halfedge*& h_next ) const
     *
     * \brief Determines the two triangles  that share a mesh edge are
     * cut through by  a geodesic that passes through  the interior of
     * the edge.
     *
     * \param cp_prev The previous vertex.
     * \param  cp_node The geodesic  vertex that  belongs to  the edge
     * interior.
     * \param cp_next The next vertex.
     * \param h_prev A reference to a pointer to a half-edge of one of
     * the two  triangles. The origin  vertex of the half-edge  is one
     * end  vertex of  the  mesh edge  containing  the given  geodesic
     * vertex in its interior.
     * \param h_next  A reference to a  pointer to a  half-edge of the
     * other triangle.  The origin vertex of the half-edge is the same
     * end  vertex of  the  mesh edge  containing  the given  geodesic
     * vertex in its interior.
     */
    void search_halfedges_on_e( const CP& cp_prev , const CP& cp_node , 
      const CP& cp_next , Halfedge*& h_prev , Halfedge*& h_next ) const ;


    /**
     * \fn void search_halfedges_on_v( CP& cp_prev, CP& cp_node, CP& cp_next, Halfedge*& h_prev, Halfedge*& h_next ) const
     *
     * \brief Determines the  two triangles of a vertex  star that are
     * cut through by a geodesic  that passes through the vertex.  The
     * triangles  contain the  predecessor and  successor of  the star
     * vertex  in the  geodesic, and  they  can possibly  be the  same
     * triangle.
     *
     * \param cp_prev The previous vertex.
     * \param cp_node The star vertex.
     * \param cp_next The next vertex.
     * \param h_prev  A reference to a  pointer to a  half-edge of the
     * star  triangle that  is cut  through  by the  geodesic when  it
     * enters the star. The origin vertex of the half-edge is the star
     * vertex.
     * \param h_next  A reference to a  pointer to a  half-edge of the
     * star  triangle that  is cut  through  by the  geodesic when  it
     * leaves the star. The origin vertex of the half-edge is the star
     * vertex.
     */
    void search_halfedges_on_v( CP& cp_prev , CP& cp_node , CP& cp_next ,
      Halfedge*& h_prev, Halfedge*& h_next ) const ;


    /**
     * \fn void compute_angles_on_v( Halfedge* h1 , Halfedge* h2 , CP& cp_prev , CP& cp_node , CP& cp_next , double& langle , double& rangle ) const
     *
     * \brief Computes the  left and right angles defined  by the left
     * and  right  sides,  respectively,  of  a geodesic  at  a  given
     * geodesic vertex. This geodesic vertex is also a mesh vertex.
     *
     * \param  h1 Pointer to  a half-edge  that identifies  one vertex
     * star triangle  cut through by  the geodesic when it  enters the
     * vertex star.
     * \param  h2 Pointer to  a half-edge  that identifies  one vertex
     * star triangle  cut through by  the geodesic when it  leaves the
     * vertex star.
     * \param cp_prev A reference to the geodesic vertex that precedes
     * the star vertex in the geodesic.
     * \param  cp_node  A  reference   to  the  geodesic  vertex  that
     * coincides with the star vertex.
     * \param cp_next A reference to the geodesic vertex that succedes
     * the star vertex in the geodesic.
     * \param langle A reference to the value of the left angle.
     * \param rangle A reference to the value of the right angle.
     */
    void compute_angles_on_v( Halfedge* h1 , Halfedge* h2 , 
      CP& cp_prev , CP& cp_node , CP& cp_next , double& langle ,
      double& rangle ) const ;


    /**
     * \fn void unfold_triangle_strip( Halfedge* h1 , Halfedge* h2 )
     *
     * \brief Unfold (in  R2) a strip of consecutives  mesh faces that
     * belong to the star of a vertex.
     *
     * \param h1 A pointer to  a half-edge with origin at the star
     * vertex and contained in the first strip face.
     * \param h2 A pointer to  a half-edge with origin at the star
     * vertex and contained in the last strip face.
     */
    void unfold_triangle_strip( Halfedge* h1 , Halfedge* h2 ) ;


    /**
     * \fn void compute_2d_coords_on_e( const CP& cp_prev , const CP& cp_node , const CP& cp_next )
     *
     * \brief  Computes  the  planar  Cartesian coordinates  of  three
     * consecutive vertices  of this geodesic.   These coordinates are
     * computed  by unfolding two  mesh triangles  in the  plane.  The
     * triangles themselves  and their coordinates are  given in terms
     * of  the input  parameters.  The  second vertex  belongs  to the
     * interior of a mesh edge.
     *
     * \param cp_prev The first of the three consecutive vertices.
     * \param cp_node The second of the three consecutive vertices.
     * \param cp_next The third of the three consecutive vertices.
     */
    void compute_2d_coords_on_e( const CP& cp_prev , 
      const CP& cp_node , const CP& cp_next ) ;


    /**
     * \fn void compute_2d_coords_on_v( const CP& cp_prev , const CP& cp_node , const CP& cp_next )
     *
     * \brief  Computes  the  planar  Cartesian coordinates  of  three
     * consecutive vertices  of this geodesic.   These coordinates are
     * computed  by  unfolding  mesh  triangles  in  the  plane.   The
     * triangles themselves  and their coordinates are  given in terms
     * of  the input parameters.  The second  vertex coincides  with a
     * mesh vertex,
     *
     * \param cp_prev The first of the three consecutive vertices.
     * \param cp_node The second of the three consecutive vertices.
     * \param cp_next The third of the three consecutive vertices.
     */
    void compute_2d_coords_on_v( const CP& cp_prev , 
      const CP& cp_node , const CP& cp_next ) ;


    /**
     * \fn bool compute_intersection_on_e( CP& cp_prev , CP& cp_node , CP& cp_next ) const
     *
     * \brief  Computes the  intersection  of two  line segments.  The
     * first one  is an  unfolded mesh edge  that contains  a geodesic
     * vertex (named the current vertex). The second one is defined by
     * the unfolded geodesics vertices  that precedes and succeeds the
     * current vertex.
     *
     * \param cp_prev A reference  to the geodesic vertex the precedes
     * the vertex that belongs to the interior of a mesh edge.
     * \param cp_node  A reference to  the vertex that belongs  to the
     * interior of a mesh edge.
     * \param cp_next A reference  to the geodesic vertex the succeeds
     * the vertex that belongs to the interior of a mesh edge.
     *
     * \return The logic value true in two cases: (a) the intersection
     * point is a point in  the interior of both segments and distinct
     * from the position of the  current geodesic vertex; (b) there is
     * no  intersection  between the  interior  of  the two  segments.
     * Otherwise, the logic value false is returned.
     */
    bool compute_intersection_on_e( CP& cp_prev , CP& cp_node , 
      CP& cp_next ) const ;


    /**
     * \fn void compute_intersection_on_v( CP& cp_prev , CP& cp_node , CP& cp_next , Path& new_path  ) const
     *
     * \brief Computes the intersection  points of a line segment with
     * several other  line segments.  The first  line segment connects
     * the previous and next vertices of a given geodesic vertex.  The
     * geodesic  vertex coincided  with a  mesh vertex.  The remaining
     * segments are  edges connecting  the given geodesic  vertex with
     * other vertices of the mesh. The intersection is computed in the
     * plane, after unfolding a strip of triangles of the given vertex
     * star.
     *
     * \param cp_prev A reference  to the geodesic vertex the precedes
     * the given star vertex.
     * \param cp_node A reference to the given star vertex.
     * \param cp_next A reference  to the geodesic vertex the succeeds
     * the given star vertex.
     * \param  new_path A reference  to a  list with  the intersection
     * points.
     */
    void compute_intersection_on_v( CP& cp_prev , CP& cp_node , 
      CP& cp_next , Path &new_path ) const ;


    /**
     * \fn void find_common_faces( Vertex* v1 , Vertex* v2 , Face*& f1 , Face*& f2 ) const
     *
     * \brief Finds  the two faces  incident to two  distinct vertices
     * (if any).
     *
     * \param v1 Pointer to a mesh vertex.
     * \param v2 Pointer to a mesh vertex.
     * \param f1 Reference to a pointer to a mesh face.
     * \param f2 Reference to a pointer to a mesh face.
     */
    void find_common_faces( Vertex* v1 , Vertex* v2 , Face*& f1 , 
      Face*& f2 ) const ;


    /**
     * \fn void find_common_faces( Vertex* v , Edge* e , Face*& f1 , Face*& f2 ) const
     *
     * \brief Finds  the face(s)  incident to a  mesh edge and  a mesh
     * vertex (if any).
     *
     * \param e Pointer to a mesh edge.
     * \param v Pointer to a mesh vertex.
     * \param f1 Reference to a pointer to a mesh face.
     * \param f2 Reference to a pointer to a mesh face.
     */
    void find_common_faces( Vertex* v , Edge* e , Face*& f1 , Face*& f2 ) const ;


    /**
     * \fn void find_common_faces( Edge* e1 , Edge* e2 , Face*& f1 , Face*& f2 ) const
     *
     * \brief Finds the face(s) incident to two mesh edges (if any).
     *
     * \param e1 Pointer to a mesh edge.
     * \param e2 Pointer to a mesh edge.
     * \param f1 Reference to a pointer to a mesh face.
     * \param f2 Reference to a pointer to a mesh face.
     */
    void find_common_faces( Edge* e1 , Edge* e2 , Face*& f1 , Face*& f2 ) const ;


    /**
     * \fn bool is_in_face( Vertex* v , Face* f ) const
     *
     * \brief Decides whether  a given mesh vertex belongs  to a given
     * mesh face.
     *
     * \param v Pointer to a mesh vertex.
     * \param f Pointer to a mesh face.
     *
     * \return The logic value true if the given vertex belongs to the
     * given face, and the logic value false otherwise.
     */
    bool is_in_face( Vertex* v , Face* f ) const ;


    /**
     * \fn bool is_in_face( Edge* e , Face* f ) const
     *
     * \brief Decides  whether a  given mesh edge  belongs to  a given
     * mesh face.
     *
     * \param e Pointer to a mesh edge.
     * \param f Pointer to a mesh face.
     *
     * \return The logic  value true if the given  edge belongs to the
     * given face, and the logic value false otherwise.
     */
    bool is_in_face( Edge* e , Face* f ) const ;


    /**
     * \fn double angle( Vertex* v1 , Vertex* v2 , Vertex* v3 ) const
     *
     * \brief Calculates the  angle defined at the first  of the three
     * given mesh vertices of a mesh triangle.
     *
     * \param v1 The first vertex of a mesh triangle.
     * \param v2 The second vertex of a mesh triangle.
     * \param v3 The third vertex of a mesh triangle.
     *
     * \return The angle at v1 defined by  the vectors ( v2 - v1 ) and
     * ( v3 - v1 ).
     */
    double angle( Vertex* v1 , Vertex* v2 , Vertex* v3 ) const ;


    /**
     * \fn double angle( Vertex* v1 , Vertex* v3 , CP& cp ) const
     *
     * \brief Calculates the  angle defined at the first  of the three
     * given  geodesic  vertices.  The  first  two  vertices are  mesh
     * vertices, and  the third one may  be either a mesh  vertex or a
     * point in the interior of a mesh edge.
     *
     * \param v1 The first geodesic vertex.
     * \param v2 The second geodesic vertex.
     * \param cp The third geodesic vertex.
     *
     * \return The angle at v1 defined by  the vectors ( v2 - v1 ) and
     * ( cp - v1 ).
     */
    double angle( Vertex* v1 , Vertex* v2 , CP& cp ) const ;


    /**
     * \fn void unfold_first_triangle( Halfedge* h , double& xb , double& xc , double& yc ) const
     *
     * \brief Unfold a mesh face  in the plane. The face is identified
     * by one of  its half-edges. The origin vertex  of this half-edge
     * is mapped to the point (0,0), the following vertex is mapped to
     * a point  in the  X axis, and  the third  vertex is mapped  to a
     * point in the upper-half of the plane.
     *
     * \param h  A pointer to a half-edge.
     * \param xb  A reference  to the X  coordinate of  the second
     * Cartesian vertex of the unfolded mesh face.
     * \param  xc A  reference to  the X  coordinate of  the third
     * Cartesian vertex of the unfolded mesh face.
     * \param  yc A  reference to  the Y  coordinate of  the third
     * Cartesian vertex of the unfolded mesh face.
     */
    void unfold_first_triangle( Halfedge* h , double& xb ,
      double& xc , double& yc ) const ;


    /**
     * \fn void unfold_next_triangle( Halfedge* h , double xb , double yb , double& xc , double& yc ) const
     *
     * \brief Unfold a mesh face  in the plane. The face is identified
     * by one of  its half-edges. The origin vertex  of this half-edge
     * is mapped to the point (0,0), the following vertex is mapped to
     * a  given point,  and  the third  vertex  is mapped  to a  point
     * computed by this method.
     *
     * \param h  A pointer to a half-edge.
     * \param xb  The X coordinate of the  second Cartesian vertex
     * of the unfolded mesh face.
     * \param yb  The Y coordinate of the  second Cartesian vertex
     * of the unfolded mesh face.
     * \param xc  A reference to the  X coordinate of  the third Cartesian
     * vertex of the unfolded mesh face.
     * \param yc  A reference to the  Y coordinate of  the third Cartesian
     * vertex of the unfolded mesh face.
     */
    void unfold_next_triangle( Halfedge* h , double xb , double yb ,
        double& xc , double& yc ) const ;


    /**
     * \fn void point_on_edge_2d( const CP& cp , Halfedge* h , double x1 , double y1 , double x2 , double y2 , double& x , double& y ) const
     *
     * \brief  Computes the  planar coordinates  of a  geodesic vertex
     * that belongs to the interior  of a mesh edge. These coordinates
     * are given with respect to  the unfold of the two mesh triangles
     * that share  the mesh edge containing the  geodesic vertex.  All
     * information regarding the unfolded triangles are already cached
     * in the geodesic vertex.
     *
     * \param cp A reference to a geodesic vertex.
     * \param h Pointer to a half-edge of the mesh edge containing the
     * geodesic vertex.
     * \param  x1 The  X planar  coordinate of  one end  point  of the
     * unfolded mesh edge containing the geodesic vertex.
     * \param  y1 The  Y planar  coordinate of  one end  point  of the
     * unfolded mesh edge containing the geodesic vertex.
     * \param x2 The X planar coordinate of the other end point of the
     * unfolded mesh edge containing the geodesic vertex.
     * \param y2 The Y planar coordinate of the other end point of the
     * unfolded mesh edge containing the geodesic vertex.
     * \param x  A reference  to the X  planar coordinate  of geodesic
     * vertex.
     * \param y  A reference  to the Y  planar coordinate  of geodesic
     * vertex.
     */
    void point_on_edge_2d( const CP & cp , Halfedge* h , double x1 , 
      double y1 , double x2 , double y2 , double& x , double& y ) 
      const ;


    /**
     * \fn void point_on_face_2d( const CP& cp , Halfedge* h , double x1 , double y1 , double x2 , double y2 , double x3 , double y3 , double& x , double& y ) const
     *
     * \brief  Computes the  planar coordinates  of a  geodesic vertex
     * that belongs to the interior  of a mesh face. These coordinates
     * are given  with respect  to the unfold  of the mesh  face.  All
     * information  regarding  the unfolded  face  are  cached in  the
     * geodesic vertex.
     *
     * \param cp A reference to a geodesic vertex.
     * \param h Pointer  to the halfedge whose origin  is the geodesic
     * vertex containing the cached unfolded face information.
     * \param x1  The X planar coordinate  of one vertex  of mesh face
     * containing the geodesic vertex.
     * \param y1 The Y planar coordinate of of one vertex of mesh face
     * containing the geodesic vertex.
     * \param x2  The X  planar coordinate of  a second vertex  of the
     * unfolded mesh face containing the geodesic vertex.
     * \param y2  The Y  planar coordinate of  a second vertex  of the
     * unfolded mesh face containing the geodesic vertex.
     * \param x3  The X  planar coordinate of  a second vertex  of the
     * unfolded mesh face containing the geodesic vertex.
     * \param y3  The Y  planar coordinate of  a second vertex  of the
     * unfolded mesh face containing the geodesic vertex.
     * \param x  A reference  to the X  planar coordinate  of geodesic
     * vertex.
     * \param y  A reference  to the Y  planar coordinate  of geodesic
     * vertex.
     */
    void point_on_face_2d( const CP& cp , Halfedge* h , double x1 , 
      double y1 , double x2 , double y2 , double x3 , double y3 , 
      double& x , double& y ) const ;


    /**
     * \fn bool is_consistent() const
     *
     * \brief Checks whether the  current geodesic is consistent.  The
     * geodesic  is   considered  consistent  if  each   pair  of  its
     * consecutive vertices are incident to  the same face, but not to
     * the same edge or vertex.
     *
     * \return The logic value true if the geodesic is consistent, and
     * the logic value false otherwise.
     */
    bool is_consistent() const ;


    // ---------------------------------------------------------------
    //
    // Protected data members.
    //
    // ---------------------------------------------------------------

  protected:

    const double _error ; ///< Upper bound  for the maximum error associated with any  geodesic  vertex.

    const unsigned _max_iter ; ///< Upper bound for  the number of iterations  of the path correction algorithm.

    Path _cp_path ; ///< The current geodesic path.

    MI* _mesh ; ///< Pointer to the mesh over which the geodesic is defined.

    unsigned _num_iter ;  ///< Total number of iterations of the path correction algorithm.

    std::vector< double > _pts ;  ///< Planar vertex coordinates of the unfolded triangle sequence corresponding to the current vertex of this geodesic.

    Halfedge* _h_prev ;     ///< Pointer to the first limiting half-edge of the current vertex of this geodesic.

    Halfedge* _h_next ;     ///< Pointer to the second limiting half-edge of the current vertex of this geodesic.

    double _prev_x ;       ///< Planar vertex X coordinate of the previous vertex preceding the current vertex of this geodesic.

    double _prev_y ;       ///< Planar vertex Y coordinate of the previous vertex preceding the current vertex of this geodesic.

    double _curr_x ;       ///< Planar vertex X coordinate of the current geodesic vertex.

    double _curr_y ;       ///< Planar vertex Y coordinate of the current geodesic vertex.

    double _next_x ;      ///< Planar vertex X coordinate of the vertex succeeding the current vertex of this geodesic.

    double _next_y ;      ///< Planar vertex Y coordinate of the vertex succeeding the current vertex of this geodesic.

  } ;


  /**
   * \fn void Approx_geodesics< Mesh >::compute( GP* v1 , GP* v2 , std::list< GP* >& lv )
   *
   * \brief  Computes a discrete  geodesic over  a triangle  mesh. The
   * initial and final vertices of the geodesic can be vertices of the
   * mesh, or points in the interior of edges or faces of the mesh.
   *
   * \param v1 The initial point of the discrete geodesic.
   * \param v2 The final point of the discrete geodesic.
   * \param lv A reference to a list of geodesic vertices.
   *
   * \return A  polygonal curve  consisting of the  geodesic vertices,
   * which are vertices of the mesh and/or points in the interior of a
   * mesh edge or face.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::compute( 
				    GP* v1 , 
				    GP* v2 , 
				    std::list< GP* >& lv 
				   ) 
  {
    /* ---------------------------------------------------------------
     *
     * The implementation of this  method follows the algorithm in the
     * paper:
     *
     * Dimas Martínez, Luiz Velho, and Paulo C. Carvalho
     * Computing geodesics on triangular meshes
     * Computer and Graphics, 29(5):667–675, 2005
     *
     * ----------------------------------------------------------------
     */

    /*
     * Make sure the output list of vertices is empty.
     */
    assert( lv.empty() ) ;

    /*
     * Compute the first approximation to the discrete geodesic.
     */
    compute_first_approximation( v1 , v2 , lv ) ;

    /*
     * Apply the path correction algorithm to fix the path.
     */
    shorten_geodesic( lv ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_first_approximation( GP* vi, GP* ve , std::list< GP* >& lv )
   *
   * \brief Computes  a fast, but  rough approximation for  a discrete
   * geodesic  defined over a  triangle mesh.   The initial  and final
   * vertices of the geodesic can be vertices of the mesh the geodesic
   * is defined on, or points in the interior of edges or faces of the
   * mesh.
   *
   * \param vi The initial vertex of the discrete geodesic.
   * \param ve The final vertex of the discrete geodesic.
   * \param lv A reference to a list of geodesic vertices.
   *
   * \return A polygonal curve  consisting of geodesic vertices, which
   * are vertices of the mesh (there is vertex inside edges or faces).
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::compute_first_approximation(
							GP* vi ,
							GP* ve ,
							std::list< GP* >& lv 
						       )
  {
    /*
     * Make all mesh faces active.
     */
    for ( 
	 FaceIterator fit = _mesh->faces_begin() ; 
	 !_mesh->is_done( fit ) ; 
	 _mesh->move_forward( fit ) 
	)
    {
      _mesh->activate( *fit ) ; 
    }

    /*
     * Compute the first approximation considering the entire mesh.
     */
    compute_restricted_first_approximation( vi , ve , lv ) ;

    /*
     * Make all mesh faces inactive.
     */
    for ( 
	 FaceIterator fit = _mesh->faces_begin() ; 
	 !_mesh->is_done( fit ) ; 
	 _mesh->move_forward( fit ) 
	)
    {
      _mesh->deactivate( *fit ) ; 
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_restricted_first_approximation( GP* vi, GP* ve , std::list< GP* >& lv )
   *
   * \brief Computes  a fast, but  rough approximation for  a discrete
   * geodesic defined over  a restricted set of faces  of the triangle
   * mesh. The  restricted faces are the active  one.  Face activation
   * must be performed before this  method is invoked. The initial and
   * final vertices  of the geodesic can  be vertices of  the mesh the
   * geodesic is  defined on,  or points in  the interior of  edges or
   * faces of the mesh.
   *
   * \param vi The initial vertex of the discrete geodesic.
   * \param ve The final vertex of the discrete geodesic.
   * \param lv A reference to a list of geodesic vertices.
   *
   * \return A polygonal curve  consisting of geodesic vertices, which
   * are vertices of the mesh (there is vertex inside edges or faces).
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::compute_restricted_first_approximation(
								   GP* vi ,
								   GP* ve ,
								   std::list< GP* >& lv 
								  )
  {
#ifdef DEBUGMODE
    assert( ( vi != 0 ) && ( ve != 0 ) ) ;
    assert( vi->on_vertex() || vi->on_edge() || vi->on_face() );
    assert( ve->on_vertex() || ve->on_edge() || ve->on_face() );
    assert( lv.empty() ) ;
#endif

    /*
     * If the initial  and final points are the  same, then place only
     * one of them in the list  of geodesic vertices ( the geodesic is
     * only one vertex) .
     */
    if ( vi == ve ) {
      lv.push_back( vi ) ;
      return ;
    }

    /*
     * If  the  initial  and  final  points  do  not  share  the  same
     * "pointers", but they  are the same vertices in  the geodesic or
     * they belong  to the  same mesh element,  we cannot  compute the
     * geodesic.
     */
    if ( in_the_same_mesh_element( vi , ve , lv ) ) {
      return ;
    }

    /*
     * This procedure  uses the Fast Marching Method  (FMM) to compute
     * the values of  a distance function at the  vertices of the mesh
     * (the solution of the  Eikonal equation).  During its execution,
     * the  procedure separates  the  mesh vertices  into three  sets:
     * alive, close,  and far.  So,  its first step  initializes those
     * three sets  by tagging the  vertices with values  that identify
     * the sets they belong to.
     */

    /*
     * Initialize the alive, close, and far sets.
     */
    initialize_all_attributes() ;

    /*
     * The procedure  sweeps the mesh simulating  a front propagation.
     * As the  front propagates,  the procedure computes  the distance
     * function  values  at  the  vertices. A  vertex  whose  distance
     * function value  has been fixed is  in the alive  set.  A vertex
     * that has been  swept by the front, but  whose distance function
     * value has not been fixed is in the close set. A vertex that has
     * not been  visited is in  the far set.  Initially,  all vertices
     * adjacent to the initial vertex are placed in the alive set, and
     * the distance  value at those  vertices is set to  the Euclidean
     * distance  from the  vertex  to the  initial  vertex (i.e.,  the
     * length of  the edge connecting  them). Those vertices  form the
     * initial front.
     */

    /*
     * Compute the initial front. 
     */
    std::vector< Vertex* > li ;
    compute_initial_front( vi , li ) ;

    /* 
     * Include the final vertex  into the polygonal curve defining the
     * geodesic.
     */
    lv.push_front( ve ) ;

    /*
     * If  the final  vertex is  a vertex  of the  mesh, verify  if it
     * belongs to the  initial front. If so, there  is nothing else to
     * be done.
     */
    bool adjacent = false ;
    if ( ve->on_vertex() ) {
      /*
       * Look for the final vertex in the initial front.
       */
      GPV* gpve = dynamic_cast< GPV* >( ve ) ;

#ifdef DEBUGMODE
      assert( gpve != 0 ) ;
#endif

      typename std::vector< Vertex* >::iterator it = find( 
							  li.begin() ,
							  li.end() ,
							  gpve->get_vertex()
							 ) ;
      adjacent = ( it != li.end() ) ;
    }

    /*
     * If the final vertex is not  a mesh vertex nor is it adjacent to
     * the initial  front, we must  compute its adjacent  vertices and
     * verify if any of them belongs to the alive set. If so, there is
     * nothing else  to do. Otherwise, we start  the front propagation
     * process.
     */
    if ( !adjacent ) {
      bool alive = false ;
      std::vector< Vertex* > le ;
      compute_final_front( ve , le , alive ) ;

      /*
       * The value  of adjacent  is false if  and only if  no adjacent
       * vertex of  the final vertex  is alive. This means  that those
       * vertices  must go to  the close  set (if  they are  not there
       * already).
       */
      if ( !alive ) {
	if ( !compute_close( li , le ) ) {
          return ;
        }
      }

      /*
       * Compute the geodesic path.
       */
      compute_path( vi , ve , li , le , lv ) ;
    }

    /*
     * Insert the first vertex into the list of geodesic vertices.
     */
    lv.push_front( vi ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::dupdate( Vertex* v )
   *
   * \brief  Update the  value of  the  distance function  at a  given
   * vertex.
   *
   * \param v A vertex of the mesh.
   */
  template< typename Mesh >
  bool 
  Approx_geodesics< Mesh >::dupdate( Vertex* v ) 
  {
    /* ---------------------------------------------------------------
     *
     * The calculations below can be found in the paper
     *
     * R. Kimmel and J. A. Sethian
     * Computing geodesic paths on manifolds
     * Proceedings of the National Academy of Science (PNAS)
     * vol. 95, pp. 8431-8435,
     * July 1998
     *
     * ----------------------------------------------------------------
     */

    /*
     * Get the vertices in the star of vertex v.
     */
    std::vector< Vertex* > lverts ;
    get_restricted_adj_vertices( v , lverts ) ;

    /*
     * Set up an iterator to this vertex.
     */
    typename std::vector< Vertex* >::iterator itv = lverts.begin() ;

    /*
     * Get the first vertex of the star of vertex v.
     */
    Vertex* v_ult = *itv ;

    /*
     * Get the second vertex of the star of vertex v.
     */
    ++itv ;
    Vertex* v_cur = *itv ;

    /*
     * Traverse the remaining vertices of the star of v.
     */
    do {
      /* 
       * If  two consecutive  vertices of  the  star of  v are  alive,
       * update the value  of the distance function at  v with respect
       * to those two consecutive vertices.
       */
      if ( _mesh->is_alive( v_ult ) && _mesh->is_alive( v_cur ) ) {
        /*
         * Find  out  which  of  the  two vertices  has  the  smallest
         * distance function value.
         */
        Vertex* va ;
        Vertex* vb ;
        if ( _mesh->get_dist( v_ult ) > _mesh->get_dist( v_cur ) ) {
          va = v_cur ;
          vb = v_ult ;
        } 
	else {
          vb = v_cur ;
          va = v_ult ;
        }

        /*
         * Get the  coordinates of the vertex va,  which is associated
         * with a  distance function value  no larger than the  one of
         * vertex vb.
         */
        double xa , ya , za ;
        _mesh->get_coords( va , xa , ya , za ) ;

        /*
         * Get the  coordinates of the vertex vb,  which is associated
         * with a distance  function value no smaller than  the one of
         * vertex va.
         */
        double xb , yb , zb ;
        _mesh->get_coords( vb , xb , yb , zb ) ;

        /*
         * Get the coordinates of vertex v.
         */
        double xv , yv , zv ;
        _mesh->get_coords( v , xv , yv , zv ) ;

        /*
         * Compute the coordinates of the vectors  ac = ( va - v ) and
         * bc = ( vb - v ).
         */
        double ac_x = ( xa - xv ) ;
        double ac_y = ( ya - yv ) ;
        double ac_z = ( za - zv ) ;

        double bc_x = ( xb - xv ) ;
        double bc_y = ( yb - yv ) ;
        double bc_z = ( zb - zv ) ;

        /*
         * Compute the difference of distance function value at va and
         * vb.
         */
        double u = _mesh->get_dist( va ) - _mesh->get_dist( vb ) ;

        /*
         * Compute the distance from v to va.
         */
        const double a = sqrt( 
			      ( bc_x * bc_x ) + 
			      ( bc_y * bc_y ) + 
			      ( bc_z * bc_z ) 
			     ) ;

        /*
         * Compute the distance from v to vb.
         */
        const double b = sqrt( 
			      ( ac_x * ac_x ) + 
			      ( ac_y * ac_y ) + 
			      ( ac_z * ac_z ) 
			     ) ;

        /*
         * Compute the cosine of the angle defined at v by the vectors
         * (va - v) and (vb -  v). Note that those vectors need not be
         * unit vectors.
         */
        double dot_prod = ( bc_x * ac_x ) + ( bc_y * ac_y ) + 
	  ( bc_z * ac_z ) ;

        double cost = ( dot_prod ) / ( a * b ) ;

        /*
         * Compute the sine of the same angle.
         */
        double sint = sin( acos( cost ) ) ;

        /*
         * Now, we must solve  a quadratic equation whose coefficients
         * are computed below. The  solution for this equation will be
         * used to  compute the distance  function value at  the given
         * vertex v.
         */
        const double AA = ( a * a ) + ( b * b ) - 2 * a * b * cost ;

        const double BB = 2 * b * u * ( a * cost - b ) ;

        const double CC = ( b * b )
            * ( ( u * u ) - ( a * a ) * ( sint * sint ) ) ;

        /*
         * Solve AA x² + BB x + CC = 0
         */
        double delta = ( BB * BB ) - 4.0 * AA * CC ;

        /*
         * Take care of loss of significance problem.
         */
        double t0 ;
        double t1 ;

        if ( delta > 0 ) {
          delta = sqrt( delta ) ;
          if ( BB > 0 ) {
            t0 = - ( BB + delta ) / ( 2 * AA ) ;
            t1 = - ( 2 * CC ) / ( BB + delta ) ;
          } 
	  else {
            t0 = ( -BB + delta ) / ( 2 * AA ) ;
            t1 = ( 2 * CC ) / ( -BB + delta ) ;
          }
        } 
	else if ( delta == 0 ) {
          t0 = ( -BB ) / ( 2.0 * AA ) ;
          t1 = t0 ;
        } 
	else {
          t0 = u ;
          t1 = t0 ;
        }

        /*
         * The distance function value "t" must satisfy u < t
         *
         * So, we pick the  value in { t0 , t1 }  that is smaller than
         * u. If both values are,  then we choose the smaller. If none
         * are, we let t be u.
         */
        double t ;

        if ( u < t0 ) {
          t = std::min( t0 , t1 ) ;
        } 
	else if ( u < t1 ) {
          t = t1 ;
        }
	else {
          t = u ;
        }

        /*
         * Re-use variable t0:
         */
        t0 = ( b * ( t - u ) ) / t ;

        /*
         * Set the distance value at v.
         */
        if ( ( u < t ) && ( ( a * cost ) < t0 ) && ( t0 < ( a / cost ) ) ) {
          _mesh->set_dist( 
			  v , 
			  std::max( 
				   _mesh->get_dist( v ) , 
				   -t + _mesh->get_dist( va ) 
				  ) 
			 ) ;
        } 
	else {
          _mesh->set_dist( 
			  v , 
			  std::max(
				   _mesh->get_dist( v ) ,
				   std::max(
					    -b + _mesh->get_dist( va ) ,
					    -a + _mesh->get_dist( vb ) 
					   ) 
				  ) 
			 ) ;
        }
      }
      // endif is_alive( va ) && is_alive( vb )

      /*
       * Get the next pair of consecutive vertices of the star of v.
       */
      v_ult = v_cur ;

      ++itv ;
      
      if ( itv == lverts.end() ) {
        itv = lverts.begin() ;
      }
      
      v_cur = *itv ;
    } 
    while ( v_ult != * ( lverts.begin() ) ) ;
    
    return ( _mesh->get_dist( v ) > -HUGE_VAL ) ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::initialize_all_attributes()
   *
   * \brief Sets the set tag of all mesh vertices.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::initialize_all_attributes()
  {
    for ( 
	 VertexIterator itv = _mesh->vertices_begin() ; 
	 !_mesh->is_done( itv ) ; 
	 _mesh->move_forward( itv ) 
	) 
    {
      Vertex* v = _mesh->get_vertex( itv ) ;
      if ( is_incident_to_active_face( v ) ) {
	_mesh->set_alive( v , false ) ;
	_mesh->set_dist( v , -HUGE_VAL ) ;
	_mesh->set_close( v , false ) ;
      }
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_initial_front( GP* vi, std::vector< Vertex* >& li )
   *
   * \brief Sets  up the initial propagation front,  which consists of
   * the  vertices connected  to the  initial vertex  of  the geodesic
   * being computed.
   *
   * \param vi Pointer to the initial vertex of the discrete geodesic.
   * \param li A reference to a list of vertices.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::compute_initial_front(
						  GP* vi , 
						  std::vector< Vertex* >& li 
						 )
  {
    /*
     * Get the coordinates of the initial vertex.
     */
    double px , py , pz ;
    get_coords( vi , px , py , pz ) ;

    /*
     * Get the list of vertices connected to the initial vertex of the
     * geodesic.
     */
    std::vector< Vertex* > lverts ;
    get_restricted_adj_vertices( vi , lverts ) ;

    /*
     * Loop over the  vertices connected to the initial  vertex of the
     * geodesic.
     */
    for ( typename std::vector< Vertex* >::iterator it = 
	  lverts.begin() ; it != lverts.end() ; ++it ) {
      /*
       * Get a vertex connected to the initial vertex.
       */
      Vertex* pV = *it ;

      /*
       * Tag the vertex as a vertex in the alive set.
       */
      _mesh->set_alive( pV , true ) ;
      
      /*
       * Get the vertex coordinates.
       */
      double ppV_x , ppV_y , ppV_z ;
      _mesh->get_coords( pV , ppV_x , ppV_y , ppV_z ) ;
      
      /*
       * Compute its distance from the initial vertex.
       */
      double vx = px - ppV_x ;
      double vy = py - ppV_y ;
      double vz = pz - ppV_z ;

      double length = sqrt( 
			   ( vx * vx ) + 
			   ( vy * vy ) + 
			   ( vz * vz ) 
			  ) ;

      /*
       * Set the distance function value.
       */
      _mesh->set_dist( pV , -length ) ;
      
      /*
       * Insert the  vertex into the  list of vertices of  the initial
       * front.
       */
      li.push_back( pV ) ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_final_front( GP* ve, std::vector< Vertex* >& le , bool& alive ) const
   *
   * \brief Computes the front around the final vertex and verifies if
   * any of the front vertices is alive.
   *
   * \param ve Pointer to the initial vertex of the discrete geodesic.
   * \param le A reference to a list of vertices adjacent to the final
   * vertex.
   * \param  alive Flag  to indicate  that  an alive  vertex has  been
   * found.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::compute_final_front(
						GP* ve , 
						std::vector< Vertex* >& le ,
						bool& alive
					       )
    const
  {
    /*
     * Get the vertices adjacent to the final vertex.
     */ 
    get_restricted_adj_vertices( ve , le ) ;

    for ( typename std::vector< Vertex* >::iterator it = le.begin() ; 
	  it != le.end() ; ++it ) {
      /*
       * Pick one vertex adjacent to the final vertex.
       */
      Vertex* v = *it ;
      
      /*
       * If this vertex is already  in the alive set, set the adjacent
       * flag to  true, and stop  the loop. Otherwise,  keep iterating
       * the loop.
       */
      if ( _mesh->is_alive( v ) ) {
	alive = true ;
	break ;
      }
    }

    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::compute_close( std::vector< Vertex* >& li, std::vector< Vertex* >& le )
   *
   * \brief  Sweeps the  initial front  over the  mesh to  compute the
   * distance function value in  the vertices belonging to the current
   * far set.
   *
   * \param li List of vertices of the initial front.
   * \param le List of vertices adjacent ot the final vertex.
   *
   * \return Return true if this computation was done successfully.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::compute_close( 
					  std::vector< Vertex* >& li ,
					  std::vector< Vertex* >& le 
					 ) 
  {
    /* 
     * Defines  a heap  to get  the close  set vertex  associated with
     * smallest function  value by the time the  vertices are selected
     * to go to the alive set.
     */
    SortedSet close ;

    /*
     * Consider one vertex of the initial front at a time.
     */
    for ( unsigned i = 0 ; i < li.size() ; i++ ) {
      /*
       * Get the  adjacent vertices of the initial  front vertex under
       * consideration.
       */
      std::vector< Vertex* > lverts ;
      get_restricted_adj_vertices( li[ i ] , lverts ) ;

      /*
       * For each vertex  in the list of adjacent  vertices, compute a
       * distance value.
       */
      for ( typename std::vector< Vertex* >::iterator it = 
	      lverts.begin() ; it != lverts.end() ; ++it ) {
        /*
         * Get an adjacent vertex.
         */
        Vertex* v1 = *it ;

        /*
	 * If the vertex is not in the alive set nor in the close set,
         * compute a distance function  value for it.  This value must
         * come from two vertices in the alive set, which are adjacent
         * to v1.
	 */
        if ( !_mesh->is_alive( v1 ) && !_mesh->is_close( v1 ) ) {
          /*
           * Compute the value of the distance function at v1.
           */
          if ( dupdate( v1 ) ) {
            /*
	     * If  the  value  was  successfully computed,  place  the
             * vertex in the close set.
	     */
            close.push( v1 , _mesh->get_dist( v1 ) ) ;
            _mesh->set_close( v1 , true ) ;
          }
        }
      }
    }

    /*
     * Get  the vertex  of the  close set  with the  smallest function
     * value.
     */
    Vertex* pV = close.top() ;

    /*
     * Remove it from the close set.
     */
    _mesh->set_close( pV , false ) ;
    close.pop() ;

    /*
     * Look for  the vertex  in the list  of adjacent vertices  of the
     * final geodesic vertex. If it is there, there is nothing else to
     * be done.
     */
    typename std::vector< Vertex* >::iterator it = 
      find( le.begin() , le.end() , pV ) ;

    /* 
     * While the vertex is not  adjacent to the final vertex, place it
     * into the  alive set, and propagate  the front up  to sweep more
     * mesh vertices.
     */
    while ( ( pV != 0 ) && ( it == le.end() ) ) {
      /*
       * Place the vertex into the alive set.
       */
      _mesh->set_alive( pV , true ) ;

      /*
       * Get its adjacent vertices. 
       */
      std::vector< Vertex* > lverts ;
      get_restricted_adj_vertices( pV , lverts ) ;

      /*
       * Consider one adjacent vertex at a time.
       */
      for ( typename std::vector< Vertex* >::iterator itv = 
	      lverts.begin() ; itv != lverts.end() ; ++itv ) {
	/*
	 * Get the current vertex.
	 */
        Vertex* v2 = *itv ;

        /*
	 * If the vertex is not in the alive set nor in the close set,
         * compute a distance function  value for it.  This value must
         * come from two vertices in the alive set, which are adjacent
         * to v2.
	 */
        if ( !_mesh->is_alive( v2 ) && !_mesh->is_close( v2 ) ) {
          /*
	   * If the value was  successfully computed, place the vertex
           * in the close set.
	   */
          if ( dupdate( v2 ) ) {
            /*
	     * If  the  value  was  successfully computed,  place  the
             * vertex in the close set.
	     */
            close.push( v2 , _mesh->get_dist( v2 ) ) ;
            _mesh->set_close( v2 , true ) ;
          }
        }
      }

      /*
       * Get  the next  close set  vertex with  the  smallest distance
       * function value.
       */
      if ( close.empty() ) {
        return false ;
      }

      pV = close.top() ;

      /*
       * Remove it from the close set.
       */
      _mesh->set_close( pV , false ) ;
      close.pop() ;

      /*
       * Look for this  vertex in the set of  vertices adjacent to the
       * final vertex.
       */
      it = find( le.begin() , le.end() , pV ) ;
    }

    /*
     * At  this point,  either  the close  set  is empty  or a  vertex
     * adjacent  to  the  final  vertex  has been  reached,  or  both.
     * Ideally,  only   the  second  and  the   third  situations  are
     * desirable.
     */

    close.clear() ;

    return true ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_path( GP* vi , GP* ve , std::vector< Vertex* >& li , std::vector< Vertex* >& le  , std::list< GP* >& lv )
   *
   * \brief Computes  the set  of mesh vertices  that forms  the first
   * approximation to the discrete geodesic being computed.
   *
   * \param vi The initial vertex of the geodesic.
   * \param ve The final vertex of the geodesic.
   * \param li The set of vertices of the initial front.
   * \param le The set of vertices adjacent to the final vertex.
   * \param lv A reference to the list of geodesic vertices.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::compute_path( 
					 GP* vi , 
					 GP* ve ,
					 std::vector< Vertex* >& li ,
					 std::vector< Vertex* >& le ,
					 std::list< GP* >& lv
					)
  {
      assert( vi );
      assert( ve );
    /*
     * Pick the vertex adjacent to  the final vertex with the smallest
     * distance function value. Make this vertex the current vertex.
     */
    Vertex* curr = NULL ;
    double dist = -HUGE_VAL ;

    for ( unsigned int i = 0 ; i < le.size() ; i++ ) {
      if ( _mesh->get_dist( le[i] ) > dist ) {
        dist = _mesh->get_dist( le[i] ) ;
        curr = le[i] ;
      }
    }

    /*
     * Find the vertices of the geodesic.
     */ 
    typename std::vector< Vertex* >::iterator it ;
    
    do {
      /*
       * Look for the current vertex in the initial propagation front.
       */
      it = find( li.begin() , li.end() , curr ) ;
      
      /*
       * Remove the current vertex from the alive set.
       */
      _mesh->set_alive( curr , false ) ;
      
      /*
       * Insert the current vertex into the set of geodesic vertices.
       */
      lv.push_front( new GPV( curr ) ) ;

      /*
       * If the initial front has not been reached, keep expanding the
       * final front. Otherwise, there is nothing else to be done.
       */
      if ( it == li.end() ) {
	/*
	 * Get the adjacent vertices of the current vertex.
	 */ 
	double dd = HUGE_VAL ;
	Vertex* minim = 0 ;
      
	std::vector< Vertex* > lverts ;
	get_restricted_adj_vertices( curr , lverts ) ;
      
	/*
	 * For each of the adjacent  vertices of the current vertex that
	 * belongs to the alive set, pick the one with smallest distance
	 * value.
	 */
	for ( typename std::vector< Vertex* >::iterator it1 = 
		lverts.begin() ; it1 != lverts.end() ; ++it1 ) {
	  Vertex* aux = *it1 ;
	  if ( _mesh->is_alive( aux ) ) {
	    double dd1 = -_mesh->get_dist( aux ) ;
	    if ( dd1 < dd ) {
	      dd = dd1 ;
	      minim = aux ;
	    }
	  }
	}
      
	/*
	 * There must be one vertex  adjacent to the current that also
	 * belongs to  the alive set. If  not, we have  a problem with
	 * the algorithm.
	 */
#ifdef DEBUGMODE
	assert( minim != 0 ) ;
#endif

	/*
	 * Update the current vertex pointer.
	 */
	curr = minim ;
      }
    }
    while ( it == li.end() ) ;

    /*
     * At  this  point, a  vertex  from  the  initial front  has  been
     * reached,  and we  have finished  the computation  of  the first
     * geodesic aproximation.
     */

    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::in_the_same_mesh_element( GP* vi, GP* ve , std::list< GP* >& lv ) const
   *
   * \brief If  the initial and final  vertices of a  geodesic are the
   * same mesh element or belong  to the same mesh element, the method
   * computes  the geodesic  and returns  true. Otherwise,  it returns
   * false.
   *
   * \param vi The initial vertex of the discrete geodesic.
   * \param ve The final vertex of the discrete geodesic.
   * \param lv A reference to a list of geodesic vertices.
   *
   * \return The  logic value true  if the initial and  final vertices
   * are  in  the  same  mesh  element,  and  the  logic  value  false
   * otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::in_the_same_mesh_element(
						     GP* vi ,
						     GP* ve ,
						     std::list< GP* >& lv 
						    )
    const
  {
#ifdef DEBUGMODE
    assert( ( vi != 0 ) && ( ve != 0 ) ) ;
    assert( 
	   ( vi->on_vertex() || vi->on_edge() || vi->on_face() ) &&
	   ( ve->on_vertex() || ve->on_edge() || ve->on_face() )
	  ) ;
    assert( lv.empty() ) ;
#endif

    return vi->in_the_same_mesh_element( _mesh , ve , lv ) ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::get_coords( GP* v, double& x , double& y , double& z ) const
   *
   * \brief Returns  the three coordinates of a  given geodesic vertex.
   *
   * \param v Pointer to a geodesic vertex.
   * \param x A reference to the X coordinate.
   * \param y A reference to the Y coordinate.
   * \param z A reference to the Z coordinate.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::get_coords( 
				       GP* v , 
				       double& x , 
				       double& y , 
				       double& z 
				      ) 
    const
  {
    /*
     * The given  geodesic vertex  may coincide with  a vertex  of the
     * mesh, or  it may belong to interior  of an edge or  face of the
     * mesh.
     */
#ifdef DEBUGMODE
    assert( v != 0 ) ;
    assert( v->on_vertex() || v->on_edge() || v->on_face() ) ;
#endif

    v->get_coords( _mesh , x , y , z ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::get_restricted_adj_vertices( Vertex* v , std::vector< Vertex* >& a ) const
   *
   * \brief Finds all  mesh vertices connected to a  given mesh vertex
   * by an  edge that is  incident to a  face belonging to  the active
   * face list.
   *
   * \param v A pointer to a vertex.
   * \param a A reference to a list of vertices.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::get_restricted_adj_vertices( 
							Vertex* v , 
							std::vector< Vertex* >& a 
						       ) 
    const
  {
    Halfedge* h1 = _mesh->get_halfedge( v ) ;
    Halfedge* h2 = h1 ;

    do {
      Face* f1 = _mesh->get_face( h2 ) ;
      Face* f2 = _mesh->get_face( _mesh->get_mate( h2 ) ) ;

      /* 
       * Insert the origin vertex of  "he2" in the list of vertices if
       * and only if the edge  connecting it to vertex "v" is incident
       * with an active face.
       */
      if ( _mesh->is_active( f1 ) || _mesh->is_active( f2 ) ) {
	a.push_back( _mesh->get_org( _mesh->get_next( h2 ) ) ) ;
      }

      h2 = _mesh->get_next( _mesh->get_mate( h2 ) ) ;
    } 
    while ( h2 != h1 ) ;
    
    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::get_restricted_adj_vertices( GP* p , std::vector< Vertex* >& a ) const
   *
   * \brief Finds all mesh vertices  adjacent to the mesh element that
   * contains a given geodesic vertex.  In particular, if the geodesic
   * vertex  is a  mesh vertex,  the  method finds  the adjacent  mesh
   * vertices.   If the  geodesic vertex  is a  mesh edge,  the method
   * finds the  vertices in  the edge star.  Finally, if  the geodesic
   * vertex is a mesh face, the method finds the face vertices.
   *
   * \param p A pointer to a geodesic vertex.
   * \param a A reference to a list of vertices.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::get_restricted_adj_vertices( 
							GP* p , 
							std::vector< Vertex* >& a 
						       ) 
    const
  {
#ifdef DEBUGMODE
    assert( p != 0 ) ;
    assert( p->on_vertex() || p->on_edge() || p->on_face() ) ;
#endif

    if ( p->on_vertex() ) {
      /* Geodesic vertex is a mesh vertex. */
      GPV* gpv = dynamic_cast< GPV* >( p ) ;

#ifdef DEBUGMODE
      assert( gpv != 0 ) ;
#endif

      get_restricted_adj_vertices( gpv->get_vertex() , a ) ;
    }
    else if ( p->on_edge() ) {
      /* Geodesic vertex is in the interior of a mesh edge. */
      GPE* gpv = dynamic_cast< GPE* >( p ) ;

#ifdef DEBUGMODE
      assert( gpv != 0 ) ;
#endif

      Halfedge* h1 ;
      Halfedge* h2 ;
      _mesh->get_halfedges( gpv->get_edge() , h1 , h2 ) ;

      Face* f1 = _mesh->get_face( h1 ) ;
      Face* f2 = _mesh->get_face( h2 ) ;

#ifdef DEBUGMODE
      assert( _mesh->is_active( f1 ) || _mesh->is_active( f2 ) ) ;
#endif

      a.push_back( _mesh->get_org( h1 ) ) ;
      if ( _mesh->is_active( f1 ) ) {
	a.push_back( _mesh->get_org( _mesh->get_prev( h2 ) ) ) ;
      }
      a.push_back( _mesh->get_org( h2 ) ) ;
      if ( _mesh->is_active( f2 ) ) {
	a.push_back( _mesh->get_org( _mesh->get_prev( h1 ) ) ) ;
      }
     
    }
    else {
      /* Geodesic vertex is in the interior of a mesh face. */
      GPF* gpv = dynamic_cast< GPF* >( p ) ;

#ifdef DEBUGMODE
      assert( gpv != 0 ) ;
      assert( _mesh->is_active( gpv->get_face() ) ) ;
#endif

      Halfedge* h1 = _mesh->get_halfedge( gpv->get_face() ) ;

      a.push_back( _mesh->get_org( h1 ) ) ;
      a.push_back( _mesh->get_org( _mesh->get_next( h1 ) ) ) ;
      a.push_back( _mesh->get_org( _mesh->get_prev( h1 ) ) ) ;
    }
    
    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::is_incident_to_active_face( Vertex* v ) const
   *
   * \brief Decides whether a given vertex is incident to an active face.
   *
   * \param v A pointer to a vertex.
   *
   * \return The logic  value true if the given  vertex is incident to
   * an active face, and the logic value false otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::is_incident_to_active_face( Vertex* v ) const
  {
    Halfedge* h1 = _mesh->get_halfedge( v ) ;
    Halfedge* h2 = h1 ;

    do {
      Face* f = _mesh->get_face( h2 ) ;

      /* 
       * Insert the origin vertex of  "he2" in the list of vertices if
       * and only if the edge  connecting it to vertex "v" is incident
       * with an active face.
       */
      if ( _mesh->is_active( f ) ) {
	return true ;
      }

      h2 = _mesh->get_next( _mesh->get_mate( h2 ) ) ;
    } 
    while ( h2 != h1 ) ;
    
    return false ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::shorten_geodesic( std::list< GP* >& lv )
   *
   * \brief Shortens a  geodesic as much as possible  using the path
   * correction algorithm.
   *
   * \param lv A  list of vertices of the curve  to be shortened. This
   * list will  contain the  list of vertices  of the  shortened curve
   * after this method is finished.
   *
   * \sa compute_first_approximation.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::shorten_geodesic( std::list< GP* >& lv )
  {
    /*
     * Make sure the current path is empty.
     */
    assert( _cp_path.empty() ) ;
    
    /*
     * Create a curve  point object for each geodesic  vertex.
     */
    while ( !lv.empty() ) {
      _cp_path.push_back( CP( lv.front() ) ) ;
      lv.pop_front() ;
    }

    /*
     * Start the shortening process.
     */
    start_shortening() ;

    /*
     * Place  the  curve  points  into  the output  list  of  geodesic
     * vertices.
     */

    while ( !_cp_path.empty() ) {
      lv.push_back( _cp_path.front().get_point() ) ;
      _cp_path.pop_front() ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::shorten_geodesic( const std::list< GP* >& lv_in , std::list< GP* >& lv_out )
   *
   * \brief Shortens a  geodesic as much as possible  using the path
   * correction algorithm.
   *
   * \param lv_in A list of vertices of the curve to be shortended.
   * \param lv_out A list of verttices of the shortended curve.
   *
   * \sa compute_first_approximation.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::shorten_geodesic(
					     const std::list< GP* >& lv_in , 
					     std::list< GP* >& lv_out 
					    )
  {
    /*
     * Make sure the current path is empty.
     */
    assert( _cp_path.empty() ) ;

    /*
     * Create a curve point object for each geodesic vertex.
     */
    typename std::list< GP* >::const_iterator lit ;
    for ( lit = lv_in.begin() ; lit != lv_in.end() ; ++lit ) {
      GP* gp = *lit ;
      if ( gp->on_vertex() ) {
	GPV* gpv = dynamic_cast< GPV* >( gp ) ;
	_cp_path.push_back( CP( new GPV( *gpv ) ) ) ;
      }
      else if ( gp->on_edge() ) {
	GPE* gpe = dynamic_cast< GPE* >( gp ) ;
	_cp_path.push_back( CP( new GPE( *gpe ) ) ) ;
      }
      else {
#ifdef DEBUGMODE
	assert( gp->on_face() ) ;
#endif
	GPF* gpf = dynamic_cast< GPF* >( gp ) ;
	_cp_path.push_back( CP( new GPF( *gpf ) ) ) ;
      }
    }

    /*
     * Start the shortening process.
     */
    start_shortening() ;

    /*
     * Place  the  curve  points  into  the output  list  of  geodesic
     * vertices.
     */

    if ( !lv_out.empty() ) {
      lv_out.clear() ;
    }

    while ( !_cp_path.empty() ) {
      lv_out.push_back( _cp_path.front().get_point() ) ;
      _cp_path.pop_front() ;
    }

    return ;
  }


  /**
   * \fn double Approx_geodesics< Mesh >::compute_length() const
   *
   * \brief  Calculates the length of this geodesic.
   *
   * \return The length of this geodesic.
   */
  template< typename Mesh >
  double
  Approx_geodesics< Mesh >::compute_length() const
  {
#ifdef DEBUGMODE
    assert( _cp_path.size() > 1 ) ;
#endif 

    typename std::list< CP >::const_iterator cit = _cp_path.begin() ;
    const CP& cp1 = *cit ;

    double xi , yi , zi ;
    cp1.get_point()->get_coords( _mesh , xi , yi , zi ) ;

    ++cit ;

    double length = 0 ;
    while ( cit != _cp_path.end() ) {
      const CP& cp2 = *cit ;

      double xe, ye, ze ;
      cp2.get_point()->get_coords( _mesh , xe , ye , ze ) ;

      length += sqrt(
		     ( xe - xi ) * ( xe - xi ) +
		     ( ye - yi ) * ( ye - yi ) +
		     ( ze - zi ) * ( ze - zi )
		    ) ;
      xi = xe ;
      yi = ye ;
      zi = ze ;

      ++cit ;
    }

    return length ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::start_shortening()
   *
   * \brief Carry on the actual geodesic shortening process.
   *
   * \sa compute_first_approximation, shorten_geodesic
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::start_shortening()
  {
    /*
     * Compute the length of this geodesic.
     */
    double old_length = compute_length() ;

    /*
     * Initialize iteration counter.
     */
    _num_iter = 0 ;

    /*
     * If  the geodesics  has at  least  three vertices,  it may  need
     * correction.  Otherwise,  there is  no need for  correction (the
     * geodesic is the true one).
     */
    if ( _cp_path.size() > 2 ) {
      /*
       * The  entire  geodesic  is  sequentially  traversed  from  its
       * initial to  its final vertex.  For each  interior vertex, its
       * "error"  is   computed  and  tested   against  a  pre-defined
       * threshold.  If  the error is no smaller  than the pre-defined
       * threshold,  the  geodesic  is  locally shortened  around  the
       * vertex.   The loop  ends as  soon as  the geodesic  cannot be
       * further  shortened, the  maximum  error is  smaller than  the
       * pre-defined threshold, or the maximum number of iterations is
       * reached.
       */

      double length_ratio = 1 ;
      do {
	/*
	 * Define  an iterator  for the  first three  vertices  of the
	 * current geodesic.
	 */
	Path_iterator it1 , it2 , it3 ;

	/*
	 * Let the iterator  variables it1, it2, and it3  point to the
	 * first, second, and third  vertices of the current geodesic,
	 * respectively.
	 */
	it1 = it2 = _cp_path.begin() ;

	++it2 ;
	
	it3 = it2 ;
	++it3 ;

	/*
	 * We now loop over the entire geodesic until we reach its end
	 * point. The method  eliminate_vertices() makes sure there is
	 * no  redundant geodesic vertices  from the  geodesic initial
	 * point  to the geodesic  vertex pointed  by iterator  it3. A
	 * vertex is  considered redundant if it shares  the same mesh
	 * face with  the previous and next  vertices.  By eliminating
	 * redundant   vertices,    the   method   may    change   the
	 * iterators. Finally, the method returns the logic value true
	 * if  it2  is  not   the  last  geodesic  vertex,  and  false
	 * otherwise.
	 */
	while ( eliminate_vertices( it1 , it2 , it3 ) ) {
	  /*
	   * The vertices  pointed by it1,  it2, and it3 do  not share
	   * the same mesh  face. If they do, there  is a problem with
	   * the code.
	   */

	  /*
	   * Compute the  error associated with the  vertex pointed by
	   * it2.
	   */
	  pre_processing( it1 , it2 , it3 ) ;

	  /*
	   * Get the geodesic vertex pointed by it2.
	   */
	  CP& cp = *it2 ;

	  /*
	   * If the  error associated  with current vertex  is smaller
	   * than the pre-defined threshold, we can ignore the current
	   * vertex.  Otherwise, we try to shorten the geodesic.
	   */

	  /*
	   * Flag  to  indicate how  to  set  up  the geodesic  vertex
	   * iterators for the next loop iteration.
	   */
	  bool move_ahead = false ;

	  if ( cp.get_error() >= _error ) {
	    /*
	     * If vertex cp is in the  interior of a mesh edge, it may
	     * be moved  along the edge or  snapped to one  of the two
	     * edge  end  points.    Otherwise,  the  geodesic  vertex
	     * coincides with a mesh vertex,  and it will be moved (as
	     * it  is  in the  heap).   Whenever  the vertex  position
	     * changes, the geodesic is shortened.
	     */
	    bool has_shortened = true ;
	    if ( cp.on_edge() ) {
	      has_shortened = pos_on_e( it1 , it2 , it3 ) ;

	      /*
	       * If the  geodesic vertex  is now too  close to  a mesh
	       * vertex,  we must verify  the possibility  of snapping
	       * the former to the latter.  Note that the test to find
	       * out whether  "cp" points to a vertex  in the interior
	       * of a mesh edge is  necessary.  The reason is that the
	       * method "pos_on_e()" can make  a point in the interior
	       * of a mesh edge into a mesh vertex.
	       */
	      bool has_snapped = false ;
	      if ( cp.on_edge() ) {
		has_snapped = snap_to_vertex( it1 , it2 , it3 ) ;
	      }
	      
	      /*
	       * It might be  the case that the geodesic  has not been
	       * moved before the test  for snapping, but the point is
	       * then moved  by snapping.  Consequently,  the geodesic
	       * is  locally  changed  and  we must  update  the  heap
	       * information  of the points  affected by  the snapping
	       * (including the snapped point itself).
	       */
	      if ( !has_shortened && has_snapped ) {
		has_shortened = true ;
	      }
	    }
	    else {
#ifdef DEBUGMODE
	      assert( cp.on_vertex() ) ;
#endif
	      pos_on_v( it1 , it2 , it3 ) ;
	    }

	    /*
	     * If the geodesic has  been shortened by the displacement
	     * of  the vertex pointed  by iterator  it2, then  we move
	     * iterator it1 to one geodesic vertex ahead.
	     */
	    if ( has_shortened ) {
	      /*
	       * It is  possible that the vertices pointed  by it1 and
	       * its two immediate following  vertices are in the same
	       * mesh face.  If so,  the second vertex must be removed
	       * before  iterator  it1  is  set  to  the  next  vertex
	       * position.
	       */
	      Path_iterator pit1 , pit2 , pit3 ;
	      pit1 = it1 ;
	      pit2 = pit1 ;
	      if ( it1 != _cp_path.begin() ) {
		--pit1 ;
	      }
	      else {
		++pit2 ;
	      }

#ifdef DEBUGMODE
	      assert( pit2 != _cp_path.end() ) ;
#endif

	      pit3 = pit2 ;
	      ++pit3 ;

	      if ( it1 != _cp_path.begin() ) {
		if ( !can_be_removed( *pit1 , *pit2 , *pit3 ) ) {
		  it1 = pit3 ;
		}
	      }
	      else if ( !can_be_removed( *pit1 , *pit2 , *pit3 ) ) {
		it1 = pit2 ;
	      }

	      move_ahead = false ;
	    }
	    else {
	      move_ahead = true ;
	    }
	  }
	  else {
	    move_ahead = true ;
	  }

	  /*
	   * set up iterators for the next loop iterations.
	   */

	  /*
	   * Move forward it1 one vertex position.
	   */
	  if ( move_ahead ) {
	    ++it1 ;
	  }

	  it2 = it1 ;
	  ++it2 ;
	  
#ifdef DEBUGMODE
	  assert( it2 != _cp_path.end() ) ;
#endif
	  it3 = it2 ;
	  ++it3 ;

	  /*
	   * Increment the iteration counter.
	   */
	  ++_num_iter ;
	}
	
	/*
	 * The current geodesic has  been traversed once. Now, we must
	 * check  whether there  is a  need for  trying to  shorten it
	 * again.
	 */
	double new_length = compute_length() ;
	length_ratio = new_length / old_length ;
	old_length = new_length ;
      }
      while (
	     ( _cp_path.size() > 2 ) &&
	     ( length_ratio < _SHORTENING_RATIO ) &&
	     ( _num_iter < _max_iter )
	    ) ;
    }

#ifdef DEBUGMODE
    /*
     * Check geodesic consistency
     */
    assert( is_consistent() ) ;
#endif
      
    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::eliminate_vertices( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief  Eliminates  all geodesic  vertices  that  follow a  given
   * vertex and share one mesh face with its predecessor and successor
   * vertices.
   *
   * \param prev The predecessor of the given vertex in the geodesic.
   * \param node The given vertex.
   * \param next The successor of the given vertex in the geodesic.
   *
   * \return The logic  value true if the end of  the geodesic has not
   * been reached during the  elimination process, and the logic value
   * false otherwise.
   */
  template< typename Mesh >
  bool 
  Approx_geodesics< Mesh >::eliminate_vertices( 
					       Path_iterator& prev ,
					       Path_iterator& node , 
					       Path_iterator& next 
					      ) 
  {
    /*
     * Make sure the iterators point to three consecutive vertices.
     */
#ifdef DEBUGMODE
    assert( prev != _cp_path.end() ) ;
    assert( node != _cp_path.begin() ) ;
    Path_iterator temp_it1 = node ;
    --temp_it1 ;
    assert( ( *temp_it1 ).get_point() == ( *prev ).get_point() ) ;
#endif

    /*
     * If the third vertex do not exist, there is nothing to be done.
     */
    if ( next == _cp_path.end() ) {
      return false ;
    }

    /*
     * Make sure the iterators point to three consecutive vertices.
     */
#ifdef DEBUGMODE
    assert( node != _cp_path.end() ) ;
    Path_iterator temp_it3 = node ;
    ++temp_it3 ;
    assert( ( *temp_it3 ).get_point() == ( *next ).get_point() ) ;
#endif

    /*
     * While  the   geodesic  vertex  pointed  by   iterator  node  is
     * redundant, remove the vertex  and move back iterator "prev" one
     * position.
     */
    while ( can_be_removed( *prev , *node , *next ) ) {
      /*
       * The vertex  pointed by "node" is  in the same  face as either
       * the previous  or the  next vertex or  both. In any  case, the
       * vertex  must be removed  from the  geodesic (the  geodesic is
       * shortened!).
       */
      CP& cp_node = *node ;
      _cp_path.erase( node ) ;
      
      /*
       * Release the memory associated with the given geodesic vertex.
       */
      if ( cp_node.get_point() != 0 ) {
	delete cp_node.get_point() ;
      }

      /*
       * Move back iterator "prev" one position (if possible).
       */
      if ( prev != _cp_path.begin() ) {
	--prev ;
      }

      /*
       * Set up a new position for iterator "node"
       */
      node = prev ;
      ++node ;

#ifdef DEBUGMODE
      assert( node != _cp_path.end() ) ;
#endif

      /*
       * Set up a new position for iterator "next" (if possible)
       */
      next = node ;
      ++next ;

      /*
       * If the end of the geodesic has been reached, there is nothing
       * else  to be  done.  Otherwise,  we keep  checking  for vertex
       * redundancy.
       */
      if ( next == _cp_path.end() ) {
	return false ;
      }
    }

    return true ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::pre_processing( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief  Calculates the  error  associated with  a given  geodesic
   * vertex.
   *
   * \param prev The predecessor of the given vertex in the geodesic.
   * \param node The given vertex.
   * \param next The successor of the given vertex in the geodesic.
   *
   * \return The logic value true if the current vertex has not been
   * removed from the geodesic, and the logic value false otherwise.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::pre_processing( 
					   Path_iterator& prev ,
					   Path_iterator& node , 
					   Path_iterator& next 
					  ) 
  {
    CP& cp = *node ;

    /*
     * Make  sure the  vertex is  a  mesh vertex  or a  vertex in  the
     * interior of an  edge. It can never be a  vertex in the interior
     * of a face.
     */
#ifdef DEBUGMODE
    assert( cp.on_vertex() || cp.on_edge() ) ;
#endif

#ifdef DEBUGMODE
    assert( !can_be_removed( *prev , *node , *next ) ) ;
#endif

    /*
     * Compute the error associated with the given geodesic vertex.
     */
    if ( cp.on_edge() ) {
      compute_error_on_e( prev , node , next ) ;
    }
    else {
      compute_error_on_v( prev , node , next ) ;
    }

    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::pos_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief Shortens  the geodesic by straightening it  in the segment
   * defined  by  three  given  vertices, which  are  named  previous,
   * current, and next. The current  vertex belongs to the interior of
   * an edge.
   *
   * \param prev The previous vertex.
   * \param node The current vertex.
   * \param next The next vertex.
   *
   * \return The logic  value true if the geodesic  has been shortened
   * by  moving   the  current  vertex   along  the  interior   of  an
   * edge. Otherwise, the logic value false is returned.
   *
   * \sa pre_on_v, pre_on_e, pre_processing, pos_on_v 
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::pos_on_e( 
				     Path_iterator& prev ,
				     Path_iterator& node ,
				     Path_iterator& next 
				    ) 
  {
#ifdef DEBUGMODE
    assert( (*node).on_edge() ) ;
#endif

    /*
     * Previous vertex.
     */
    CP& cp_prev = *prev ;

    /*
     * Current vertex 
     */
    CP& cp_node = *node ;

    /*
     * Next vertex
     */
    CP& cp_next = *next ;
    
    return compute_intersection_on_e( cp_prev , cp_node , cp_next ) ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::pos_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief Shortens  the geodesic by straightening it  in the segment
   * defined  by  three  given  vertices, which  are  named  previous,
   * current, and next. The current  vertex belongs to the interior of
   * an edge.
   *
   * \param prev The previous vertex.
   * \param node The current vertex.
   * \param next The next vertex.
   *
   * \sa pos_on_e 
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::pos_on_v( 
				     Path_iterator& prev ,
				     Path_iterator& node ,
				     Path_iterator& next 
				    ) 
  {
#ifdef DEBUGMODE
    assert( ( *node ).on_vertex() ) ;
#endif

    /*
     * Previous vertex.
     */
    CP& cp_prev = *prev ;

    /*
     * Current vertex 
     */
    CP& cp_node = *node ;

    /*
     * Next vertex
     */
    CP& cp_next = *next ;

    /*
     * Get the cached information  to compute the intersection between
     * the line segment  from the previous to the  next vertex and the
     * edges  of the  unfolded triangles  of the  star of  the current
     * vertex.
     */
    Halfedge* h_prev ;
    Halfedge* h_next ;
    get_limits( h_prev , h_next ) ;

    SIDE side ;
    side = cp_node.get_side() ;

#ifdef DEBUGMODE
    assert(
        ( side != CP::NO_SIDE ) &&
        ( side != CP::NON_INITIALIZED )
    ) ;
#endif

    /*
     * Find the intersection points.
     */
    Path new_path ;
    if ( side == CP::LEFT_SIDE ) {
      compute_intersection_on_v( 
				cp_next , 
				cp_node , 
				cp_prev , 
				new_path 
			       ) ;
      new_path.reverse() ;
    } 
    else {
      compute_intersection_on_v(
				cp_prev , 
				cp_node , 
				cp_next , 
				new_path 
			       ) ;
    }

    /*
     * The current vertex must always be eliminated from the geodesic,
     * as the geodesic  has been shortened to pass on  one side of the
     * current vertex.
     */
    _cp_path.erase( node ) ;

    /*
     * Release the memory associated with the removed vertex.
     */
    if ( cp_node.get_point() != 0 ) {
      delete cp_node.get_point() ;
    }

    /*
     * Insert into  the current geodesic  its new vertices  (which are
     * the intersection points between a line segment and the unfolded
     * triangle edges).
     */
    _cp_path.insert( next , new_path.begin() , new_path.end() ) ;

    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::snap_to_vertex( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief Snaps a given geodesic  vertex, which lies in the interior
   * of a mesh  edge, to a mesh vertex if the  geodesic vertex and its
   * immediate neighbors are too close  to the mesh vertex and if they
   * belong to the concave local side of the geodesic.
   *
   * \param prev The previous vertex.
   * \param node The current vertex.
   * \param next The next vertex.
   *
   * \return The  logic value  true if the  given geodesic  vertex has
   * been  snapped  to  a  mesh  vertex, and  the  logic  value  false
   * otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::snap_to_vertex(
					   Path_iterator& prev ,
					   Path_iterator& node , 
					   Path_iterator& next 
					  ) 
  {
    /*
     * Get a reference to the given node.
     */
    CP& cp_node = *node ;

#ifdef DEBUGMODE
    /*
     * Given vertex must be a geodesic vertex lying in the interior of
     * a mesh edge.
     */
    assert( cp_node.on_edge() ) ;

#endif

    /*
     * Is the given geodesic vertex too close to one of the end points
     * of the mesh edge it lies  on? If so, we must check the distance
     * from  its neighbors  (in the  geodesic)  to the  same edge  end
     * point.
     */
    if ( 
	(       cp_node.get_t()   < _TOL_E_NEAR_V ) || 
	( ( 1 - cp_node.get_t() ) < _TOL_E_NEAR_V ) 
       ) {
      /*
       * Get the mesh edge the given geodesic vertex lies on.
       */
      Edge* e = cp_node.get_edge() ;

      /*
       * Get the two half-edges of the mesh edge.
       */
      Halfedge* h_temp[ 2 ] ;
      _mesh->get_halfedges( e , h_temp[0] , h_temp[1] ) ;

      /*
       * Find out the closest edge end point (a mesh vertex).
       */
      Halfedge* h ;
      if ( cp_node.get_t() < 0.5 ) {
        h = h_temp[ 0 ] ;
      }
      else {
        h = h_temp[ 1 ] ;
      }

      /*
       * Starting from  the previous vertex, we  traverse the geodesic
       * backward to find a geodesic vertex on a mesh vertex or in the
       * interior  of a  mesh edge  that does  not contain  the origin
       * vertex of h.
       */
      bool done = false ;
      Path_iterator it1 = prev ;
      while ( !done ) {
	if ( ( it1 == _cp_path.begin() ) || ( *it1 ).on_vertex() ) {
	  done = true ; 
	}
	else {
#ifdef DEBUGMODE
	  assert( ( *it1 ).on_edge() ) ;
#endif
	  _mesh->get_halfedges( 
			       ( *it1 ).get_edge() , 
			       h_temp[ 0 ] , 
			       h_temp[ 1 ] 
			      ) ;

	  if ( 
	      ( _mesh->get_org( h_temp[ 0 ] ) != _mesh->get_org( h ) ) &&
	      ( _mesh->get_org( h_temp[ 1 ] ) != _mesh->get_org( h ) )
	     ) {
	    done = true ;
	  }
	  else {
	    --it1 ;
	  }
	}
      }


      /*
       * Starting  from  the next  vertex,  we  traverse the  geodesic
       * forward to find a geodesic vertex  on a mesh vertex or in the
       * interior  of a  mesh edge  that does  not contain  the origin
       * vertex of h.
       */
      done = false ;
      Path_iterator it3 = next ;
      while ( !done ) {
	Path_iterator it4 = it3 ;
	++it4 ;
	if ( ( it4 == _cp_path.end() ) || ( *it3 ).on_vertex() ) {
	  done = true ; 
	}
	else {
#ifdef DEBUGMODE
	  assert( ( *it3 ).on_edge() ) ;
#endif
	  _mesh->get_halfedges( 
			       ( *it3 ).get_edge() , 
			       h_temp[ 0 ] , 
			       h_temp[ 1 ] 
			      ) ;

	  if ( 
	      ( _mesh->get_org( h_temp[ 0 ] ) != _mesh->get_org( h ) ) &&
	      ( _mesh->get_org( h_temp[ 1 ] ) != _mesh->get_org( h ) )
	     ) {
	    done = true ;
	  }
	  else {
	    ++it3 ;
	  }
	}
      }

      /*
       * We now must  find out whether the geodesic  vertex pointed by
       * iterator node is on the  concave or convex side of the sector
       * defined by the origin vertex of half-edge h, and the vertices
       * pointed by iterators it1 and it3. If it is in the concave, we
       * can snap the geodesic vertex  pointed by iterator node to the
       * origin vertex of half-edge  h. Otherwise, we leave the vertex
       * where it is.
       */
      CP& cp1 = *it1 ;
      CP& cp3 = *it3 ;

      /*
       * Creates the  possible new geodesic  vertex (one end  point of
       * the edge).
       */
      CP cp_new( new GPV( _mesh->get_org( h ) ) ) ;

      /*
       * Compute the left  and right angles at the  vertex cp_new with
       * respect to the two  geodesic segments defined by cp1, cp_new,
       * and cp3.
       */
      Halfedge* h_prev ;
      Halfedge* h_next ;
      search_halfedges_on_v( 
			    cp1 , 
			    cp_new , 
			    cp3 , 
			    h_prev , 
			    h_next 
			   ) ;

      double rangle , langle ;
      compute_angles_on_v( 
			  h_prev , 
			  h_next , 
			  cp1 , 
			  cp_new , 
			  cp3 , 
			  langle ,
			  rangle 
			 ) ;

      /*
       * Find out the geodesic side the given geodesic vertex is.
       */
      Path_iterator it2 = it1 ;
      ++it2 ;
      CP& cp2 = *it2 ;

      Edge* e2 = cp2.get_edge() ;
      Halfedge* h2[ 2 ] ;
      _mesh->get_halfedges( e2 , h2[ 0 ] , h2[ 1 ] ) ;

      bool snap = false ;
      if ( 
	  ( h2[ 0 ] == _mesh->get_prev( h_prev ) ) || 
	  ( h2[ 1 ] == _mesh->get_prev( h_prev ) ) 
	 ) {

	/*
	 * The  geodesic vertex  pointed by  iterator node  is  on the
	 * right  side of  the geodesic.   Now, if  the right  side is
	 * concave,  then we can  definitely snap  this vertex  to the
	 * origin of h.
	 */
        if ( langle < _MYPI ) {
          snap = true ;
        }
      } 
      else {
	/*
	 * The geodesic vertex pointed by iterator node is on the left
	 * side of  the geodesic.  Now,  if the left side  is concave,
	 * then we can definitely snap this vertex to the origin of h.
	 */
        if ( rangle < _MYPI ) {
          snap = true ;
        }
      }

      /*
       * If the snap flag has been set to true, the snapping operation
       * takes place.  This amount  to the elimination of all vertices
       * in  between the vertices  pointed by  iterators it1  and it4,
       * including  the  vertex  pointed  by iterator  node  (the  one
       * snapped to a mesh vertex).
       */
      if ( snap ) {
        for ( Path_iterator it = it2 ; it != it3 ; ++it ) {
          CP &cp = *it ;
	  /*
	   * Release the memory associated with the actual vertex.
	   */
	  if ( cp.get_point() != 0 ) {
	    delete cp.get_point() ;
	  }
        }

        _cp_path.erase( it2 , it3 ) ;

        _cp_path.insert( it3 , cp_new ) ;

        prev = it1 ;
        next = it3 ;
	
	/*
	 * The geodesic has been modified!
	 */
        return true ;
      }
      else {
	/*
	 * The snap operation  did not take place, so  we must release
	 * the memory allocated for the new vertex.
	 */
#ifdef DEBUGMODE
	assert( cp_new.get_point() != 0 ) ;
#endif
	delete cp_new.get_point() ;
      }
    }

    /*
     * Given vertex is not too close to a mesh vertex.
     */
    return false ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::can_be_removed( const CP& cp_prev , const CP& cp_node , const CP& cp_next ) const
   *
   * \brief  Decides whether a  given geodesic  vertex can  be removed
   * from the  current geodesic path. The decision  takes into account
   * the  predecessor  and successor  of  the  vertex  in the  current
   * geodesic.
   *
   * \param cp_prev Reference to the  previous vertex.
   * \param cp_node Reference to a vertex.
   * \param cp_next Reference to the next vertex.
   *
   * \return The  logic value true if  a given geodesic  vertex can be
   * removed, and the logic value false otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::can_be_removed(
					   const CP& cp_prev , 
					   const CP& cp_node ,
					   const CP& cp_next
					  ) 
    const 
  {
#ifdef DEBUGMODE
    /*
     * The vertex to be tested against removal cannot belong to a mesh
     * face.
     */
    assert( cp_node.on_vertex() || cp_node.on_edge() ) ;
#endif

    /*
     * Find the  common face(s) of  the previous and  current geodesic
     * vertices. There  must be  at least one  and at most  two common
     * faces.
     */
    Face* f1 = 0 ;
    Face* f2 = 0 ;
    if ( cp_prev.on_vertex() ) {
      /*
       * Previous geodesic vertex is a mesh vertex.
       */
      Vertex* v_prev = cp_prev.get_vertex() ;
      if ( cp_node.on_vertex() ) {
	/*
	 * Current geodesic vertex is also a mesh vertex.
	 */
	Vertex* v_node = cp_node.get_vertex() ;

	/*
	 * Find the common faces of two mesh vertices.
	 */
	find_common_faces( v_prev , v_node , f1 , f2 ) ;
      }
      else {
	/*
	 * Current geodesic vertex is in the interior of a mesh edge.
	 */
	Edge* e_node = cp_node.get_edge() ;

	/*
	 * Find the common faces of a mesh vertex and a mesh edge.
	 */
	find_common_faces( v_prev , e_node , f1 , f2 ) ;
      }
    }
    else if ( cp_prev.on_edge() ) {
      /*
       * Previous geodesic vertex is a mesh vertex.
       */
      Edge* e_prev = cp_prev.get_edge() ;
      if ( cp_node.on_vertex() ) {
	/*
	 * Current geodesic vertex is a mesh vertex.
	 */
	Vertex* v_node = cp_node.get_vertex() ;

	/*
	 * Find the common faces of a mesh vertex and a mesh edge.
	 */
	find_common_faces( v_node , e_prev , f1 , f2 ) ;
      }
      else {
	/*
	 * Current geodesic vertex is also a mesh edge.
	 */
	Edge* e_node = cp_node.get_edge() ;

	/*
	 * Find the common faces of two mesh edges.
	 */
	find_common_faces( e_prev , e_node , f1 , f2 ) ;
      }
    }
    else {
      /*
       * Previous geodesic vertex is in the interior of a mesh face.
       */

#ifdef DEBUGMODE
      assert( cp_prev.on_face() ) ;
#endif

      Face* f_prev = cp_prev.get_face() ;
      f2 = 0 ;
      if ( cp_node.on_vertex() ) {
	/*
	 * Current geodesic vertex is a mesh vertex.
	 */
	Vertex* v_node = cp_node.get_vertex() ;

	/*
	 * Find out if the mesh vertex belongs to face f_prev.
	 */
	if ( is_in_face ( v_node , f_prev ) ) {
	  f1 = f_prev ;
	}
      }
      else {
	/*
	 * Current geodesic vertex is a mesh edge.
	 */
	Edge* e_node = cp_node.get_edge() ;

	/*
	 * Find out if mesh edge belongs to face f_prev.
	 */
	if ( is_in_face ( e_node , f_prev ) ) {
	  f1 = f_prev ;
	}
      }
    }

 #ifdef DEBUGMODE
    assert( f1 != 0 ) ;
 #endif

    /*
     * Find out if one of f1 and f2 contains the next geodesic vertex.
     */
    if ( cp_next.on_vertex() ) {
      /*
       * Next geodesic vertex is a mesh vertex.
       */
      Vertex* v_next = cp_next.get_vertex() ;

      /*
       * Find out if the geodesic vertex belongs to one of f1 or f2.
       */
      if ( is_in_face ( v_next , f1 ) ) {
	return true ;
      }
      else if ( f2 != 0 ) {
	return is_in_face ( v_next , f2 ) ;
      }

      return false ;
    }
    else if ( cp_next.on_edge() ) {
      /*
       * Next geodesic vertex is a mesh edge.
       */
      Edge* e_next = cp_next.get_edge() ;

      /*
       * Find out if the mesh edge belongs to one of f1 or f2.
       */
      if ( is_in_face ( e_next , f1 ) ) {
	return true ;
      }
      else if ( f2 != 0 ) {
	return is_in_face ( e_next , f2 ) ;
      }

      return false ;
    }
    else {
      /*
       * Next geodesic vertex is a mesh face.
       */

#ifdef DEBUGMODE
      assert( cp_next.on_face() ) ;
#endif

      Face* f_next = cp_next.get_face() ;

#ifdef DEBUGMODE
      assert( f_next != 0 ) ;
#endif

      /*
       * Find out if the mesh face is either f1 or f2.
       */
      if ( f1 == f_next ) {
	return true ;
      }
      else if ( f2 != 0 ) {
	return f_next == f2 ;
      }

      return false ;
    }

    return false ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::update_error_info_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief  Updates the  information needed  for computing  the error
   * associated  with a given  geodesic vertex,  which belongs  to the
   * interior of a mesh edge.
   *
   * \param prev The predecessor of the given vertex in the geodesic.
   * \param node The given vertex.
   * \param next The successor of the given vertex in the geodesic.
   *
   * \sa pre_processing
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::update_error_info_on_e(
						   Path_iterator& prev ,
						   Path_iterator& node ,
						   Path_iterator& next
						  )
  {
#ifdef DEBUGMODE
    assert( (*node).on_edge() ) ;
#endif

    /*
     * Previous vertex.
     */
    CP& cp_prev = *prev ;

    /*
     * Current vertex - the one we wish to compute information about.
     */
    CP& cp_node = *node ;

    /*
     * Next vertex
     */
    CP& cp_next = *next ;

    /*
     * To  compute the  2D coordinates  of  a geodesic  vertex in  the
     * interior of an edge, we  must find the two triangles that share
     * the  mesh  edge the  current  vertex  of  the geodesic  belongs
     * to.
     */
    Halfedge* h_prev ;
    Halfedge* h_next ;
    search_halfedges_on_e(
			  cp_prev ,
			  cp_node ,
			  cp_next ,
			  h_prev ,
			  h_next 
			 ) ;

#ifdef DEBUGMODE
    assert( h_prev != 0 ) ;
    assert( h_next != 0 ) ;
#endif 

#ifdef DEBUGMODE
    Halfedge* h_temp = _mesh->get_prev( h_prev ) ;
    assert( _mesh->get_edge( h_temp ) == cp_node.get_edge() ) ;
    assert( _mesh->get_edge( h_next ) == cp_node.get_edge() ) ;
#endif 

    set_limits( h_prev , h_next ) ;

    /*
     * Unfold the two triangles delimited by h_prev and h_next.
     */
    unfold_triangle_strip( h_prev , h_next ) ;

#ifdef DEBUGMODE
    assert( _pts.size() == 6 ) ;
#endif 

    /*
     * Compute the 2D oordinates of the unfolded triangles.
     */
    compute_2d_coords_on_e( cp_prev , cp_node , cp_next ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_error_on_e( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief  Computes  the  error  associated with  a  given  geodesic
   * vertex that belongs to the interior of a mesh edge.
   *
   * \param prev The predecessor of the given vertex in the geodesic.
   * \param node The given vertex.
   * \param next The successor of the given vertex in the geodesic.
   *
   * \sa pre_processing, update_error_info_on_e
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::compute_error_on_e(
					       Path_iterator& prev ,
					       Path_iterator& node ,
					       Path_iterator& next
					      )
  {
    /*
     * Make sure the current vertex error information is updated.
     */ 
    update_error_info_on_e( prev , node , next ) ;

    /*
     * Current vertex - the one we wish to compute information about.
     */
    CP& cp_node = *node ;

    /*
     * Get the 2D coordinates of the unfolded triangles related to the
     * current vertex.
     */
    double prev_x , prev_y , node_x , node_y , next_x , next_y ;

    get_2d_coords(
		  prev_x ,
		  prev_y ,
		  node_x ,
		  node_y ,
		  next_x ,
		  next_y
		 ) ;

    /*
     * Compute the angle defined at the given vertex.
     */
    double ang = Geometric::get_angle( 
				      node_x ,
				      node_y ,
				      0 ,
				      prev_x ,
				      prev_y ,
				      0 ,
				      next_x ,
				      next_y ,
				      0 
				     ) ;

    /*
     * Update the error value associated with the given vertex.
     */
    cp_node.set_error( fabs( _MYPI - ang ) ) ;
    
    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::update_error_info_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next , double& langle , double& rangle )
   *
   * \brief  Updates the  information needed  for computing  the error
   * associated with  a geodesic given vertex, which  coincides with a
   * mesh vertex.
   *
   * \param prev The predecessor of the given vertex in the geodesic.
   * \param node The given vertex.
   * \param next The successor of the given vertex in the geodesic.
   * \param langle The angle defined  by the left side of the geodesic
   * at the given vertex.
   * \param rangle The angle defined by the right side of the geodesic
   * at the given vertex.
   *
   * \sa pre_processing
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::update_error_info_on_v(
					    Path_iterator& prev ,
					    Path_iterator& node ,
					    Path_iterator& next ,
					    double& langle ,
					    double& rangle
					   )
  {
#ifdef DEBUGMODE
    assert( (*node).on_vertex() ) ;
#endif

    /*
     * Previous vertex.
     */
    CP& cp_prev = *prev ;

    /*
     * Current vertex - the one we wish to compute information about.
     */
    CP& cp_node = *node ;

    /*
     * Next vertex
     */
    CP& cp_next = *next ;

    /*
     * Find the  two limiting half-edges of the  current vertex. These
     * half-edges  identify the  initial  and final  triangles of  the
     * triangle  strip  that contains  the  predecessor and  successor
     * vertices of the current geodesic vertex.
     */
    Halfedge* h_prev ;
    Halfedge* h_next ;
    search_halfedges_on_v(
			  cp_prev ,
			  cp_node ,
			  cp_next ,
			  h_prev ,
			  h_next 
			 ) ;

#ifdef DEBUGMODE
    assert( h_prev != 0 ) ;
    assert( h_next != 0 ) ;
    assert( _mesh->get_org( h_prev ) == cp_node.get_vertex() ) ;
    assert( _mesh->get_org( h_next ) == cp_node.get_vertex() ) ;
#endif 

    /*
     * Compute the left and right angles at the current vertex. 
     */
    compute_angles_on_v( 
			h_prev , 
			h_next , 
			cp_prev , 
			cp_node , 
			cp_next ,
			langle , 
			rangle 
		       ) ;

    /*
     * If  the left angle  is smaller  than the  right angle  and also
     * smaller than  PI, we must  reorient the limiting  half-edges of
     * the vertex.
     */
    SIDE side = CP::NO_SIDE ;
    if ( rangle <= langle ) {
      if ( rangle < _MYPI ) {
	side = CP::RIGHT_SIDE ;
      }
    }
    else {
      if ( langle < _MYPI ) {
	/*
	 * When the  left angle  is the smallest,  we must  update the
	 * half-edges that delimit the  first and last triangle cut by
	 * the geodesic.
	 */
	if ( cp_prev.on_edge() ) {
	  if ( _mesh->get_edge( h_prev ) == cp_prev.get_edge() ) {
	    h_prev = _mesh->get_next( _mesh->get_mate( h_prev ) ) ;
	  }
	} 
	else if ( cp_prev.on_vertex() ) {
	  h_prev = _mesh->get_next( _mesh->get_mate( h_prev ) ) ;
	}

	if ( cp_next.on_edge() ) {
	  if ( _mesh->get_edge( _mesh->get_prev( h_next ) ) == 
	       cp_next.get_edge() ) {
	    h_next = _mesh->get_mate( _mesh->get_prev( h_next ) ) ;
	  }
	} 
	else if ( cp_next.on_vertex() ) {
	  h_next = _mesh->get_mate( _mesh->get_prev( h_next ) ) ;
	}

	/*
	 * Note   that  nothing   has  been   done  if   the  previous
	 * (resp. next) vertex is in the  interior of a face or in the
	 * interior  of  an  edge  opposite to  the  current  geodesic
	 * vertex.
	 */ 
	Halfedge* h_temp = h_prev ;
	h_prev = h_next ;
	h_next = h_temp ;

	side = CP::LEFT_SIDE ;
      }
    }

    /*
     * Assign side information to the current vertex.
     */
    cp_node.set_side( side ) ;   

    /*
     * Set up the limiting half-edges.
     */
    set_limits( h_prev , h_next ) ;

    /*
     * Compute the unfolded triangle coordinates.
     */
    if ( cp_node.get_side() != CP::NO_SIDE ) {
      unfold_triangle_strip( 
			    h_prev ,
			    h_next
			   ) ;

#ifdef DEBUGMODE
      assert(
	     (   _pts.size()       >= 6 ) &&
	     ( ( _pts.size() % 2 ) == 0 )
	    ) ;
#endif

      if ( cp_node.get_side() == CP::RIGHT_SIDE ) {
	/*
	 * Compute the planar coordinates of the current vertex in the
	 * unfolded strip of triangles. This is only need if the right
	 * and the left angles at the current vertex are distinct, and
	 * the smallest of  the two is not larger  than PI. Otherwise,
	 * the geodesic cannot be shortened in the neighborhood of the
	 * current vertex.
	 */
	compute_2d_coords_on_v( cp_prev , cp_node , cp_next ) ;
      }
      else {
	compute_2d_coords_on_v( cp_next , cp_node , cp_prev ) ;
      }
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_error_on_v( Path_iterator& prev , Path_iterator& node , Path_iterator& next )
   *
   * \brief Computes the error associated with a given geodesic vertex
   * that coincides with a mesh vertex.
   *
   * \param prev The predecessor of the given vertex in the geodesic.
   * \param node The given vertex.
   * \param next The successor of the given vertex in the geodesic.
   *
   * \sa pre_processing, update_error_info_on_v
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::
  compute_error_on_v(
		     Path_iterator& prev ,
		     Path_iterator& node ,
		     Path_iterator& next
		    )
  {
    /*
     * Make sure the current vertex error information is updated.
     */ 
    double langle ;
    double rangle ;
    update_error_info_on_v( 
			   prev , 
			   node , 
			   next , 
			   langle , 
			   rangle 
			  ) ;

    /*
     * Current vertex - the one we wish to compute information about.
     */
    CP& cp_node = *node ;

    /*
     * Get the limiting half-edges of the current vertex.
     */
    Halfedge* h_prev ;
    Halfedge* h_next ;
    get_limits( h_prev , h_next ) ;

    /*
     * If  both  angles  are  larger  than  PI, there  is  no  way  of
     * shortening the geodesic by moving  the given vertex. So, we set
     * the  error to zero.  Otherwise, an  error value  is set  to the
     * vertex.
     */
    if ( cp_node.get_side() == CP::NO_SIDE ) {
      cp_node.set_error( 0 ) ;
    }
    else {
      /*
       * The error is defined as  half the absolute difference of both
       * angles.
       */
      cp_node.set_error( 0.5 * fabs( rangle - langle ) ) ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::search_halfedges_on_e( const CP& cp_prev, const CP& cp_node, const CP& cp_next, Halfedge*& h_prev, Halfedge*& h_next ) const
   *
   * \brief Determines  the two triangles  that share a mesh  edge are
   * cut through by a geodesic that passes through the interior of the
   * edge.
   *
   * \param cp_prev The previous vertex.
   * \param  cp_node The  geodesic  vertex that  belongs  to the  edge
   * interior.
   * \param cp_next The next vertex.
   * \param h_prev A  reference to a pointer to a  half-edge of one of
   * the two triangles. The origin  vertex of the half-edge is one end
   * vertex of the  mesh edge containing the given  geodesic vertex in
   * its interior.
   * \param  h_next A reference  to a  pointer to  a half-edge  of the
   * other triangle.  The  origin vertex of the half-edge  is the same
   * end vertex of the mesh  edge containing the given geodesic vertex
   * in its interior.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::search_halfedges_on_e(
						  const CP& cp_prev ,
						  const CP& cp_node ,
						  const CP& cp_next ,
						  Halfedge*& h_prev ,
						  Halfedge*& h_next 
						 ) 
    const 
  {
#ifdef DEBUGMODE
    /*
     * Make sure the geodesic vertex is a mesh vertex.
     */
    assert( cp_node.on_edge() ) ;
#endif

    Edge* e = cp_node.get_edge() ;

    Halfedge* h[ 2 ] ;
    _mesh->get_halfedges( e , h[ 0 ] , h[ 1 ] ) ;

    if ( cp_prev.on_vertex() ) {
      if ( cp_prev.get_vertex() != 
	   _mesh->get_org( _mesh->get_prev( h[ 0 ] ) ) ) {
        Halfedge* h_temp = h[ 0 ] ;
        h[ 0 ] = h[ 1 ] ;
        h[ 1 ] = h_temp ;
      }
    } 
    else if ( cp_prev.on_edge() ) {
      if ( 
	  ( cp_prev.get_edge() != 
	    _mesh->get_edge( _mesh->get_next( h[ 0 ] ) ) ) &&
	  ( cp_prev.get_edge() != 
	    _mesh->get_edge( _mesh->get_prev( h[ 0 ] ) ) ) 
	 ) 
      {
        Halfedge* h_temp = h[ 0 ] ;
        h[ 0 ] = h[ 1 ] ;
        h[ 1 ] = h_temp ;
      }
    }
    else {

#ifdef DEBUGMODE
      assert( cp_prev.on_face() ) ;
#endif

      if ( cp_prev.get_face() != _mesh->get_face( h[ 0 ] ) ) {
        Halfedge* h_temp = h[ 0 ] ;
        h[ 0 ] = h[ 1 ] ;
        h[ 1 ] = h_temp ;
      }
    }

    /*
     * At this  point, half-edge h[ 0  ] is the half-edge  of the edge
     * that contains  the current geodesic vertex in  its interior. In
     * addition, h[ 0  ] is in the same mesh face  as the previous and
     * the  current geodesic  vertices,  while  h[ 1  ]  is the  other
     * half-edge of the edge to which half-edge h[ 0 ] belongs.
     */

    h[ 0 ] = _mesh->get_next( h[ 0 ] ) ;

    h_prev = h[ 0 ] ;
    h_next = h[ 1 ] ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::search_halfedges_on_v( CP& cp_prev, CP& cp_node, CP& cp_next, Halfedge*& h_prev, Halfedge*& h_next ) const
   *
   * \brief Determines the two triangles of a vertex star that are cut
   * through  by  a geodesic  that  passes  through  the vertex.   The
   * triangles  contain  the predecessor  and  successor  of the  star
   * vertex  in  the geodesic,  and  they  can  possibly be  the  same
   * triangle.
   *
   * \param cp_prev The previous vertex.
   * \param cp_node The star vertex.
   * \param cp_next The next vertex.
   * \param h_prev A reference to a pointer to a half-edge of the star
   * triangle that is  cut through by the geodesic  when it enters the
   * star. The origin vertex of the half-edge is the star vertex.
   * \param h_next A reference to a pointer to a half-edge of the star
   * triangle that is  cut through by the geodesic  when it leaves the
   * star. The origin vertex of the half-edge is the star vertex.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::search_halfedges_on_v( 
						  CP& cp_prev ,
						  CP& cp_node ,
						  CP& cp_next ,
						  Halfedge*& h_prev ,
						  Halfedge*& h_next 
						 ) 
    const 
  {
#ifdef DEBUGMODE
    /*
     * Make sure the geodesic vertex is a mesh vertex.
     */
    assert( cp_node.on_vertex() ) ;
#endif

    Vertex* v = cp_node.get_vertex() ;

    if ( cp_prev.on_edge() ) {
      /*
       * The previous  geodesic vertex belongs  to the interior  of an
       * edge. This edge  is either incident to the  current vertex or
       * is the  opposite side of  the face that contains  the current
       * vertex.
       */
      Halfedge* h[ 2 ] ;
      _mesh->get_halfedges( cp_prev.get_edge() , h[ 0 ]  , h[ 1 ] ) ;

      if ( _mesh->get_org( h[ 0 ] ) == v ) {
        h_prev = h[ 0 ] ;
      } 
      else if ( _mesh->get_org( h[ 1 ] ) == v ) {
        h_prev = h[ 1 ] ;
      } 
      else if ( _mesh->get_org( _mesh->get_prev( h[ 0 ] ) ) == v ) {
        h_prev = _mesh->get_prev( h[ 0 ] ) ;
      } 
      else {
        h_prev = _mesh->get_prev( h[ 1 ] ) ;
      }
    } 
    else if ( cp_prev.on_vertex() ) {
      Halfedge* h = _mesh->get_halfedge( v ) ;

      while ( _mesh->get_org( _mesh->get_next( h ) ) != cp_prev.get_vertex() ) {
        h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
      }
      
      h_prev = h ;
    }
    else {

#ifdef DEBUGMODE
      assert( cp_prev.on_face() ) ;
#endif
      Halfedge* h = _mesh->get_halfedge( v ) ;

      while ( _mesh->get_face( h ) != cp_prev.get_face() ) {
        h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
      }
      
      h_prev = h ;
    }

    /*
     * The  pointer h_prev must  point to  the half-edge  whose origin
     * vertex is the  current vertex and whose face is  the last one (
     * (in  a counterclockwise traversal  of the  star of  the current
     * vertex) that contains the previous vertex.
     */

    if ( cp_next.on_edge() ) {
      /*
       * The  next  geodesic vertex  belongs  to  the  interior of  an
       * edge. This edge  is either incident to the  current vertex or
       * is the  opposite side of  the face that contains  the current
       * vertex.
       */
      Halfedge* h[ 2 ] ;
      _mesh->get_halfedges( cp_next.get_edge() , h[ 0 ] , h[ 1 ] ) ;
      
      if ( _mesh->get_org( h[ 0 ] ) == v ) {
        h_next = _mesh->get_next( h[ 1 ] ) ;
      } 
      else if ( _mesh->get_org( h[ 1 ] ) == v ) {
        h_next = _mesh->get_next( h[ 0 ] ) ;
      } 
      else if ( _mesh->get_org( _mesh->get_prev( h[ 0 ] ) ) == v ) {
        h_next = _mesh->get_prev( h[ 0 ] ) ;
      } 
      else {
        h_next = _mesh->get_prev( h[ 1 ] ) ;
      }
    } 
    else if ( cp_next.on_vertex() ) {
      Halfedge* h = _mesh->get_halfedge( v ) ;

      while ( _mesh->get_org( _mesh->get_prev( h ) ) != cp_next.get_vertex() ) {
        h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
      }
      
      h_next = h ;
    }
    else {

#ifdef DEBUGMODE
      assert( cp_next.on_face() ) ;
#endif
      Halfedge* h = _mesh->get_halfedge( v ) ;

      while ( _mesh->get_face( h ) != cp_next.get_face() ) {
        h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
      }
      
      h_next = h ;
    }

    /*
     * The  pointer h_next must  point to  the half-edge  whose origin
     * vertex is  the current vertex and  whose face is  the first one
     * (in  a counterclockwise traversal  of the  star of  the current
     * vertex) that contains the next vertex.
     */

#ifdef DEBUGMODE
    assert( h_prev != 0 ) ;
    assert( h_next != 0 ) ;
    assert( _mesh->get_org( h_prev ) == v ) ;
    assert( _mesh->get_org( h_next ) == v ) ;
#endif

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_angles_on_v( Halfedge* h1 , Halfedge* h2 , CP& cp_prev , CP& cp_node , CP& cp_next , double& langle , double& rangle ) const
   *
   * \brief Computes the left and right angles defined by the left and
   * right  sides, respectively,  of a  geodesic at  a  given geodesic
   * vertex. This geodesic vertex is also a mesh vertex.
   *
   * \param h1 Pointer to a  half-edge that identifies one vertex star
   * triangle cut  through by the  geodesic when it enters  the vertex
   * star.
   * \param h2 Pointer to a  half-edge that identifies one vertex star
   * triangle cut  through by the  geodesic when it leaves  the vertex
   * star.
   * \param cp_prev  A reference to the geodesic  vertex that precedes
   * the star vertex in the geodesic.
   * \param cp_node A reference  to the geodesic vertex that coincides
   * with the star vertex.
   * \param cp_next  A reference to the geodesic  vertex that succedes
   * the star vertex in the geodesic.
   * \param langle A reference to the value of the left angle.
   * \param rangle A reference to the value of the right angle.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::compute_angles_on_v( 
						Halfedge* h1 ,
						Halfedge* h2 ,
						CP& cp_prev ,
						CP& cp_node ,
						CP& cp_next ,
						double& langle ,
						double& rangle 
					       )
    const 
  {
    /*
     * The goal is to compute the left and right angles defined by the
     * geodesic  at a mesh  vertex.  This  is a  matter of  adding the
     * angles at the vertex with  respect to each triangle incident to
     * it.  However,  at most  two of these  triangles may  have their
     * interior  cut through  by the  geodesic. When  this  happens, a
     * portion of the angle at the  vertex belongs to the left side of
     * the geodesic and another portion belong to the right side. Both
     * portions are computed first in the following code, and then the
     * angles  at  the  vertex  with  respect to  the  remaining  star
     * triangles  are  found  and   added  to  the  angle  accumulator
     * variables.
     */

    /*
     * Inicialize the values of the left and right angle.
     */
    langle = rangle = 0 ;

    /*
     * Get a pointer to the star mesh vertex.
     */
    Vertex* v = cp_node.get_vertex() ;

    /*
     * If the geodesic vertex preceding the star vertex belongs to the
     * interior of a mesh edge, we determine whether the geodesic cuts
     * through the interior  of the edge or coincides  with it. In the
     * first case, we  compute the left and right  angles the geodesic
     * defines at the star vertex  with respect to the triangle cut by
     * the geodesic.   In the second case,  we let the  right angle be
     * angle at  the star vertex with  respect to the  triangle on the
     * right side of  the geodesic --- the geodesic  coincides with an
     * edge of this triangle.
     */
    if ( cp_prev.on_edge() ) {
      Edge* e_prev = cp_prev.get_edge() ;
      Halfedge* h_prev[ 2 ] ;
      _mesh->get_halfedges( e_prev , h_prev[ 0 ] , h_prev[ 1 ] ) ;

      if ( 
	  ( _mesh->get_org( h_prev[ 0 ] ) == v ) || 
	  ( _mesh->get_org( h_prev[ 1 ] ) == v ) 
	 ) 
      {
        /*
         * The geodesic coincides with the mesh edge that contains the
         * star vertex and the vertex that precedes the star vertex in
         * the geodesic.
         */
        rangle += angle( 
			v , 
			_mesh->get_org( _mesh->get_next( h1 ) ) ,
			_mesh->get_org( _mesh->get_prev( h1 ) ) 
		       ) ;
      }
      else {
        /*
         * The geodesic cuts through the  interior of a mesh edge when
         * entering the star  of the given vertex. So,  we compute the
         * left and  right angles at  the star vertex with  respect to
         * the  star triangle  cut  through by  the  geodesic when  it
         * enters the star.
         */
        langle += angle( 
			v , 
			_mesh->get_org( _mesh->get_next( h1 ) ) ,
			cp_prev 
		       ) ;

        rangle += angle(
			v , 
			_mesh->get_org( _mesh->get_prev( h1 ) ) ,
			cp_prev
		       ) ;
      }
    } 
    else if ( cp_prev.on_vertex() ) {
      /*
       * The vertex  preceding the star  vertex also coincides  with a
       * mesh vertex. So,  we let the right angle be  the angle at the
       * star vertex with respect to  the triangle on the right of the
       * geodesic and  that contains the star vertex  and its previous
       * vertex.
       */
      rangle += angle( 
		      v , 
		      _mesh->get_org( _mesh->get_next( h1 ) ) , 
		      _mesh->get_org( _mesh->get_prev( h1 ) ) 
		     ) ;
    }
    else {

#ifdef DEBUGMODE
      assert( cp_prev.on_face() ) ;
#endif

      /*
       * The geodesic starts at the  interior of a mesh face, and then
       * passes through the given vertex.  So, we compute the left and
       * right angles at the star vertex with respect to two triangles
       * resulting from  splitting the mesh  face by the  line passing
       * through the  geodesic vertex and its  predecessor (inside the
       * face).
       */
      langle += angle( 
		      v ,
		      _mesh->get_org( _mesh->get_next( h1 ) ) ,
		      cp_prev 
		     ) ;

      rangle += angle( 
		      v ,
		      _mesh->get_org( _mesh->get_prev( h1 ) ) ,
		      cp_prev 
		     ) ;
    }

    /*
     * In the following we repeat the same steps as before, but now we
     * consider the  left and  right angles with  respect to  the star
     * vertex and its successor in  the geodesic, and the triangle cut
     * through by  the geodesic  when it leaves  the star  vertex. All
     * cases are the same, but treated with respect to the star vertex
     * successor.
     */
    if ( cp_next.on_edge() ) {
      Edge* e_next = cp_next.get_edge() ;
      Halfedge* h_next[ 2 ] ;
      _mesh->get_halfedges( e_next , h_next[ 0 ] , h_next[ 1 ] ) ;

      if ( 
	  ( _mesh->get_org( h_next[ 0 ] ) == v ) || 
	  ( _mesh->get_org( h_next[ 1 ] ) == v ) 
	 ) 
      {
        rangle += angle( 
			v , 
			_mesh->get_org( _mesh->get_next( h2 ) ) ,
			_mesh->get_org( _mesh->get_prev( h2 ) ) 
		       ) ;
      } 
      else {
        rangle += angle(
			v ,
			_mesh->get_org( _mesh->get_next( h2 ) ) ,
			cp_next 
		       ) ;

        langle += angle(
			v ,
			_mesh->get_org( _mesh->get_prev( h2 ) ) ,
			cp_next 
		       ) ;
      }
    } 
    else if ( cp_next.on_vertex() ) {
      rangle += angle(
		      v ,
		      _mesh->get_org( _mesh->get_next( h2 ) ) , 
		      _mesh->get_org( _mesh->get_prev( h2 ) ) 
		     ) ;
    }
    else {

#ifdef DEBUGMODE
      assert( cp_next.on_face() ) ;
#endif

      /*
       * The geodesic  passes throught the given  geodesic vertex, and
       * then ends at its successor, which is a geodesic vertex in the
       * interior of a  face of the mesh. So, we  compute the left and
       * right angles at the star vertex with respect to two triangles
       * resulting from  splitting the mesh  face by the  line passing
       * through the  geodesic vertex and its  predecessor (inside the
       * face).
       */
      langle += angle( 
		      v ,
		      _mesh->get_org( _mesh->get_next( h2 ) ) ,
		      cp_next
		     ) ;

      rangle += angle( 
		      v ,
		      _mesh->get_org( _mesh->get_prev( h2 ) ) ,
		      cp_next 
		     ) ;
    }

    /*
     * Finally,  we  add  to  the  left and  right  angle  accumulator
     * variables the remaining angles.  Those remaining angles are the
     * ones defined at the star  vertex in the star triangles that are
     * not cut through by the geodesic.
     */
    Halfedge* h = _mesh->get_mate( _mesh->get_prev( h1 ) ) ;
    while ( h != h2 ) {
      rangle += angle( 
		      v , 
		      _mesh->get_org( _mesh->get_next( h ) ) , 
		      _mesh->get_org( _mesh->get_prev( h ) ) 
		     ) ;
      
      h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
    }

    h = _mesh->get_mate( _mesh->get_prev( h2 ) ) ;
    while ( h != h1 ) {
      langle += angle( 
		      v ,
		      _mesh->get_org( _mesh->get_next( h ) ) ,
		      _mesh->get_org( _mesh->get_prev( h ) ) 
		     ) ;

      h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::unfold_triangle_strip( Halfedge* h1 , Halfedge* h2 )
   *
   * \brief Unfold  (in R2)  a strip of  consecutives mesh  faces that
   * belong to the star of a vertex.
   * \param h1 A pointer to a half-edge with origin at the star vertex
   * and contained in the first strip face.
   * \param h2 A pointer to a half-edge with origin at the star vertex
   * and contained in the last strip face.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::unfold_triangle_strip( 
						  Halfedge* h1 ,
						  Halfedge* h2
						 )
  {
#ifdef DEBUGMODE
    /*
     * The half-edges h1 and h2 must share the same origin vertex.
     */
    assert( _mesh->get_org( h1 ) == _mesh->get_org( h2 ) ) ;
#endif

    /*
     * Clear the current list of coordinates.
     */
    _pts.clear() ;

    /*
     * Unfold the first face. The star vertex is mapped to the point A =
     * ( 0 , 0  ), the following vertex is mapped to a point  B = ( xb ,
     * 0), and  the third vertex is  mapped to a point  C = (  xc, yc ),
     * where yc  is positive. We  do not store  A in the output  list of
     * unfolded vertices  (i.e., the parameter pts),  as the coordinates
     * of A are known.
     */
    double xb ;
    double yb = 0 ;
    double xc , yc ;
    unfold_first_triangle( h1 , xb , xc , yc ) ;

    /* 
     * Store the coordinates (xb,yb) and (xc,yc). 
     */
    _pts.push_back( xb ) ;
    _pts.push_back( yb ) ;
    _pts.push_back( xc ) ;
    _pts.push_back( yc ) ;

    /*
     * Unfold the  remaining faces.  Note that the  next unfolded face
     * must share an edge with the previous unfolded face.  So, all we
     * have to do is to compute  the coordinates of only one vertex of
     * the unfolded face.
     */
    Halfedge* h = h1 ;

    while ( h != h2 ) {
      /*
       * Get the next half-edge with the same origin vertex as h1.
       */
      h = _mesh->get_mate( _mesh->get_prev( h ) ) ;

      /*
       * Get  the  planar  coordinates  of  the two  vertices  of  the
       * previous unfolded face that  also belong to the next unfolded
       * face (the one being computed).
       */
      xb = xc ;
      yb = yc ;

      /*
       * Unfold the face containing the half-edge "h".
       */
      unfold_next_triangle( h , xb , yb , xc , yc ) ;

      /*
       * Store the coordinates ( xb , yb ) and ( xc , yc ).
       */
      _pts.push_back( xc ) ;
      _pts.push_back( yc ) ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_2d_coords_on_e( const CP& cp_prev , const CP& cp_node , const CP& cp_next )
   *
   * \brief  Computes  the   planar  Cartesian  coordinates  of  three
   * consecutive  vertices of  this geodesic.   These  coordinates are
   * computed  by unfolding  two  mesh triangles  in  the plane.   The
   * triangles themselves and their  coordinates are given in terms of
   * the input  parameters. The second vertex belongs  to the interior
   * of a mesh edge.
   *
   * \param cp_prev The first of the three consecutive vertices.
   * \param cp_node The second of the three consecutive vertices.
   * \param cp_next The third of the three consecutive vertices.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::compute_2d_coords_on_e( 
						   const CP& cp_prev ,
						   const CP& cp_node ,
						   const CP& cp_next
						  ) 
  {
    Halfedge* h_prev ;
    Halfedge* h_next ;
    get_limits( h_prev , h_next ) ;

#ifdef DEBUGMODE
    assert( h_prev != 0 ) ;
    assert( h_next != 0 ) ;
#endif

#ifdef DEBUGMODE
    assert( _pts.size() == 6 ) ;
#endif

    Halfedge* h[ 2 ] = { _mesh->get_prev( h_prev ) , h_next } ;

    double prev_x ;
    double prev_y ;
    if ( cp_prev.on_vertex() ) {

#ifdef DEBUGMODE
      Halfedge* h_temp = _mesh->get_next( h_prev ) ;
      assert( _mesh->get_org( h_temp ) == cp_prev.get_vertex() ) ;
#endif

      prev_x = _pts[ 0 ] ;
      prev_y = _pts[ 1 ] ;
    } 
    else if ( cp_prev.on_edge() ) {
      if ( _mesh->get_edge( _mesh->get_next( h[ 0 ] ) ) == 
	   cp_prev.get_edge() ) {
        point_on_edge_2d( 
			 cp_prev , 
			 _mesh->get_next( h[ 0 ] ) ,
			 0 ,
			 0 ,
			 _pts[ 0 ] ,
			 _pts[ 1 ] ,
			 prev_x , 
			 prev_y 
			) ;
      } 
      else {
	Halfedge* h_temp = _mesh->get_prev( h[ 0 ] ) ;

#ifdef DEBUGMODE
	assert( _mesh->get_edge( h_temp ) == cp_prev.get_edge() ) ;
#endif

        point_on_edge_2d(
			 cp_prev ,
			 h_temp ,
			 _pts[ 0 ] ,
			 _pts[ 1 ] ,
			 _pts[ 2 ] ,
			 _pts[ 3 ] ,
			 prev_x ,
			 prev_y 
		       ) ;
      }
    }
    else {
#ifdef DEBUGMODE
      assert( cp_prev.on_face() ) ;
      assert( _mesh->get_face( h[ 0 ] ) == cp_prev.get_face() ) ;
#endif

      point_on_face_2d( 
		       cp_prev , 
		       _mesh->get_next( h[ 0 ] ) ,
		       0 ,
		       0 ,
		       _pts[ 0 ] ,
		       _pts[ 1 ] ,
		       _pts[ 2 ] ,
		       _pts[ 3 ] ,
		       prev_x , 
		       prev_y 
		      ) ;
    }

    double next_x ;
    double next_y ;
    if ( cp_next.on_vertex() ) {
#ifdef DEBUGMODE
      Halfedge* h_temp = _mesh->get_prev( h_next ) ;
      assert( _mesh->get_org( h_temp ) == cp_next.get_vertex() ) ;
#endif

      next_x = _pts[ 4 ] ;
      next_y = _pts[ 5 ] ;
    } 
    else if ( cp_next.on_edge() ) {
      if ( _mesh->get_edge( _mesh->get_prev( h[ 1 ] ) ) == 
	   cp_next.get_edge() ) {
        point_on_edge_2d( 
			 cp_next ,
			 _mesh->get_prev( h[ 1 ] ) ,
			 _pts[ 4 ] ,
			 _pts[ 5 ] ,
			 0 ,
			 0 ,
			 next_x , 
			 next_y 
			) ;
      } 
      else {
	Halfedge* h_temp = _mesh->get_next( h[ 1 ] ) ;

#ifdef DEBUGMODE
	assert( _mesh->get_edge( h_temp ) == cp_next.get_edge() ) ;
#endif

        point_on_edge_2d( 
			 cp_next , 
			 h_temp ,
			 _pts[ 2 ] ,
			 _pts[ 3 ] ,
			 _pts[ 4 ] ,
			 _pts[ 5 ] ,
			 next_x ,
			 next_y 
		        ) ;
      }
    }
    else {
#ifdef DEBUGMODE
      assert( cp_next.on_face() ) ;
      assert( _mesh->get_face( h[ 1 ] ) == cp_next.get_face() ) ;
#endif

      point_on_face_2d( 
		       cp_next , 
		       h[ 1 ] ,
		       0 ,
		       0 ,
		       _pts[ 2 ] ,
		       _pts[ 3 ] ,
		       _pts[ 4 ] ,
		       _pts[ 5 ] ,
		       next_x , 
		       next_y 
		      ) ;
    }

    double node_x ;
    double node_y ;
    point_on_edge_2d( 
		     cp_node ,
		     h[ 1 ] ,
		     0 ,
		     0 ,
		     _pts[ 2 ] ,
		     _pts[ 3 ] ,
		     node_x ,
		     node_y 
		    ) ;

    set_2d_coords(
		  prev_x ,
		  prev_y ,
		  node_x ,
		  node_y ,
		  next_x ,
		  next_y
		 ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_2d_coords_on_v( const CP& cp_prev , const CP& cp_node , const CP& cp_next )
   *
   * \brief  Computes  the   planar  Cartesian  coordinates  of  three
   * consecutive  vertices of  this geodesic.   These  coordinates are
   * computed by unfolding mesh triangles in the plane.  The triangles
   * themselves and their coordinates are  given in terms of the input
   * parameters. The second vertex coincides with a mesh vertex,
   *
   * \param cp_prev The first of the three consecutive vertices.
   * \param cp_node The second of the three consecutive vertices.
   * \param cp_next The third of the three consecutive vertices.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::compute_2d_coords_on_v( 
						   const CP& cp_prev , 
						   const CP& cp_node ,
						   const CP& cp_next
						  ) 
  {
    Halfedge* h_prev ;
    Halfedge* h_next ;
    get_limits( h_prev , h_next ) ;

#ifdef DEBUGMODE
    assert( h_prev != 0 ) ;
    assert( h_next != 0 ) ;
    assert( _mesh->get_org( h_prev ) == cp_node.get_vertex() ) ;
    assert( _mesh->get_org( h_next ) == cp_node.get_vertex() ) ;
#endif

#ifdef DEBUGMODE
    assert(   _pts.size()       >= 6 ) ;
    assert( ( _pts.size() % 2 ) == 0 ) ;
#endif

    double prev_x ;
    double prev_y ;
    if ( cp_prev.on_vertex() ) {
#ifdef DEBUGMODE
      Halfedge* h_temp = _mesh->get_next( h_prev ) ;
      assert( _mesh->get_org( h_temp ) == cp_prev.get_vertex() ) ;
#endif
      prev_x = _pts[ 0 ] ;
      prev_y = _pts[ 1 ] ;
    } 
    else if ( cp_prev.on_edge() ) {
      if ( _mesh->get_edge( h_prev ) == cp_prev.get_edge() ) {
        point_on_edge_2d( 
			 cp_prev ,
			 h_prev ,
			 0 ,
			 0 ,
			 _pts[ 0 ] ,
			 _pts[ 1 ] ,
			 prev_x ,
			 prev_y 
			) ;
      } 
      else {
	Halfedge* h_temp = _mesh->get_next( h_prev ) ;

#ifdef DEBUGMODE
	assert( _mesh->get_edge( h_temp ) == cp_prev.get_edge() ) ;
#endif

        point_on_edge_2d( 
			 cp_prev , 
			 h_temp , 
			 _pts[ 0 ] , 
			 _pts[ 1 ] ,
			 _pts[ 2 ] , 
			 _pts[ 3 ] , 
			 prev_x , 
			 prev_y 
			) ;
      }
    }
    else {
#ifdef DEBUGMODE
      assert( cp_prev.on_face() ) ;
      assert( _mesh->get_face( h_prev ) == cp_prev.get_face() ) ;
#endif

      point_on_face_2d( 
		       cp_prev ,
		       h_prev ,
		       0 ,
		       0 ,
		       _pts[ 0 ] , 
		       _pts[ 1 ] ,
		       _pts[ 2 ] , 
		       _pts[ 3 ] , 
		       prev_x , 
		       prev_y 
		      ) ;      
    }

    double next_x ;
    double next_y ;
    if ( cp_next.on_vertex() ) {
#ifdef DEBUGMODE
      Halfedge* h_temp = _mesh->get_prev( h_next ) ;
      assert( _mesh->get_org( h_temp ) == cp_next.get_vertex() ) ;
#endif

      next_x = _pts[ _pts.size() - 2 ] ;
      next_y = _pts[ _pts.size() - 1 ] ;
    } 
    else if ( cp_next.on_edge() ) {
      if ( _mesh->get_edge( _mesh->get_prev( h_next ) ) == 
	   cp_next.get_edge() ) {
        point_on_edge_2d( 
			 cp_next ,
			 _mesh->get_prev( h_next ) ,
			 _pts[ _pts.size() - 2 ] ,
			 _pts[ _pts.size() - 1 ] ,
			 0 ,
			 0 ,
			 next_x ,
			 next_y 
			) ;
      } 
      else {
	Halfedge* h_temp = _mesh->get_next( h_next ) ;

#ifdef DEBUGMODE
	assert( _mesh->get_edge( h_temp ) == cp_next.get_edge() ) ;
#endif

        point_on_edge_2d(
			 cp_next , 
			 h_temp ,
			 _pts[ _pts.size() - 4 ] ,
			 _pts[ _pts.size() - 3 ] ,
			 _pts[ _pts.size() - 2 ] ,
			 _pts[ _pts.size() - 1 ] ,
			 next_x , 
			 next_y 
		       ) ;
      }
    }
    else {
#ifdef DEBUGMODE
      assert( cp_next.on_face() ) ;
      assert( _mesh->get_face( h_next ) == cp_prev.get_face() ) ;
#endif

      point_on_face_2d( 
		       cp_next ,
		       h_next ,
		       0 ,
		       0 ,
		       _pts[ _pts.size() - 4 ] ,
		       _pts[ _pts.size() - 3 ] ,
		       _pts[ _pts.size() - 2 ] ,
		       _pts[ _pts.size() - 1 ] ,
		       next_x , 
		       next_y 
		      ) ;
      
    }

    set_2d_coords(
		  prev_x ,
		  prev_y ,
		  0 ,
		  0 ,
		  next_x ,
		  next_y
		 ) ;
		  
    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::compute_intersection_on_e( CP& cp_prev , CP& cp_node , CP& cp_next ) const
   *
   * \brief Computes the intersection  of two line segments. The first
   * one  is an  unfolded mesh  edge that  contains a  geodesic vertex
   * (named  the current  vertex). The  second one  is defined  by the
   * unfolded  geodesics  vertices  that  precedes  and  succeeds  the
   * current vertex.
   *
   * \param cp_prev  A reference to  the geodesic vertex  the precedes
   * the vertex that belongs to the interior of a mesh edge.
   * \param  cp_node A  reference to  the vertex  that belongs  to the
   * interior of a mesh edge.
   * \param cp_next  A reference to  the geodesic vertex  the succeeds
   * the vertex that belongs to the interior of a mesh edge.
   *
   * \return The logic  value true in two cases:  (a) the intersection
   * point is  a point in the  interior of both  segments and distinct
   * from the position of the current geodesic vertex; (b) there is no
   * intersection   between  the   interior  of   the   two  segments.
   * Otherwise, the logic value false is returned.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::compute_intersection_on_e(
						      CP& cp_prev ,
						      CP& cp_node ,
						      CP& cp_next 
						     ) 
    const
  {

#ifdef DEBUGMODE
    assert( cp_node.on_edge() ) ;
#endif

    Halfedge* h_prev ;
    Halfedge* h_next ;
    get_limits( h_prev , h_next ) ;

    double prev_x , 
           prev_y , 
           node_x , 
           node_y , 
           next_x , 
           next_y ;

    get_2d_coords( 
		  prev_x ,
		  prev_y ,
		  node_x ,
		  node_y , 
		  next_x ,
		  next_y 
		 ) ;

#ifdef DEBUGMODE
    double inter , temp ;
#else
    double inter = 0 , temp = 0 ;
#endif

    INTERSECTION_TYPE ret = 
      Geometric::compute_segment_intersection( 
					      prev_x , 
					      prev_y , 
					      next_x , 
					      next_y , 
					      0 , 
					      0 , 
					      _pts[ 2 ] ,
					      _pts[ 3 ] ,
					      temp ,
					      inter 
					     ) ;

#ifdef DEBUGMODE
    assert(
        ( ret == Geometric::POINT ) ||
        ( ret == Geometric::NOINTERSECTION )
    ) ;
#endif
    
    Halfedge* h[ 2 ] = { _mesh->get_prev( h_prev ) , h_next } ;
    
    if ( ret == Geometric::POINT ) {
      Halfedge* h_node[ 2 ] ;
      _mesh->get_halfedges( 
			   cp_node.get_edge() , 
			   h_node[ 0 ] , 
			   h_node[ 1 ] 
			  ) ;
      
      if ( h_node[ 0 ] == h[ 0 ] ) {
	/*
	 * We check whether  the vertex has moved up  to a pre-defined
	 * tolerance. If so, we update the vertex position. Otherwise,
	 * we return  the logic value  false, so that we  don't update
	 * the heap unnecessarily.
	 */
        if ( fabs( cp_node.get_t() - ( 1 - inter ) ) < _TOL_T_EDGE ) {
          return false ;
        }

        cp_node.set_edge( cp_node.get_edge() , 1 - inter ) ;
      } 
      else {
	/*
	 * We check whether  the vertex has moved up  to a pre-defined
	 * tolerance. If so, we update the vertex position. Otherwise,
	 * we return  the logic value  false, so that we  don't update
	 * the heap unnecessarily.
	 */
        if ( fabs( cp_node.get_t() - inter ) < _TOL_T_EDGE ) {
          return false ;
        }

        cp_node.set_edge( cp_node.get_edge() , inter ) ;
      }
    } 
    else {
      /*
       * The geodesic vertex now coincides with a mesh vertex.
       */
      if ( inter > 1.0 ) {
	GP* gp = cp_node.get_point() ;
#ifdef DEBUGMODE
	assert( gp != 0 ) ;
#endif
	delete gp ;
        cp_node.set_point( new GPV( _mesh->get_org( h[ 0 ] ) ) ) ;
      } 
      else {
	GP* gp = cp_node.get_point() ;
#ifdef DEBUGMODE
	assert( gp != 0 ) ;
#endif
	delete gp ;
        cp_node.set_point( new GPV( _mesh->get_org( h[ 1 ] ) ) ) ;
      }
    }

    return true ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::compute_intersection_on_v( CP& cp_prev , CP& cp_node , CP& cp_next , Path& new_path  ) const
   *
   * \brief Computes  the intersection points  of a line  segment with
   * several other line segments.  The first line segment connects the
   * previous  and next  vertices  of a  given  geodesic vertex.   The
   * geodesic  vertex  coincided with  a  mesh  vertex. The  remaining
   * segments  are edges  connecting  the given  geodesic vertex  with
   * other vertices of  the mesh. The intersection is  computed in the
   * plane, after unfolding  a strip of triangles of  the given vertex
   * star.
   *
   * \param cp_prev  A reference to  the geodesic vertex  the precedes
   * the given star vertex.
   * \param cp_node A reference to the given star vertex.
   * \param cp_next  A reference to  the geodesic vertex  the succeeds
   * the given star vertex.
   * \param  new_path A  reference  to a  list  with the  intersection
   * points.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::compute_intersection_on_v(
						      CP& cp_prev ,
						      CP& cp_node ,
						      CP& cp_next , 
						      Path& new_path 
						     ) 
    const 
  {
    Halfedge* h_prev ;
    Halfedge* h_next ;
    get_limits( h_prev , h_next ) ;

    double prev_x , 
           prev_y , 
           node_x , 
           node_y , 
           next_x , 
           next_y ;

    get_2d_coords( 
		  prev_x ,
		  prev_y ,
		  node_x ,
		  node_y , 
		  next_x ,
		  next_y 
		 ) ;

    Halfedge* h = _mesh->get_mate( _mesh->get_prev( h_prev ) ) ;
    int i = 0 ;
    do {
      double inter , temp ;

      INTERSECTION_TYPE ret = 
	Geometric::compute_segment_intersection( 
						prev_x ,
						prev_y , 
						next_x , 
						next_y , 
						0 , 
						0 , 
						_pts[ 2 * i + 2 ] ,
						_pts[ 2 * i + 3 ] ,
						temp , 
						inter 
					       ) ;

#ifdef DEBUGMODE
      assert(
          ( ret == Geometric::POINT ) ||
          ( ret == Geometric::NOINTERSECTION )
      ) ;
#endif
      
      CP cp_new ;
      if ( ret == Geometric::POINT ) {
        Halfedge* h_node[ 2 ] ;
        _mesh->get_halfedges( 
			     _mesh->get_edge( h ) , 
			     h_node[ 0 ] , 
			     h_node[ 1 ] 
			    ) ;

        if ( h_node[ 0 ] == h ) {
	  GPE* gpe = new GPE( _mesh->get_edge( h ) , inter ) ;
          cp_new.set_point( gpe ) ;
        } 
	else {
	  GPE* gpe = new GPE( _mesh->get_edge( h ) , 1 - inter ) ;
          cp_new.set_point( gpe ) ;
        }
      } 
      else {
	GPV* gpv = new GPV( _mesh->get_org( _mesh->get_next( h ) ) ) ;
        cp_new.set_point( gpv ) ;
      }

      new_path.push_back( cp_new ) ;

      h = _mesh->get_mate( _mesh->get_prev( h ) ) ;
      i++ ;
    } 
    while ( h != _mesh->get_mate( _mesh->get_prev( h_next ) ) ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::find_common_faces( Vertex* v1 , Vertex* v2 , Face*& f1 , Face*& f2 ) const
   *
   * \brief Finds the two faces  incident to two distinct vertices (if
   * any).
   *
   * \param v1 Pointer to a mesh vertex.
   * \param v2 Pointer to a mesh vertex.
   * \param f1 Reference to a pointer to a mesh face.
   * \param f2 Reference to a pointer to a mesh face.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::find_common_faces(
					      Vertex* v1 , 
					      Vertex* v2 ,
					      Face*& f1 ,
					      Face*& f2
					     ) 
    const 
  {
#ifdef DEBUGMODE
    /*
     * Vertices cannot be the same.
     */
    assert( v1 != v2 ) ;
#endif

    f1 = 0 ;
    f2 = 0 ;

    Halfedge* h = _mesh->get_halfedge( v1 ) ;
    Halfedge* haux = h ;
    do {
      if ( _mesh->get_org( _mesh->get_next( haux ) ) == v2 ) {
	f1 = _mesh->get_face( haux ) ;
	f2 = _mesh->get_face( _mesh->get_mate( haux ) ) ;
	return ;
      }

      haux = _mesh->get_mate( _mesh->get_prev( haux ) ) ;
    }
    while ( haux != h ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::find_common_faces( Vertex* v , Edge* e , Face*& f1 , Face*& f2 ) const
   *
   * \brief  Finds the  face(s) incident  to a  mesh edge  and  a mesh
   * vertex (if any).
   *
   * \param e Pointer to a mesh edge.
   * \param v Pointer to a mesh vertex.
   * \param f1 Reference to a pointer to a mesh face.
   * \param f2 Reference to a pointer to a mesh face.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::find_common_faces(
					      Vertex* v ,
					      Edge* e , 
					      Face*& f1 ,
					      Face*& f2 
					     ) 
    const 
  {
#ifdef DEBUGMODE
    assert( v != 0 ) ;
    assert( e != 0 ) ;
#endif

    /*
     * Initialize the pointer to the face we want to find.
     */
    f1 = 0 ;
    f2 = 0 ;

    Halfedge* he1 ;
    Halfedge* he2 ;
    _mesh->get_halfedges( e , he1 , he2 ) ;

    if ( 
	( _mesh->get_org( he1 ) == v ) || 
	( _mesh->get_org( he2 ) == v ) 
       ) {
      f1 = _mesh->get_face( he1 ) ;
      f2 = _mesh->get_face( he2 ) ;
    }
    else {
      Halfedge* h = _mesh->get_halfedge( v ) ;
      Halfedge* haux = h ;
      do {
	if ( _mesh->get_edge( _mesh->get_next( haux) ) == e ) {
	  f1 = _mesh->get_face( haux ) ;
	  return ;
	}

	haux = _mesh->get_mate( _mesh->get_prev( haux ) ) ;
      }
      while ( haux != h ) ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::find_common_faces( Edge* e1 , Edge* e2 , Face*& f1, Face*& f2 ) const
   *
   * \brief Finds the face(s) incident to two mesh edges (if any).
   *
   * \param e1 Pointer to a mesh edge.
   * \param e2 Pointer to a mesh edge.
   * \param f1 Reference to a pointer to a mesh face.
   * \param f2 Reference to a pointer to a mesh face.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::find_common_faces(
					      Edge* e1 ,
					      Edge* e2 , 
					      Face*& f1 ,
					      Face*& f2 
					     ) 
    const 
  {
#ifdef DEBUGMODE
    assert( e1 != 0 ) ;
    assert( e2 != 0 ) ;
#endif

    /*
     * Initialize the pointer to the face we want to find.
     */
    f1 = 0 ;
    f2 = 0 ;

    Halfedge* he1 ;
    Halfedge* he2 ;
    _mesh->get_halfedges( e1 , he1 , he2 ) ;

    if ( e1 == e2 ) {
      f1 = _mesh->get_face( he1 ) ;
      f2 = _mesh->get_face( he2 ) ;
    }
    else {
      Halfedge* he3 ;
      Halfedge* he4 ;
      _mesh->get_halfedges( e2 , he3 , he4 ) ;
	
      Face* f3 = _mesh->get_face( he1 ) ;
      Face* f4 = _mesh->get_face( he2 ) ;
      Face* f5 = _mesh->get_face( he3 ) ;
      Face* f6 = _mesh->get_face( he4 ) ;

      if ( ( f3 == f5 ) || ( f3 == f6 ) ) {
	f1 = f3 ;
      }
      else if ( ( f4 == f5 ) || ( f4 == f6 ) ) {
	f1 = f4 ;
      }
    }

    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::is_in_face( Vertex* v , Face* f ) const
   *
   * \brief Decides  whether a  given mesh vertex  belongs to  a given
   * mesh face.
   *
   * \param v Pointer to a mesh vertex.
   * \param f Pointer to a mesh face.
   *
   * \return The logic  value true if the given  vertex belongs to the
   * given face, and the logic value false otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::is_in_face(
				       Vertex* v , 
				       Face* f
				      ) 
    const 
  {
#ifdef DEBUGMODE
    assert( v != 0 ) ;
    assert( f != 0 ) ;
#endif

    Halfedge* h1 = _mesh->get_halfedge( f ) ;
    if ( _mesh->get_org( h1 ) == v ) {
      return true ;
    }

    h1 = _mesh->get_next( h1 ) ;
    if ( _mesh->get_org( h1 ) == v ) {
      return true ;
    }

    h1 = _mesh->get_next( h1 ) ;
    if ( _mesh->get_org( h1 ) == v ) {
      return true ;
    }

    return false ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::is_in_face( Edge* e , Face* f ) const
   *
   * \brief Decides whether a given  mesh edge belongs to a given mesh
   * face.
   *
   * \param e Pointer to a mesh edge.
   * \param f Pointer to a mesh face.
   *
   * \return The  logic value  true if the  given edge belongs  to the
   * given face, and the logic value false otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::is_in_face(
				       Edge* e , 
				       Face* f
				      ) 
    const 
  {
#ifdef DEBUGMODE
    assert( e != 0 ) ;
    assert( f != 0 ) ;
#endif

    Halfedge* h1 = _mesh->get_halfedge( f ) ;
    if ( _mesh->get_edge( h1 ) == e ) {
      return true ;
    }

    h1 = _mesh->get_next( h1 ) ;
    if ( _mesh->get_edge( h1 ) == e ) {
      return true ;
    }

    h1 = _mesh->get_next( h1 ) ;
    if ( _mesh->get_edge( h1 ) == e ) {
      return true ;
    }

    return false ;
  }


  /**
   * \fn double Approx_geodesics< Mesh >::angle( Vertex* v1 , Vertex* v2 , Vertex* v3 ) const
   *
   * \brief Calculates  the angle  defined at the  first of  the three
   * given mesh vertices of a mesh triangle.
   *
   * \param v1 The first vertex of a mesh triangle.
   * \param v2 The second vertex of a mesh triangle.
   * \param v3 The third vertex of a mesh triangle.
   *
   * \return The angle at v1 defined by  the vectors ( v2 - v1 ) and (
   * v3 - v1 ).
   */
  template< typename Mesh >
  double 
  Approx_geodesics< Mesh >::angle( 
				  Vertex* v1 ,
				  Vertex* v2 ,
				  Vertex* v3 
				 ) 
    const 
  {
    double x[ 3 ] , y[ 3 ] , z[ 3 ] ;
    
    _mesh->get_coords( v1 , x[ 0 ] , y[ 0 ] , z[ 0 ] ) ;
    _mesh->get_coords( v2 , x[ 1 ] , y[ 1 ] , z[ 1 ] ) ;
    _mesh->get_coords( v3 , x[ 2 ] , y[ 2 ] , z[ 2 ] ) ;

    return Geometric::get_angle(
				x[ 0 ] ,
				y[ 0 ] ,
				z[ 0 ] ,
				x[ 1 ] ,
				y[ 1 ] ,
				z[ 1 ] ,
				x[ 2 ] ,
				y[ 2 ] ,
				z[ 2 ] 
			       ) ;
  }


  /**
   * \fn double Approx_geodesics< Mesh >::angle( Vertex* v1 , Vertex* v2 , CP& cp ) const
   *
   * \brief Calculates  the angle  defined at the  first of  the three
   * given  geodesic  vertices.   The  first  two  vertices  are  mesh
   * vertices, and  the third  one may  be either a  mesh vertex  or a
   * point in the interior of a mesh edge.
   *
   * \param v1 The first geodesic vertex.
   * \param v2 The second geodesic vertex.
   * \param cp The third geodesic vertex.
   *
   * \return The angle at  v1 defined by the vector ( v2  - v1 ) and (
   * cp - v1 ).
   */
  template< typename Mesh >
  double 
  Approx_geodesics< Mesh >::angle( 
				  Vertex* v1 ,
				  Vertex* v2 ,
				  CP& cp 
				 ) 
    const 
  {
    double x[ 3 ] , y[ 3 ] , z[ 3 ] ;

    _mesh->get_coords( v1 , x[ 0 ] , y[ 0 ] , z[ 0 ] ) ;
    _mesh->get_coords( v2 , x[ 1 ] , y[ 1 ] , z[ 1 ] ) ;

    cp.get_point()->get_coords( _mesh , x[ 2 ] , y[ 2 ] , z[ 2 ] ) ;
    
    return Geometric::get_angle( 
				x[ 0 ] ,
				y[ 0 ] ,
				z[ 0 ] ,
				x[ 1 ] ,
				y[ 1 ] ,
				z[ 1 ] ,
				x[ 2 ] ,
				y[ 2 ] ,
				z[ 2 ] 
			       ) ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::unfold_first_triangle( Halfedge* h , double& xb , double& xc , double& yc ) const
   *
   * \brief Unfold a mesh face in the plane. The face is identified by
   * one of  its half-edges.  The origin vertex  of this  half-edge is
   * mapped to  the point (0,0), the  following vertex is  mapped to a
   * point in the X axis, and the third vertex is mapped to a point in
   * the upper-half of the plane.
   *
   * \param h  A pointer to a half-edge.
   * \param xb A reference to the X coordinate of the second Cartesian
   * vertex of the unfolded mesh face.
   * \param xc A reference to  the X coordinate of the third Cartesian
   * vertex of the unfolded mesh face.
   * \param yc A reference to  the Y coordinate of the third Cartesian
   * vertex of the unfolded mesh face.
   */
  template< typename Mesh >
  void
  Approx_geodesics< Mesh >::unfold_first_triangle(
						  Halfedge* h ,
						  double& xb ,
						  double& xc ,
						  double& yc 
						 ) 
    const 
  {
    /*
     * Get the 3D coordinates of the first vertex of the mesh face.
     */
    double x1 , y1 , z1 ;
    _mesh->get_coords( _mesh->get_org( h ) , x1 , y1 , z1 ) ;

    /*
     * Get the 3D coordinates of the second vertex of the mesh face.
     */
    double x2 , y2 , z2 ;
    _mesh->get_coords( 
		      _mesh->get_org( _mesh->get_next( h ) ) , 
		      x2 , 
		      y2 , 
		      z2 
		     ) ;

    /*
     * Get the 3D coordinates of the third vertex of the mesh face.
     */
    double x3 , y3 , z3 ;
    _mesh->get_coords(
		      _mesh->get_org( _mesh->get_prev( h ) ) , 
		      x3 , 
		      y3 , 
		      z3 
		     ) ;

    /*
     * Compute the distance from the first to the third vertex.
     */
    double l13 = sqrt( 
		      ( x3 - x1 ) * ( x3 - x1 ) + 
		      ( y3 - y1 ) * ( y3 - y1 ) + 
		      ( z3 - z1 ) * ( z3 - z1 ) 
		     ) ;

    /*
     * Compute the distance from the second to the third vertex.
     */
    double l23 = sqrt( 
		      ( x3 - x2 ) * ( x3 - x2 ) + 
		      ( y3 - y2 ) * ( y3 - y2 ) + 
		      ( z3 - z2 ) * ( z3 - z2 ) 
		     ) ;

    /*
     * Set the value  of the X coordinate of the  second vertex of the
     * unfolded triangle.
     */
    xb = sqrt( 
	      ( x2 - x1 ) * ( x2 - x1 ) + 
	      ( y2 - y1 ) * ( y2 - y1 ) + 
	      ( z2 - z1 ) * ( z2 - z1 ) 
	     ) ;

    /*
     * Compute  the X and  Y coordinates  of the  third vertex  of the
     * unfolded triangle.
     */
    Geometric::compute_third_vertex_position( 
					     0 ,
					     0 , 
					     xb ,
					     0 ,
					     l13 ,
					     l23 ,
					     xc ,
					     yc 
					    ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::unfold_next_triangle( Halfedge* h , double xb , double yb , double& xc , double& yc ) const
   *
   * \brief Unfold a mesh face in the plane. The face is identified by
   * one of  its half-edges.  The origin vertex  of this  half-edge is
   * mapped to  the point (0,0), the  following vertex is  mapped to a
   * given point, and  the third vertex is mapped  to a point computed
   * by this method.
   *
   * \param h  A pointer to a half-edge.
   * \param xb The X coordinate  of the second Cartesian vertex of the
   * unfolded mesh face.
   * \param yb The Y coordinate  of the second Cartesian vertex of the
   * unfolded mesh face.
   * \param xc A reference to  the X coordinate of the third Cartesian
   * vertex of the unfolded mesh face.
   * \param yc A reference to  the Y coordinate of the third Cartesian
   * vertex of the unfolded mesh face.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::unfold_next_triangle( 
						 Halfedge* h ,
						 double xb ,
						 double yb ,
						 double& xc ,
						 double& yc 
						) 
    const 
  {
    /*
     * Get the 3D coordinates of the first vertex of the mesh face.
     */
    double x1 , y1 , z1 ;
    _mesh->get_coords( _mesh->get_org( h ) , x1 , y1 , z1 ) ;

    /*
     * Get the 3D coordinates of the second vertex of the mesh face.
     */
    double x2 , y2 , z2 ;
    _mesh->get_coords( 
		      _mesh->get_org( _mesh->get_next( h ) ) ,
		      x2 ,
		      y2 ,
		      z2 
		     ) ;

    /*
     * Get the 3D coordinates of the third vertex of the mesh face.
     */
    double x3 , y3 , z3 ;
    _mesh->get_coords( 
		      _mesh->get_org( _mesh->get_prev( h ) ) , 
		      x3 , 
		      y3 , 
		      z3 
		     ) ;

    /*
     * Compute the distance from the first to the third vertex.
     */
    double l13 = sqrt( 
		      ( x3 - x1 ) * ( x3 - x1 ) + 
		      ( y3 - y1 ) * ( y3 - y1 ) +
		      ( z3 - z1 ) * ( z3 - z1 ) 
		     ) ;

    /*
     * Compute the distance from the second to the third vertex.
     */
    double l23 = sqrt( 
		      ( x3 - x2 ) * ( x3 - x2 ) + 
		      ( y3 - y2 ) * ( y3 - y2 ) + 
		      ( z3 - z2 ) * ( z3 - z2 ) 
		     ) ;

    /*
     * Compute  the X and  Y coordinates  of the  third vertex  of the
     * unfolded triangle.
     */
    Geometric::compute_third_vertex_position( 
					     0 ,
					     0 ,
					     xb ,
					     yb ,
					     l13 ,
					     l23 ,
					     xc ,
					     yc 
					    ) ;

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::point_on_edge_2d( const CP& cp , Halfedge* h , double x1 , double y1 , double x2 , double y2 , double& x , double& y ) const
   *
   * \brief Computes the planar  coordinates of a geodesic vertex that
   * belongs to  the interior  of a mesh  edge. These  coordinates are
   * given with respect  to the unfold of the  two mesh triangles that
   * share  the  mesh  edge   containing  the  geodesic  vertex.   All
   * information regarding  the unfolded triangles  are already cached
   * in the geodesic vertex.
   *
   * \param cp A reference to a geodesic vertex.
   * \param h Pointer  to a half-edge of the  mesh edge containing the
   * geodesic vertex.
   * \param  x1  The X  planar  coordinate of  one  end  point of  the
   * unfolded mesh edge containing the geodesic vertex.
   * \param  y1  The Y  planar  coordinate of  one  end  point of  the
   * unfolded mesh edge containing the geodesic vertex.
   * \param x2 The  X planar coordinate of the other  end point of the
   * unfolded mesh edge containing the geodesic vertex.
   * \param y2 The  Y planar coordinate of the other  end point of the
   * unfolded mesh edge containing the geodesic vertex.
   * \param  x A  reference to  the  X planar  coordinate of  geodesic
   * vertex.
   * \param  y A  reference to  the  Y planar  coordinate of  geodesic
   * vertex.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::point_on_edge_2d( 
					     const CP& cp ,
					     Halfedge* h ,
					     double x1 ,
					     double y1 , 
					     double x2 ,
					     double y2 ,
					     double& x ,
					     double& y 
					    ) 
    const 
  {
    Halfedge* h_cp[ 2 ] ;
    _mesh->get_halfedges( cp.get_edge() , h_cp[ 0 ] , h_cp[ 1 ] ) ;

    if ( h_cp[ 0 ] == h ) {
      x = x1 + cp.get_t() * ( x2 - x1 ) ;
      y = y1 + cp.get_t() * ( y2 - y1 ) ;
    } else {
      x = x2 + cp.get_t() * ( x1 - x2 ) ;
      y = y2 + cp.get_t() * ( y1 - y2 ) ;
    }

    return ;
  }


  /**
   * \fn void Approx_geodesics< Mesh >::point_on_face_2d( const CP& cp , Halfedge* h , double x1 , double y1 , double x2 , double y2 , double x3 , double y3 , double& x , double& y ) const
   *
   * \brief Computes the planar  coordinates of a geodesic vertex that
   * belongs to  the interior  of a mesh  face. These  coordinates are
   * given  with  respect  to  the  unfold  of  the  mesh  face.   All
   * information  regarding  the  unfolded  face  are  cached  in  the
   * geodesic vertex.
   *
   * \param cp A reference to a geodesic vertex.
   * \param h  Pointer to  the halfedge whose  origin is  the geodesic
   * vertex containing the cached unfolded face information.
   * \param  x1 The X  planar coordinate  of one  vertex of  mesh face
   * containing the geodesic vertex.
   * \param y1 The  Y planar coordinate of of one  vertex of mesh face
   * containing the geodesic vertex.
   * \param  x2 The  X planar  coordinate of  a second  vertex  of the
   * unfolded mesh face containing the geodesic vertex.
   * \param  y2 The  Y planar  coordinate of  a second  vertex  of the
   * unfolded mesh face containing the geodesic vertex.
   * \param  x3 The  X planar  coordinate of  a second  vertex  of the
   * unfolded mesh face containing the geodesic vertex.
   * \param  y3 The  Y planar  coordinate of  a second  vertex  of the
   * unfolded mesh face containing the geodesic vertex.
   * \param  x A  reference to  the  X planar  coordinate of  geodesic
   * vertex.
   * \param  y A  reference to  the  Y planar  coordinate of  geodesic
   * vertex.
   */
  template< typename Mesh >
  void 
  Approx_geodesics< Mesh >::point_on_face_2d( 
					     const CP& cp ,
					     Halfedge* h ,
					     double x1 ,
					     double y1 , 
					     double x2 ,
					     double y2 ,
					     double x3 ,
					     double y3 ,
					     double& x ,
					     double& y 
					    ) 
    const 
  {
#ifdef DEBUGMODE
    assert( cp.on_face() ) ;
#endif

    Face* f = cp.get_face() ;

    double u ;
    double v ;
    double w ;
    if ( h == _mesh->get_halfedge( f ) ) {
      u = cp.get_u() ;
      v = cp.get_v() ;
      w = 1 - u - v ;
      if ( fabs(     w ) < 1e-14 ) w = 0 ;
      if ( fabs( 1 - w ) < 1e-14 ) w = 1 ;
    }
    else if ( _mesh->get_next( h ) == _mesh->get_halfedge( f ) ) {
      v = cp.get_u() ;
      w = cp.get_v() ;
      u = 1 - v - w ;
      if ( fabs(     u ) < 1e-14 ) u = 0 ;
      if ( fabs( 1 - u ) < 1e-14 ) u = 1 ;
    }
    else {
      w = cp.get_u() ;
      u = cp.get_v() ;
      v = 1 - w - u ;
      if ( fabs(     v ) < 1e-14 ) v = 0 ;
      if ( fabs( 1 - v ) < 1e-14 ) v = 1 ;
    }
    
    x = u * x1 + v * x2 + w * x3 ;
    y = u * y1 + v * y2 + w * y3 ;

    return ;
  }


  /**
   * \fn bool Approx_geodesics< Mesh >::is_consistent() const
   *
   * \brief Checks  whether the  current geodesic is  consistent.  The
   * geodesic is considered consistent if each pair of its consecutive
   * vertices are incident to the same  face, but not to the same edge
   * or vertex.
   *
   * \return The logic  value true if the geodesic  is consistent, and
   * the logic value false otherwise.
   */
  template< typename Mesh >
  bool
  Approx_geodesics< Mesh >::is_consistent() const
  {    
    typename std::list< CP >::const_iterator pit1 = _cp_path.begin() ;
    while ( pit1 != _cp_path.end() ) {
      typename std::list< CP >::const_iterator pit2 = pit1 ;
      ++pit2 ;
      if ( pit2 != _cp_path.end() ) {
	const CP& node = *pit1 ;
	const CP& next = *pit2 ;
	if ( node.on_vertex() ) {
	  Vertex* v_node = node.get_vertex() ;
	  if ( next.on_vertex() ) {
	    Face* f1 ;
	    Face* f2 ;
	    find_common_faces( v_node , next.get_vertex() , f1 , f2 ) ;
	    if ( f1 == 0 ) {
	      return false ;
	    }
	  }
	  else if ( next.on_edge() ) {
	    Face* f1 ;
	    Face* f2 ;
	    find_common_faces( v_node , next.get_edge() , f1 , f2 ) ;
	    if ( f1 == 0 ) {
	      return false ;
	    }
	  }
	  else if ( !is_in_face( v_node , next.get_face() ) ) {
	    return false ;
	  }
	}
	else if ( node.on_edge() ) {
	  Edge* e_node = node.get_edge() ;
	  if ( next.on_vertex() ) {
	    Face* f1 ;
	    Face* f2 ;
	    find_common_faces( next.get_vertex() , e_node , f1 , f2 ) ;
	    if ( f1 == 0 ) {
	      return false ;
	    }
	  }
	  else if ( next.on_edge() ) {
	    Face* f1 ;
	    Face* f2 ;
	    find_common_faces( next.get_edge() , e_node , f1 , f2 ) ;
	    if ( f1 == 0 ) {
	      return false ;
	    }
	  }
	  else if ( !is_in_face( e_node , next.get_face() ) ) {
	    return false ;
	  }
	}
	else {
	  Face* f_node = node.get_face() ;
	  if ( next.on_vertex() ) {
	    if ( !is_in_face( next.get_vertex() , f_node ) ) {
	      return false ;
	    }
	  }
	  else if ( next.on_edge() ) {
	    if ( !is_in_face( next.get_edge() , f_node ) ) {
	      return false ;
	    }
	  }
	  else if ( f_node != next.get_face() ) {
	    return false ;
	  }
	}
      }
      else {
	break ;
      }
      
      ++pit1 ;
    }

    return true ;
  }

}

/** @} *///end of group class.

#endif	/* ADG_GEODESICS_HPP */
