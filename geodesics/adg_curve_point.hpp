/** 
 * \file adg_curve_point.hpp
 *  
 * \brief Definition  of class Curve_point, which  represents a vertex
 * of a discrete geodesic (an  open, simple polygonal line) during the
 * execution of the discrete geodesic finder algorithm first step. The
 * geodesic vertex can be a vertex of the mesh over which the geodesic
 * is defined,  or it  may belong to  the interior  of a mesh  edge or
 * face.
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

#ifndef ADG_CURVE_POINT_HPP
#define ADG_CURVE_POINT_HPP

#include "geodesic_point.hpp"          // Geodesic_point
#include "geodesic_point_vertex.hpp"   // Geodesic_point_vertex
#include "geodesic_point_edge.hpp"     // Geodesic_point_edge
#include "geodesic_point_face.hpp"     // Geodesic_point_face


#include <cmath>                       // HUGE_VAL
#include <vector>                      // std::vector


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
   * \class Curve_point
   *
   * \brief This class represents a  vertex of a discrete geodesic (an
   * open, simple polygonal line) during the execution of the discrete
   * geodesic finder algorithm first step.  The geodesic vertex can be
   * a vertex  of the mesh over  which the geodesic is  defined, or it
   * may belong to the interior of a mesh edge or face.
   */
  template< typename Mesh >
  class Curve_point {
  public:
    // ---------------------------------------------------------------
    //
    // Type name definitions.
    //
    // ---------------------------------------------------------------

    /**
     * \typedef SIDE
     *
     * \brief Definition of a  type name for classifying this geodesic
     * vertex with respect to a shorter version of the geodesic (i.e.,
     * the vertex  may be  either on the  shorter geodesic, or  on its
     * left or right  side. If the vertex has  not been classified yet
     * or  it  lies on  an  edge  of the  mesh,  it  is classified  as
     * NON_INITIALIZED.
     */
    typedef enum {
      NON_INITIALIZED = 0 ,
      NO_SIDE = 1 ,
      LEFT_SIDE = 2 ,
      RIGHT_SIDE = 3
    } 
    SIDE ;
    

    /**
     * \typedef Vertex
     *
     * \brief Definition of a type name for the mesh vertices.
     */
    typedef typename Mesh::Vertex Vertex ;


    /**
     * \typedef Edge
     *
     * \brief Definition of a type name for the mesh edges.
     */
    typedef typename Mesh::Edge Edge ;


    /**
     * \typedef Face
     *
     * \brief Definition of a type name for the face edges.
     */
    typedef typename Mesh::Face Face ;


    /**
     * \typedef GP
     *
     * \brief Definition of a type name for a geodesic point.
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


    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn Curve_point()
     *
     * \brief Creates an  instance of this class.
     */
    Curve_point() 
    {
      _error = HUGE_VAL ;
      _side = NON_INITIALIZED ;
      _point = 0 ;

      return ;
    }


    /**
     * \fn Curve_point( GP* p )
     *
     * \brief Creates an  instance of this class.
     *
     * \param p Pointer to a geodesic  vertex (it can be either a mesh
     * vertex or  a point in  the interior of  an edge or face  of the
     * mesh.
     */
    Curve_point( GP* p ) 
    {
#ifdef DEBUGMODE
      assert( p != 0 ) ;
      assert( p->on_vertex() || p->on_edge() || p->on_face() ) ;
#endif
      _error = HUGE_VAL ;
      _side = NON_INITIALIZED ;
      _point = p ;

      return ;
    }


   /**
     * \fn Curve_point( const Curve_point& cp )
     *
     * \brief Creates an instance  of this class form another instance
     * of this class.
     *
     * \param cp A reference to an instance of this class.
     */
    Curve_point( const Curve_point& cp )
    {
      _error = cp._error ;
      _side = cp._side ;
      _point = cp._point ;

      return ;
    }


    /**
     * \fn ~Curve_point()
     *
     * \brief Destroys an instance of this class.
     */    
    ~Curve_point()
    {}


    /**
     * \fn inline SIDE get_side() const
     *
     * \brief Returns  the \a SIDE value associated  with this vertex.
     * This value indicates the type  of this vertex with respect to a
     * shorter  geodesic being  computed by  the geodesic  finder path
     * correction algorithm.
     *
     * \return The \a SIDE value associated with this vertex.
     */
    inline SIDE get_side() const
    {
      return _side ;
    }


    /**
     * \fn inline void set_side( SIDE side )
     *
     * \brief  Assigns a  \a SIDE  value  to this  vertex. This  value
     * indicates the  type of  this vertex with  respect to  a shorter
     * geodesic being computed by  the geodesic finder path correction
     * algorithm.
     *
     * \param side A \a SIDE value to be assigned with this vertex.
     */
    inline void set_side( SIDE side )
    {
      _side = side ;
    }


    /**
     * \fn inline double get_error() const
     *
     * \brief  Returns the  error associated  with this  vertex.  This
     * value is a measure of how straight the geodesic currently being
     * computed by the path  correction algorithm is in a neighborhood
     * of this vertex.
     *
     * \return The error associated with this vertex.
     */
    inline double get_error() const
    {
      return _error ;
    }


    /**
     * \fn inline void set_error( double error )
     *
     * \brief Assigns an error value with this vertex. This value is a
     * measure of  how straight the geodesic  currently being computed
     * by the path  correction algorithm is in a  neighborhood of this
     * vertex.
     *
     * \param error The value to  the assigned to the error associated
     * with this vertex.
     */
    inline void set_error( double error ) 
    {
      _error = error ;
    }


    /**
     * \fn inline GP* get_point() const
     *
     * \brief Returns a pointer to the actual geodesic point.
     *
     * \return A pointer to the actual geodesic point.
     */
    inline GP* get_point() const
    {
      return _point ;
    }


    /**
     * \fn inline void set_point( GP* point )
     *
     * \brief Assigns the address of  a geodesic point to the geodesic
     * point pointer of this instance.
     *
     * \param point A pointer to a geodesic point.
     */
    inline void set_point( GP* point )
    {
#ifdef DEBUGMODE
      assert( point != 0 ) ;
#endif

      _point = point ;
    }


    /**
     * \fn inline bool on_vertex() const
     *
     * \brief Determines  whether this geodesic vertex is  a vertex of
     * the mesh over which the geodesic is defined.
     *
     * \return  The logic  value true  if  this geodesic  vertex is  a
     * vertex of the mesh over  which the geodesic is defined, and the
     * logic value false otherwise.
     */    
    inline bool on_vertex() const
    {
      assert( _point != 0 ) ;

      return _point->on_vertex() ;
    }


    /**
     * \fn inline bool on_edge() const
     *
     * \brief Determines  whether this geodesic vertex  belongs to the
     * interior  of an edge  of the  mesh over  which the  geodesic is
     * defined.
     *
     * \return The logic value true  if this geodesic vertex is in the
     * interior  of an edge  of the  mesh over  which the  geodesic is
     * defined, and the logic value false otherwise.
     */  
    inline bool on_edge() const
    {
      assert( _point != 0 ) ;

      return _point->on_edge() ;
    }


    /**
     * \fn inline bool on_face() const
     *
     * \brief Determines  whether this geodesic vertex  belongs to the
     * interior  of a  face of  the mesh  over which  the  geodesic is
     * defined.
     *
     * \return The logic value true  if this geodesic vertex is in the
     * interior  of a  face of  the mesh  over which  the  geodesic is
     * defined, and the logic value false otherwise.
     */  
    inline bool on_face() const
    {
      assert( _point != 0 ) ;

      return _point->on_face() ;
    }


    /**
     * \fn inline Vertex* get_vertex() const
     *
     * \brief  Returns  the  mesh  vertex  that  coincides  with  this
     * geodesic point.  If this geodesic  point is not a  mesh vertex,
     * throws an exception.
     *
     * \return A pointer  to the mesh vertex that  corresponds to this
     * geodesic point.
     */
    inline Vertex* get_vertex() const
    {
#ifdef DEBUGMODE
      assert( on_vertex() ) ;
#endif

      GPV* gpv = dynamic_cast< GPV* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpv != 0 ) ;
#endif

      return gpv->get_vertex() ;
    }


    /**
     * \fn inline Edge* get_edge() const
     *
     * \brief Returns the mesh  edge that contains this geodesic point
     * in its interior.  If this geodesic point is not in a mesh edge,
     * throws an exception.
     *
     * \return A pointer to the  mesh edge that contains this geodesic
     * vertex.
     */
    inline Edge* get_edge() const
    {
#ifdef DEBUGMODE
      assert( on_edge() ) ;
#endif

      GPE* gpe = dynamic_cast< GPE* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpe != 0 ) ;
#endif

      return gpe->get_edge() ;
    }


    /**
     * \fn inline double get_t() const
     *
     * \brief  Returns  the ratio  that  indicates  how  far away  the
     * geodesic vertex  is from  the origin mesh  vertex of  the first
     * half-edge of the edge on which this geodesic vertex lies.
     *
     * \return  The ratio  that indicates  how far  away  the geodesic
     * vertex is from the origin mesh vertex of the first half-edge of
     * the edge this vertex lies on.
     */
    inline double get_t() const
    {
#ifdef DEBUGMODE
      assert( on_edge() ) ;
#endif

      GPE* gpe = dynamic_cast< GPE* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpe != 0 ) ;
#endif

      return gpe->get_t() ;
    }


    /**
     * \fn inline Face* get_face() const
     *
     * \brief Returns the mesh  face that contains this geodesic point
     * in its interior.  If this geodesic point is not in a mesh face,
     * throws an exception.
     *
     * \return A pointer to the  mesh face that contains this geodesic
     * vertex.
     */
    inline Face* get_face() const
    {
#ifdef DEBUGMODE
      assert( on_face() ) ;
#endif

      GPF* gpf = dynamic_cast< GPF* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpf != 0 ) ;
#endif

      return gpf->get_face() ;
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
#ifdef DEBUGMODE
      assert( on_face() ) ;
#endif

      GPF* gpf = dynamic_cast< GPF* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpf != 0 ) ;
#endif

      return gpf->get_u() ;
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
#ifdef DEBUGMODE
      assert( on_face() ) ;
#endif

      GPF* gpf = dynamic_cast< GPF* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpf != 0 ) ;
#endif

      return gpf->get_v() ;
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
      assert( on_vertex() ) ;
      assert( v != 0 ) ;
#endif

      GPV* gpv = dynamic_cast< GPV* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpv != 0 ) ;
#endif

      gpv->set_vertex( v ) ;

      return ;
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
      assert( on_edge() ) ;
      assert( 
	     ( e != 0 ) && 
	     ( fabs( t ) > 0 ) &&  ( fabs( t ) < 1 ) 
	    );
#endif

      GPE* gpe = dynamic_cast< GPE* >( _point ) ;

#ifdef DEBUGMODE
      assert( gpe != 0 ) ;
#endif

      gpe->set_edge( e , t ) ;

      return ;
    }

private:

    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    double _error ;         ///< Error value of this vertex. 

    SIDE _side ;            ///< Side value of this vertex.

    GP* _point ;            ///< Pointer to the actual geodesic vertex.

  } ;

}

/** @} */ //end of group class.

#endif	/* ADG_CURVE_POINT_HPP */

