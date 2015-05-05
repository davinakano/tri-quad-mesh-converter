/** 
 * \file floaterpar.hpp
 *  
 * \brief Definition  of class  FloaterPar, which implements  a method
 * for  parametrizing  a  disk-like  surface  patch  in  \f$R^3\f$  to
 * \f$R^2\f$.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * mfsiqueira at gmail (dot) com
 *
 * \version 1.0
 * \date March 2011
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef FLOATERPAR_HPP
#define FLOATERPAR_HPP

#include "meshpatch/patch.hpp"                          // meshpatch::Patch
#include "meshpatch/patchedge.hpp"                      // meshpatch::PatchEdge
#include "meshpatch/patchvertex.hpp"                    // meshpatch::PatchVertex

#include "floater/floaterpar_vertex_attributes.hpp"   // FloaterPar_vertex_attributes

#include <vector>                             // std::vector
#include <cmath>                              // acos

/** 
 * \defgroup FLOATERPARNameSpace Namespace floaterpar.
 * @{
 */

/**
 * \namespace floaterpar
 *
 * \brief The namespace floaterpar  contains all classes related to an
 * implementation  of  the mesh  parametrization  method developed  by
 * Floater.
 */

namespace floaterpar {

  /**
   * \class FloaterPar
   *
   * \brief  This  class  contains  a  function  for  parametrizing  a
   * triangular surface  patch in \f$R^3\f$ to  \f$R^2\f$ by Floater's
   * method.
   */
  class FloaterPar {
  public:

    // ---------------------------------------------------------------
    //
    // Type definitions.
    //
    // ---------------------------------------------------------------

    /**
     * \typedef PMesh
     *
     * \brief  Definition of a  type name  for the  triangular surface
     * patch.
     */
    typedef meshpatch::Patch< FloaterPar_vertex_attributes > PMesh ;


    /**
     * \typedef PVertex
     *
     * \brief Definition  of a type name for  triangular surface patch
     * vertices.
     */
    typedef meshpatch::PatchVertex< FloaterPar_vertex_attributes > PVertex ;


    /**
     * \typedef PEdge
     *
     * \brief Definition  of a type name for  triangular surface patch
     * edges.
     */
    typedef meshpatch::PatchEdge< FloaterPar_vertex_attributes > PEdge ;


    // ---------------------------------------------------------------
    //
    // Public methods.
    //
    // ---------------------------------------------------------------

    /**
     * \fn FloaterPar( unsigned nv , double* vset , unsigned nf , unsigned* fset , double* tv )
     *
     * \brief Creates an  instance of this class.
     *
     * \param nv The number of vertices of the surface patch.
     * \param  vset The  Cartesian  coordinates of  the surface  patch
     * vertices.
     * \param nf number of faces of the surface patch.
     * \param fset The indices of the vertices of each face.
     * \param tv an array with the t-values of the patch vertices. 
     */    
    FloaterPar(
	       unsigned   nv,
	       double*    vset,
	       unsigned   nf,  
	       unsigned*  fset,
	       double*    tv
	      )
    :
      _MYPI( acos( -1 ) ) ,
      _is_done( false )
    {
      //
      // Builds the surface patch to be parametrized.
      //
      _patch = new PMesh( nv , nf , fset , vset ) ;

      //
      // Initialize vertex attributes.
      //
      initialize_vertex_attributes( tv ) ;

      return ;
    }


    /**
     * \fn FloaterPar( PMesh* patch )
     *
     * \brief Creates an  instance of this class.
     *
     * \param patch A pointer to a surface patch to be parametrized by
     * Floater's method.
     */    
    FloaterPar( PMesh* patch ) 
      : 
      _MYPI( acos( -1 ) ) , 
      _is_done( false )
    {
      //
      // Initialize the patch pointer.
      //
      _patch = new PMesh( *patch ) ;

      //
      // Initialize vertex attributes.
      //
      for ( unsigned i = 0 ; i < _patch->get_number_of_vertices() ; ++i ) {
	/*
	 * Set dummy coordinate values.
	 */
	_patch->get_vertex_attributes( i ).set_ucoord( HUGE_VAL ) ;
	_patch->get_vertex_attributes( i ).set_vcoord( HUGE_VAL ) ;
      }

      return ;
    }


    /**
     * \fn ~FloaterPar()
     *
     * \brief Destroys this object.
     */
    ~FloaterPar()
    {
      if ( patch() != 0 ) {
	delete patch() ;
      }

      return ;
    }


    /** 
     * \fn void parametrize() 
     *
     * \brief  computes a  parametrization of  the  associated surface
     * patch using Floater's method as  described in the paper "M.  S.
     * Floater,  Parametrization and  smooth approximation  of surface
     * triangulations,  Computer  Aided  Geometric Design  14  (1997),
     * 231-250."
     */
    void parametrize() ;


    /**
     * \fn inline PMesh* patch() const 
     *
     * \brief Returns  a pointer to the surface  patch associated with
     * this object.
     *
     * @\return A  pointer to the  surface patch associated  with this
     * object.
     */
    inline PMesh* patch() const
    { 
      return _patch ; 
    }


    /**
     * \fn inline unsigned get_number_of_vertices() const 
     *
     * \brief Returns the number of vertices of the associated surface
     * patch.
     *
     * @\return  The  number of  vertices  of  the associated  surface
     * patch.
     */
    inline unsigned get_number_of_vertices() const
    { 
      return patch()->get_number_of_vertices() ; 
    }


    /**
     * \fn inline unsigned get_number_of_edges() const 
     *
     * \brief Returns  the number of  edges of the  associated surface
     * patch.
     *
     * @\return The number of edges of the associated surface patch.
     */
    inline unsigned get_number_of_edges() const
    { 
      return patch()->get_number_of_edges() ; 
    }


    /**
     * \fn inline unsigned get_number_of_faces() const 
     *
     * \brief Returns  the number of  faces of the  associated surface
     * patch.
     *
     * @\return The number of faces of the associated surface patch.
     */
    inline unsigned get_number_of_faces() const
    { 
      return patch()->get_number_of_faces() ; 
    }


    /**
     *  \fn inline double get_tval( unsigned i ) const
     *
     * \brief Returns  the t value of  the i-th vertex  of the surface
     * patch.
     *
     * \param i The index of the i-th vertex of the associated surface
     * patch.
     *
     * \return The t value of the i-th vertex of the surface patch.
     */
    inline double get_tval( unsigned i ) const
    {
      return patch()->get_vertex_attributes( i ).get_tval() ; 
    }


    /**
     * \fn inline double get_ucoord( unsigned i ) const
     *
     * \brief Returns the first barycentric (planar) coordinate of the
     * i-th vertex of the associated surface patch.
     *
     * \param i The index of the i-th vertex of the associated surface
     * patch.
     *
     * \return The  first barycentric (planar) coordinate  of the i-th
     * vertex of the surface patch.
     */
    inline double get_ucoord( unsigned i ) const
    {
      return patch()->get_vertex_attributes( i ).get_ucoord() ;
    }


    /**
     *  \fn inline double get_vcoord( unsigned i ) const
     *
     * \brief  Returns the second  barycentric (planar)  coordinate of
     * the i-th vertex of the surface patch.
     *
     * \param i The index of the i-th vertex of the associated surface
     * patch.
     *
     * \return The second barycentric  (planar) coordinate of the i-th
     * vertex of the surface patch.
     */
    inline double get_vcoord( unsigned i ) const
    {
      return patch()->get_vertex_attributes( i ).get_vcoord() ;
    }


    /**
     * \fn void sample( unsigned nv , double* param_pts , double*& image_pts )
     *
     * \brief Sample  the parametrization associated  with this object
     * to create a set of 3d points, which are the image of the sample
     * points.
     *
     * \param nv The number of sample points.
     * \param  param_pts  The  Cartesian  coordinates  of  the  sample
     * (parameter) points.
     * \param image_pts The Cartesian  coordinates of the image of the
     * parameter points.
     */
    void sample( unsigned nv , double* param_pts , 
      double*& image_pts ) ;


  private:
    // ---------------------------------------------------------------
    //
    // Private data types
    //
    // ---------------------------------------------------------------

    /**
     * \struct BoundingBox
     *
     * \brief A struct representing a bounding box in 2D.
     */
    struct BoundingBox {
      double _xmin ;    ///< The X coordinate of the leftmost-lower point of the bounding box.
      double _ymin ;    ///< The Y coordinate of the leftmost-lower point of the bounding box.
      double _xmax ;    ///< The X coordinate of the rightmost-upper point of the bounding box.
      double _ymax ;    ///< The Y coordinate of the rightmost-upper point of the bounding box.
    } ;


    // ---------------------------------------------------------------
    //
    // Private member functions
    //
    // ---------------------------------------------------------------

    /**
     * \fn void initialize_vertex_attributes( double* tv ) 
     *
     * \brief A method for creating vertex attributes.
     *
     * \param tv An array with  the t values associated with the patch
     * vertices.
     */
    void initialize_vertex_attributes( double* tv ) ;


    /**
     * \fn void separate_vertices( std::vector< unsigned >& iv , std::vector< unsigned >& bv )
     *
     * \brief  Separates  inner  from  boundary  vertices.   During  the
     * separation process, the consistency of the t-values are checked.
     *
     * \param iv A reference for an array with the inner vertices.
     * \param bv A reference for an array with the boundary vertices.
     */
    void separate_vertices( std::vector< unsigned >& iv ,
      std::vector< unsigned >& bv ) ;


    /** 
     * \fn void comp_bv_pcoords( const std::vector< unsigned >& bv ) const
     *
     * \brief  Computes  the   parameter  coordinates  of  the  boundary
     * vertices.
     *
     * \param bv An array with the boundary vertices.
     */  
    void comp_bv_pcoords( const std::vector< unsigned >& bv ) const ;


    /** 
     * \fn void comp_mv_coord( unsigned vi , const std::vector< unsigned >& star , const std::vector< unsigned >& vtab ,  double*& A , double& U , double& V )
     *
     * \brief Computes  the mean value  coordinates of a  given vertex
     * with  respect  to  the  vertices   in  its  star.  This  is  an
     * implementation of the method in
     *
     * M.  S.    Floater,  Mean  value   coordinates,  Computer  Aided
     * Geometric Design 20 (2003), 19-27.
     *
     * \param vi The index of an inner patch vertex.
     * \param star A  list with the vertices in the  star of the given
     * inner vertex.
     * \param vtab  A lookup table used  to speed up  access to vertex
     * ID.
     * \param A  An array to store  the mean value  coordinates of the
     * inner vertices of the star.
     * \param U A  variable to store the convex sum  of the mean value
     * coordinates times the U  coordinate of the boundary vertices of
     * the star.
     * \param V A  variable to store the convex sum  of the mean value
     * coordinates times the U  coordinate of the boundary vertices of
     * the star.
     */
    void comp_mv_coord( unsigned vi ,
      const std::vector< unsigned >& star ,
      const std::vector< unsigned >& vtab , double*& A , double& U ,
      double& V ) ; 


    /**
     * \fn bool is_inner_vertex( unsigned i ) const
     *
     * \brief Decides whether a given patch vertex is an inner vertex.
     *
     * \param i The index of a patch vertex.
     *
     * \return  The logic value  true if  the the  vertex is  an inner
     * vertex and the logic value false otherwise.
     */
    bool is_inner_vertex( unsigned i ) const ;


    /**
     * \fn void compute_bounding_box( std::vector< BoundingBox >& bbset )
     *
     * \brief Computes a bounding box for each parametrization triangle.
     *
     * \param bbset An array of bounding boxes (one for each triangle).
     */
    void compute_bounding_box( std::vector< BoundingBox >& bbset ) ;


    /**
     * \fn void compute_image_points( unsigned nv , double* param_pts , const std::vector< BoundingBox >& bbset , std::vector< double >& pset ) 
     *
     * \brief Compute  the image  points of a  given set  of parameter
     * points.
     *
     * \param nv The number of parameter points.
     * \param param_pts Coordinates of the parameter points.
     * \param bbset An array of bounding boxes (one for each triangle).
     * \param pset A reference to an array of image points.
     */
    void compute_image_points( unsigned nv , double* param_pts ,
      const std::vector< BoundingBox >& bbset ,
      std::vector< double >& pset ) ;


    /**
     * \fn inline bool is_point_in_bounding_box( double x , double y , const BoundingBox& bb ) const
     *
     * \brief Decides whether a point belongs to a given bounding box.
     *
     * \param x The first Cartesian coordinate of the point.
     * \param y The second Cartesian coordinate of the point.
     * \param bb A bounding box.
     */
    inline bool is_point_in_bounding_box(
					 double x ,
					 double y ,
					 const BoundingBox& bb
					)
      const
    {
      return ( x >= bb._xmin ) &&
	     ( x <= bb._xmax ) && 
	     ( y >= bb._ymin ) && 
             ( y <= bb._ymax ) ;
    }


    /**
     * \fn void compute_barycentric_coordinates( double xp , double yp , double x0 , double y0 , double x1 , double y1 , double x2 , double y2 , double& u , double& v , double& w ) const
     * 
     * \brief Computes  the barycentric  coordinates of a  given point
     * (in Cartesian  coordinates) with  respect to a  given reference
     * triangle.
     *
     * \param xp First Cartesian coordinate of the point.
     * \param yp Second Cartesian coordinate of the point.
     * \param x0 First Cartesian coordinate of the first vertex of the
     * reference triangle.
     * \param y0  Second Cartesian coordinate  of the first  vertex of
     * the reference triangle.
     * \param x1  First Cartesian coordinate  of the second  vertex of
     * the reference triangle.
     * \param y1  Second Cartesian coordinate of the  second vertex of
     * the reference triangle.
     * \param x2 First Cartesian coordinate of the third vertex of the
     * reference triangle.
     * \param y2  Second Cartesian coordinate  of the third  vertex of
     * the reference triangle.
     * \param u First barycentric coordinate of the point.
     * \param v Second barycentric coordinate of the point.
     * \param w Third barycentric coordinate of the point.
     *
     */
    void compute_barycentric_coordinates( double xp , double yp , 
      double x0 , double y0 , double x1 , double y1 , double x2 ,
      double y2 , double& u , double& v , double& w ) const ; 


    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    const double _MYPI ;   ///< The constant PI.

    bool _is_done ;        ///< A flag to indicate whether a parametrization has been computed. 

    PMesh* _patch ;        ///< A pointer to a triangular surface patch.

  } ;

}

/** @} */ //end of group class.

#endif	/* FLOATERPAR_HPP */
