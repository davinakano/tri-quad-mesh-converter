/** 
 * \file floaterpar_vertex_attributes.hpp
 *  
 * \brief  Definition  of  class  FloaterPar_vertex_attributes,  which
 * represents  the set of  vertex attributes  of a  triangular surface
 * patch   in  \f$R^3\f$   parametrized  to   a  plane   by  Floater's
 * parametrization.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
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

#ifndef FLOATERPAR_VERTEX_ATTRIBUTES_HPP
#define	FLOATERPAR_VERTEX_ATTRIBUTES_HPP

#include <cmath>    // HUGE_VAL

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
   * \class FloaterPar_vertex_attributes
   *
   * \brief This  class represents the  set of vertex attributes  of a
   * triangular surface patch in  \f$R^3\f$ parametrized to a plane by
   * Floater's method.
   */
  class FloaterPar_vertex_attributes {
  public:
    //
    // Public methods
    //
    // ---------------------------------------------------------------


    /**
     * \fn FloaterPar_vertex_attributes()
     *
     * \brief Creates an instance of this class.
     */
    FloaterPar_vertex_attributes()
    {
      set_tval( -1 ) ;
      set_ucoord( -HUGE_VAL ) ;
      set_vcoord( -HUGE_VAL ) ;
    }


    /**
     * \fn FloaterPar_vertex_attributes( double t , double u , double v )
     *
     * \brief Creates an instance of this class.
     *
     * \param t The t-value associated  with the vertex that owns this
     * attribute.
     * \param u  First barycentric coordinate of the  vertex that owns
     * this attribute.
     * \param v Second barycentric  coordinate of the vertex that owns
     * this attribute.
     */
    FloaterPar_vertex_attributes( double t , double u , double v )
    {
      set_tval( t ) ;
      set_ucoord( u ) ;
      set_vcoord( v ) ;
    }


    /**
     * \fn ~FloaterPar_vertex_attributes()
     *
     * \brief Destroys this object.
     */
    ~FloaterPar_vertex_attributes()
    {}


    /**
     * \fn inline double get_tval() const
     *
     * \brief  Returns  the  t-value  of  the vertex  that  owns  this
     * attribute,  i.e.,  the  parameter  that indicates  whether  the
     * vertex that owns  this attribute is a boundary  vertex.  If so,
     * the parameter  also indicates the position of  the vertex along
     * the patch boundary.
     *
     * \return The  t-value associated with the vertex  that owns this
     * attribute.
     */
    inline double get_tval() const
    { 
      return _t ; 
    }


    /**
     * \fn inline double get_ucoord() const
     *
     * \brief Returns  the first barycentric coordinate  of the vertex
     * that owns this attribute.
     *
     * \return  The first  barycentric coordinate  of the  vertex that
     * owns this attribute.
     */
    inline double get_ucoord() const
    { 
      return _u ;
    }


    /**
     * \fn inline double get_vcoord() const
     *
     * \brief Returns the second  barycentric coordinate of the vertex
     * that owns this attribute.
     * 
     * \return The  second barycentric  coordinate of the  vertex that
     * owns this attribute.
     */
    inline double get_vcoord() const
    { 
      return _v ; 
    }


    /**
     * \fn inline void set_tval( double t )
     *
     * \brief Sets the t-value of the vertex that owns this attribute.
     *
     * \param t The  t-value to be assigned with  the vertex that owns
     * this attribute.
     */
    inline void set_tval( double t )
    { 
      _t = t ; 
    }


    /**
     * \fn inline void set_ucoord( double u )
     *
     * \brief Sets  the value of  the first barycentric  coordinate of
     * the vertex that owns this attribute.
     *
     * \param u  The first barycentric  coordinate of the  vertex that
     * owns this attribute.
     */
    inline void set_ucoord( double u )
    {
      _u = u ;
    }


    /**
     * \fn inline void set_vcoord( double v )
     *
     * \brief Sets  the value of the second  barycentric coordinate of
     * the vertex that owns this attribute.
     *
     * \param v  The second barycentric coordinate of  the vertex that
     * owns this attribute.
     */
    inline void set_vcoord( double v )
    { 
      _v = v ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    double _t ;  ///< T-value associated with the vertex that owns this attribute.

    double _u ;  ///< First barycentric coordinate of the vertex that owns this attribute.

    double _v ;  ///< Second barycentric coordinate of the vertex that owns this attribute.

  } ;

}

/** @} */ //end of group class.

#endif   // FLOATERPAR_VERTEX_ATTRIBUTES_HPP

