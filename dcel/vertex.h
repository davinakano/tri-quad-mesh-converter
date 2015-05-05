/** 
 * \file vertex.h
 *  
 * \brief Definition  of the class  Vertex, which represents  a vertex
 * (i.e., a  point in R3)  from a surface  mesh represented by  a DCEL
 * data structure.
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

#ifndef VERTEX_H
#define VERTEX_H


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
   * Forward definition of classes.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib
           >
  class Halfedge ;


  /**
   * \class Vertex
   *
   * \brief This class represents a  vertex (i.e., a point in R3) from
   * a triangle surface mesh represented by the DCEL data structure.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib
           >
  class Vertex {
  public:
    // ---------------------------------------------------------------
    //
    // Type definitions
    //
    // ---------------------------------------------------------------

    /**
     * \typedef Halfedge
     *
     * \brief  Defines  Halfedge   as  an  alias  for  dcel::Halfedge<
     * VAttrib, FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Halfedge< VAttrib, FAttrib , EAttrib , HAttrib > 
      Halfedge ;


    // ---------------------------------------------------------------
    //
    // Public methods.
    //
    // ---------------------------------------------------------------

    /**
     * \fn Vertex( double x , double y , double z , Halfedge* h )
     *
     * \brief Creates an instance of this class.
     *
     * \param x The first Cartesian coordinate of the vertex location.
     * \param  y  The  second   Cartesian  coordinate  of  the  vertex
     * location.
     * \param z The third Cartesian coordinate of the vertex location.
     * \param h A pointer to a half-edge with origin at this vertex.
     *
     */
    Vertex(
	   double x ,
	   double y ,
	   double z ,
	   Halfedge* h
	  )
    {
      set_x_coord( x ) ;
      set_y_coord( y ) ;
      set_z_coord( z ) ;
      set_halfedge( h ) ;
    }


    /**
     * \fn ~Vertex()
     *
     * \brief Destroys an instance of this class.
     *
     */
    ~Vertex() 
    {
      set_halfedge( 0 ) ;
    }


    /**
     * \fn double x() const
     *
     * \brief Returns  the first Cartesian coordinate  of the location
     * of this vertex.
     *
     * \return The first Cartesian  coordinate of the location of this
     * vertex.
     */
    double x() const
    {
      return _x ; 
    }


    /**
     * \fn void set_x_coord( double x ) 
     *
     * \brief Assigns a value to the first Cartesian coordinate of the
     * location to this vertex.
     *
     * \param  x The  value  to  be assigned  to  the first  Cartesian
     * coordinate of the location to this vertex.
     */
    void set_x_coord( double x )
    {
      _x = x ;
    }


    /**
     * \fn double y() const
     *
     * \brief Returns the second  Cartesian coordinate of the location
     * of this vertex.
     *
     * \return The second Cartesian coordinate of the location of this
     * vertex.
     */
    double y() const
    {
      return _y ; 
    }


    /**
     * \fn void set_y_coord( double y ) 
     *
     * \brief Assigns  a value to  the second Cartesian  coordinate of
     * the location to this vertex.
     *
     * \param  y The  value to  be  assigned to  the second  Cartesian
     * coordinate of the location to this vertex.
     */
    void set_y_coord( double y )
    {
      _y = y ;
    }


    /**
     * \fn double z() const
     *
     * \brief Returns  the third Cartesian coordinate  of the location
     * of this vertex.
     *
     * \return The third Cartesian  coordinate of the location of this
     * vertex.
     */
    double z() const
    {
      return _z ;
    }


    /**
     * \fn void set_z_coord( double z ) 
     *
     * \brief Assigns a value to the third Cartesian coordinate of the
     * location to this vertex.
     *
     * \param  z The  value  to  be assigned  to  the third  Cartesian
     * coordinate of the location to this vertex.
     */
    void set_z_coord( double z )
    {
      _z = z ;
    }


    /**
     * \fn Halfedge* get_halfedge() const 
     *
     * \brief Returns a pointer to  a half-edge whose origin vertex is
     * this vertex.
     *
     * \return A  pointer to a  half-edge whose origin vertex  is this
     * vertex.
     */
    Halfedge* get_halfedge() const
    { 
      return _halfedge ; 
    }


    /**
     * \fn void set_halfedge( Halfedge* h ) 
     *
     * \brief Assigns an address to the halfedge pointer of this vertex.
     *
     * \param h The address of a halfedge with origin at this vertex.
     */
    void set_halfedge( Halfedge* h )
    {
      _halfedge = h ;
    }


    /**
     * \fn VAttrib& get_attributes()
     *
     * \brief  Returns  the set  of  attributes  associated with  this
     * vertex.
     *
     * \return A reference to the set of attributes associated with this
     * vertex.
     */
    VAttrib& get_attributes()
    {
      return _attributes ; 
    }

 

  private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    double _x ;   ///< First Cartesian coordinate of the location of this vertex.
    double _y ;   ///< Second Cartesian coordinate of the location of this vertex.
    double _z ;   ///< Third Cartesian coordinate of the location of this vertex.

    Halfedge* _halfedge ;  ///< Pointer to a half-edge with origin at this vertex.

    VAttrib _attributes ;  ///< Set of attributes associated with this vertex. 

  } ;

}

/** @} */ //end of group class.

#endif  // VERTEX_H
