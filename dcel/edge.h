/** 
 * \file edge.h
 *  
 * \brief  Definition of  the  class Edge,  which  represents an  edge
 * (i.e., a  line segment) from a  surface mesh represented  by a DCEL
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

#ifndef EDGE_H
#define EDGE_H


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
   * \class Edge
   *
   * \brief This class represents an  edge (i.e., a line segment) from
   * a triangle surface mesh represented by the DCEL data structure.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib
           >
  class Edge {
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
    // Public methods
    //
    // ---------------------------------------------------------------

    /** 
     * \fn Edge( Halfedge* h1 , Halfedge* h2 )
     *
     * \brief Creates an instance of this class.
     *
     * \param h1 A pointer to one half-edge of this edge.
     * \param h2 A pointer to the other half-edge of this edge.
     */
    Edge(
	 Halfedge* h1 ,
	 Halfedge* h2
	)
    {
      set_first_halfedge( h1 ) ;
      set_second_halfedge( h2 ) ;
    }


    /**
     * \fn ~Edge()
     *
     * \brief Destroys an instance of this class.
     */
    ~Edge()
    {
      set_first_halfedge( 0 ) ;
      set_second_halfedge( 0 ) ;
    }


    /**
     * \fn Halfedge* get_first_halfedge() const 
     *
     * \brief Returns a pointer to the first half-edge of this edge.
     *
     * \return A pointer to the first half-edge of this edge.
     */    
    Halfedge* get_first_halfedge() const
    {
      return _h1 ; 
    }


    /**
     * \fn Halfedge* get_second_halfedge() const 
     *
     * \brief Returns a pointer to the second half-edge of this edge.
     *
     * \return A pointer to the second half-edge of this edge.
     */    
    Halfedge* get_second_halfedge() const
    {
      return _h2 ; 
    }


    /**
     * \fn bool is_boundary_edge() const 
     *
     * \brief Returns the logic value true if this edge belongs to the
     * boundary of the surface and the logic value false otherwise.
     *
     * \return  The logic  value  true  if this  edge  belongs to  the
     * boundary of the surface and the logic value false otherwise.
     */
    bool is_boundary_edge() const
    {
      return ( _h1 == 0 ) || ( _h2 == 0 ) ;
    }


    /**
     * \fn void set_first_halfedge( Halfedge* h )
     *
     * \brief  Assigns a  given address  to the  pointer to  the first
     * half-edge of this edge.
     *
     * \param h The address of the first half-edge of this edge.
     */
    void set_first_halfedge( Halfedge* h ) 
    {
      _h1 = h ; 
    }


    /**
     * \fn void set_second_halfedge( Halfedge* h )
     *
     * \brief Assigns  a given  address to the  pointer to  the second
     * half-edge of this edge.
     *
     * \param h The address of the second half-edge of this edge.
     */
    void set_second_halfedge( Halfedge* h ) 
    {
      _h2 = h ; 
    }


    /**
     * \fn EAttrib& get_attributes()
     *
     * \brief Returns the set of attributes associated with this edge.
     *
     * \return A reference to the set of attributes associated with this
     * edge.
     */
    EAttrib& get_attributes()
    {
      return _attributes ; 
    }

 

  private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    Halfedge* _h1 ;  ///< Pointer to the first half-edge of this edge.

    Halfedge* _h2 ;  ///< Pointer to the second half-edge of this edge.

    EAttrib _attributes ;  ///< Set of attributes associated with this edge. 

  } ;

}

/** @} */ //end of group class.

#endif   // EDGE_H
