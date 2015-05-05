/** 
 * \file halfedge.h
 *  
 * \brief  Definition  of  the  class  Halfedge,  which  represents  a
 * half-edge from a triangle surface mesh represented by the DCEL data
 * structure.
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

#ifndef HALFEDGE_H
#define HALFEDGE_H

#include "vertex.h"   // Vertex
#include "edge.h"     // Edge


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
  class Face ;


  /**
   * \class Halfedge
   *
   * \brief This class represents  a half-edge (i.e., an oriented line
   * segment)  from a triangle  surface mesh  represented by  the DCEL
   * data structure.
   */
  template < 
             typename VAttrib ,
             typename FAttrib ,
             typename EAttrib ,
             typename HAttrib
           >
  class Halfedge {
  public:
    // ---------------------------------------------------------------
    //
    // Type definitions
    //
    // ---------------------------------------------------------------

    /**
     * \typedef Vertex
     *
     * \brief Defines  Vertex as  an alias for  dcel::Vertex< VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Vertex< VAttrib, FAttrib , EAttrib , HAttrib > Vertex ;

    /**
     * \typedef Edge
     *
     * \brief  Defines  Edge  as  an alias  for  dcel::Edge<  VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Edge< VAttrib, FAttrib , EAttrib , HAttrib > Edge ;

    /**
     * \typedef Face
     *
     * \brief  Defines  Face  as  an alias  for  dcel::Face<  VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef dcel::Face< VAttrib, FAttrib , EAttrib , HAttrib > Face ;


    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /** 
     * \fn Halfedge( Vertex* vertex , Edge* edge , Face* face , Halfedge* next , Halfedge* prev )
     *
     * \brief Creates an instance of this class.
     *
     * \param vertex A pointer to the origin vertex of this halfedge.
     * \param edge A pointer to the edge this halfedge belongs to.
     * \param face  A pointer to  the face whose half-edge  cycle this
     * halfedge belongs to.
     * \param next  A pointer  to the next  halfedge in  the half-edge
     * cycle containing this halfedge.
     * \param prev A pointer to the previous halfedge in the half-edge
     * cycle containing this halfedge.
     */
    Halfedge(
	     Vertex* vertex ,
	     Edge* edge ,
	     Face* face ,
	     Halfedge* next ,
	     Halfedge* prev
	    )
    {
      set_origin( vertex ) ;
      set_edge( edge ) ;
      set_face( face ) ;
      set_next( next ) ;
      set_prev( prev ) ;
    }


    /**
     * \fn ~Halfedge()
     *
     * \brief Destroys an instance of this class.
     */
    ~Halfedge()
    {
      set_origin( 0 ) ;
      set_edge( 0 ) ;
      set_face( 0 ) ;
      set_next( 0 ) ;
      set_prev( 0 ) ;
    }


    /**
     * \fn Vertex* get_origin() const
     *
     * \brief  Returns  a  pointer   to  the  origin  vertex  of  this
     * half-edge.
     *
     * \return A pointer to the origin vertex of this half-edge.
     */
    Vertex* get_origin() const
    {
      return _vertex ; 
    }


    /**
     * \fn void set_origin( Vertex* vertex )
     *
     * \brief Assigns an  address to the pointer to  the origin vertex
     * of this half-edge.
     *
     * \param  vertex  The  address  of  the  origin  vertex  of  this
     * half-edge.
     */
    void set_origin( Vertex* vertex )
    {
      _vertex = vertex ; 
    }


    /**
     * \fn Edge* get_edge() const
     *
     * \brief Returns a pointer to the edge this half-edge belongs to.
     *
     * \return A pointer to the edge this half-edge belongs to.
     */
    Edge* get_edge() const
    {
      return _edge ; 
    }


    /**
     * \fn void set_edge( Edge* edge )
     *
     * \brief  Assigns an  address to  the  pointer to  the edge  this
     * half-edge belongs to.
     *
     * \param edge The address of the edge this half-edge belongs to.
     */
    void set_edge( Edge* edge )
    {
      _edge = edge ; 
    }


    /**
     * \fn Face* get_face() const
     *
     * \brief Returns a pointer to the face whose half-edge cycle this
     * half-edge belongs to.
     *
     * \return  A  pointer to  the  face  whose  half-edge cycle  this
     * half-edge belongs to.
     */
    Face* get_face() const
    {
      return _face ; 
    }


    /**
     * \fn void set_face( Face* face )
     *
     * \brief  Assigns an  address to  the pointer  to the  face whose
     * half-edge cycle this half-edge belongs to.
     *
     * \param face The address of  the face whose half-edge cycle this
     * half-edge belongs to.
     */
    void set_face( Face* face )
    {
      _face = face ; 
    }


    /**
     * \fn Halfedge* get_next() const
     *
     * \brief Returns a pointer to the next half-edge in the half-edge
     * cycle this half-edge belongs to.
     *
     * \return A pointer to the  next half-edge in the half-edge cycle
     * this half-edge belongs to.
     */
    Halfedge* get_next() const
    {
      return _next ; 
    }


    /**
     * \fn void set_next( Halfedge* next )
     *
     * \brief Assigns an address to  the pointer to the next half-edge
     * in the face half-edge cycle this half-edge belongs to.
     *
     * \param next The  address of  the  next half-edge  in the  face
     * half-edge cycle this half-edge belongs to.
     */
    void set_next( Halfedge* next )
    {
      _next = next ; 
    }


    /**
     * \fn Halfedge* get_prev() const
     *
     * \brief  Returns a  pointer  to the  previous  half-edge in  the
     * half-edge cycle this half-edge belongs to.
     *
     * \return A  pointer to the  previous half-edge in  the half-edge
     * cycle this half-edge belongs to.
     */
    Halfedge* get_prev() const
    {
      return _prev ; 
    }


    /**
     * \fn void set_prev( Halfedge* prev )
     *
     * \brief  Assigns  an address  to  the  pointer  to the  previous
     * half-edge in  the face  half-edge cycle this  half-edge belongs
     * to.
     *
     * \param prev The  address of the previous half-edge  in the face
     * half-edge cycle this half-edge belongs to.
     */
    void set_prev( Halfedge* prev )
    {
      _prev = prev ; 
    }


    /**
     * \fn Halfedge* get_mate() const
     *
     * \brief  Returns a  pointer to  the mate  of this  half-edge (if
     * any).
     *
     * \return A pointer to the mate of this half-edge.
     */
    Halfedge* get_mate() const
    {
      if ( get_edge()->get_first_halfedge() == this ) {
	return get_edge()->get_second_halfedge() ;
      }

      return get_edge()->get_first_halfedge() ;
    }

     /**
     * \fn HAttrib& get_attributes()
     *
     * \brief Returns the set of attributes associated with this halfedge.
     *
     * \return A reference to the set of attributes associated with this
     * halfedge.
     */
    HAttrib& get_attributes()
    {
      return _attributes ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    Vertex* _vertex ;  ///< Pointer to the origin vertex of this half-edge.

    Edge* _edge ;  ///< Pointer to the edge this half-edge belongs to.

    Face* _face ;  ///< Pointer to the face whose half-edge cycle this half-edge belongs to.

    Halfedge* _next ;  ///< Pointer to the next half-edge in the cycle this half-edge belongs to.

    Halfedge* _prev ;  ///< Pointer to the previous half-edge in the cycle this half-edge belongs to.

    HAttrib _attributes ;   ///< Set of attributes associated with this halfedge.

  } ;

}

/** @} */ //end of group class.

#endif   // HALFEDGE_H
