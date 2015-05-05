/** 
 * \file qhalfedge.h
 *  
 * \brief  Definition  of  the  class  QHalfedge,  which  represents  a
 * half-edge from a quad surface mesh represented by the QDCEL data
 * structure.
 *
 * \author
 * Mario Augusto de Souza Lizier \n
 * Universidade Federal de Sao Carlos \n
 * Departamento de Computacao \n
 * lizier at dc (dot) ufscar (dot) br \n
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date July 2011
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef QHALFEDGE_H
#define QHALFEDGE_H

#include "qvertex.h"   // QVertex
#include "qedge.h"     // QEdge


/**
* \defgroup QDCELNameSpace Namespace qdcel.
* @{
*/

/**
* \namespace qdcel
*
* \brief The namespace dcel contains the definition and implementation
* of  all classes  of the  data structure  Doubly Connected  Edge List for quads
* (QDCEL), which  can be  used to represent  surface meshes  with empty
* boundary.
*
*/

namespace qdcel {

/**
   * Forward definition of classes.
   */
template <
        typename VAttrib ,
        typename FAttrib ,
        typename EAttrib ,
        typename HAttrib
        >
class QFace ;


/**
   * \class QHalfedge
   *
   * \brief This class represents  a half-edge (i.e., an oriented line
   * segment)  from a quad  surface mesh  represented by  the QDCEL
   * data structure.
   */
template <
        typename VAttrib ,
        typename FAttrib ,
        typename EAttrib ,
        typename HAttrib
        >
class QHalfedge {
public:
    // ---------------------------------------------------------------
    //
    // Type definitions
    //
    // ---------------------------------------------------------------

    /**
     * \typedef QVertex
     *
     * \brief Defines  Vertex as  an alias for  dcel::Vertex< VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */
    typedef qdcel::QVertex< VAttrib, FAttrib , EAttrib , HAttrib > QVertex ;

    /**
     * \typedef QEdge
     *
     * \brief  Defines  Edge  as  an alias  for  dcel::Edge<  VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef qdcel::QEdge< VAttrib, FAttrib , EAttrib , HAttrib > QEdge ;

    /**
     * \typedef QFace
     *
     * \brief  Defines  Face  as  an alias  for  dcel::Face<  VAttrib,
     * FAttrib , EAttrib , HAttrib >.
     */    
    typedef qdcel::QFace< VAttrib, FAttrib , EAttrib , HAttrib > QFace ;


    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /** 
     * \fn QHalfedge( QVertex* vertex , QEdge* edge , QFace* face , QHalfedge* next , QHalfedge* prev )
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
    QHalfedge(
        QVertex* vertex ,
        QEdge* edge ,
        QFace* face ,
        QHalfedge* next ,
        QHalfedge* prev
        )
    {
        set_origin( vertex ) ;
        set_edge( edge ) ;
        set_face( face ) ;
        set_next( next ) ;
        set_prev( prev ) ;
    }


    /**
     * \fn ~QHalfedge()
     *
     * \brief Destroys an instance of this class.
     */
    ~QHalfedge()
    {
        set_origin( 0 ) ;
        set_edge( 0 ) ;
        set_face( 0 ) ;
        set_next( 0 ) ;
        set_prev( 0 ) ;
    }


    /**
     * \fn QVertex* get_origin() const
     *
     * \brief  Returns  a  pointer   to  the  origin  vertex  of  this
     * half-edge.
     *
     * \return A pointer to the origin vertex of this half-edge.
     */
    QVertex* get_origin() const
    {
        return _vertex ;
    }


    /**
     * \fn void set_origin( QVertex* vertex )
     *
     * \brief Assigns an  address to the pointer to  the origin vertex
     * of this half-edge.
     *
     * \param  vertex  The  address  of  the  origin  vertex  of  this
     * half-edge.
     */
    void set_origin( QVertex* vertex )
    {
        _vertex = vertex ;
    }


    /**
     * \fn QEdge* get_edge() const
     *
     * \brief Returns a pointer to the edge this half-edge belongs to.
     *
     * \return A pointer to the edge this half-edge belongs to.
     */
    QEdge* get_edge() const
    {
        return _edge ;
    }


    /**
     * \fn void set_edge( QEdge* edge )
     *
     * \brief  Assigns an  address to  the  pointer to  the edge  this
     * half-edge belongs to.
     *
     * \param edge The address of the edge this half-edge belongs to.
     */
    void set_edge( QEdge* edge )
    {
        _edge = edge ;
    }


    /**
     * \fn QFace* get_face() const
     *
     * \brief Returns a pointer to the face whose half-edge cycle this
     * half-edge belongs to.
     *
     * \return  A  pointer to  the  face  whose  half-edge cycle  this
     * half-edge belongs to.
     */
    QFace* get_face() const
    {
        return _face ;
    }


    /**
     * \fn void set_face( QFace* face )
     *
     * \brief  Assigns an  address to  the pointer  to the  face whose
     * half-edge cycle this half-edge belongs to.
     *
     * \param face The address of  the face whose half-edge cycle this
     * half-edge belongs to.
     */
    void set_face( QFace* face )
    {
        _face = face ;
    }


    /**
     * \fn QHalfedge* get_next() const
     *
     * \brief Returns a pointer to the next half-edge in the half-edge
     * cycle this half-edge belongs to.
     *
     * \return A pointer to the  next half-edge in the half-edge cycle
     * this half-edge belongs to.
     */
    QHalfedge* get_next() const
    {
        return _next ;
    }


    /**
     * \fn void set_next( QHalfedge* next )
     *
     * \brief Assigns an address to  the pointer to the next half-edge
     * in the face half-edge cycle this half-edge belongs to.
     *
     * \param next The  address of  the  next half-edge  in the  face
     * half-edge cycle this half-edge belongs to.
     */
    void set_next( QHalfedge* next )
    {
        _next = next ;
    }


    /**
     * \fn QHalfedge* get_prev() const
     *
     * \brief  Returns a  pointer  to the  previous  half-edge in  the
     * half-edge cycle this half-edge belongs to.
     *
     * \return A  pointer to the  previous half-edge in  the half-edge
     * cycle this half-edge belongs to.
     */
    QHalfedge* get_prev() const
    {
        return _prev ;
    }


    /**
     * \fn void set_prev( QHalfedge* prev )
     *
     * \brief  Assigns  an address  to  the  pointer  to the  previous
     * half-edge in  the face  half-edge cycle this  half-edge belongs
     * to.
     *
     * \param prev The  address of the previous half-edge  in the face
     * half-edge cycle this half-edge belongs to.
     */
    void set_prev( QHalfedge* prev )
    {
        _prev = prev ;
    }


    /**
     * \fn QHalfedge* get_mate() const
     *
     * \brief  Returns a  pointer to  the mate  of this  half-edge (if
     * any).
     *
     * \return A pointer to the mate of this half-edge.
     */
    QHalfedge* get_mate() const
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

    QVertex* _vertex ;  ///< Pointer to the origin vertex of this half-edge.

    QEdge* _edge ;  ///< Pointer to the edge this half-edge belongs to.

    QFace* _face ;  ///< Pointer to the face whose half-edge cycle this half-edge belongs to.

    QHalfedge* _next ;  ///< Pointer to the next half-edge in the cycle this half-edge belongs to.

    QHalfedge* _prev ;  ///< Pointer to the previous half-edge in the cycle this half-edge belongs to.

    HAttrib _attributes ;   ///< Set of attributes associated with this halfedge.

} ;

}

/** @} */ //end of group class.

#endif   // QHALFEDGE_H
