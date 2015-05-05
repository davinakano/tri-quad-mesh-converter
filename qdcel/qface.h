/** 
 * \file qface.h
 *  
 * \brief Definition of the class QFace, which represents a face (i.e.,
 * a  quadrilateral)  from  a  surface  mesh represented  by  a  QDCEL  data
 * structure.
 *
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

#ifndef QFACE_H
#define QFACE_H

#include "qhalfedge.h"   // Halfedge


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
   * \class QFace
   *
   * \brief This  class represents  a face (i.e.,  a quadrilateral)  from a
   * surface mesh represented by the QDCEL data structure.
   */
template <
        typename VAttrib ,
        typename FAttrib ,
        typename EAttrib ,
        typename HAttrib
        >
class QFace {
public:
    // ---------------------------------------------------------------
    //
    // Type definitions
    //
    // ---------------------------------------------------------------

    /**
     * \typedef QHalfedge
     *
     * \brief  Defines  QHalfedge   as  an  alias  for  qdcel::QHalfedge<
     * VAttrib, FAttrib , EAttrib , HAttrib >.
     */    
    typedef qdcel::QHalfedge< VAttrib, FAttrib , EAttrib , HAttrib >
    QHalfedge ;


    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /** 
     * \fn QFace( QHalfedge* h )
     *
     * \brief Creates an instance of this class.
     *
     * \param h A pointer to the first halfedge of this face.
     */
    QFace( QHalfedge* h ) : _halfedge( h )
    {}


    /**
     * \fn ~QFace()
     *
     * \brief Destroys an instance of this class.
     */
    ~QFace()
    {
        set_halfedge( 0 ) ;
    }


    /**
     * \fn QHalfedge* get_halfedge() const
     *
     * \brief Returns a pointer to the first half-edge of this face.
     *
     * \return A pointer to the first half-edge of this face.
     */
    QHalfedge* get_halfedge() const
    {
        return _halfedge ;
    }


    /**
     * \fn void set_halfedge( QHalfedge* h )
     *
     * \brief Assigns an address to the pointer to the first half-edge
     * of this face.
     *
     * \param h The address of the first half-edge of this face.
     */
    void set_halfedge( QHalfedge* h )
    {
        _halfedge = h ;
    }


    /**
     * \fn FAttrib& get_attributes()
     *
     * \brief Returns the set of attributes associated with this face.
     *
     * \return A reference to the set of attributes associated with this
     * face.
     */
    FAttrib& get_attributes()
    {
        return _attributes ;
    }



private:
    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    QHalfedge* _halfedge ;   ///< Pointer to the first half-edge of the half-edge cycle of this face.


    FAttrib _attributes ;   ///< Set of attributes associated with this face. 

} ;

}

/** @} */ //end of group class.

#endif   // QFACE_H

