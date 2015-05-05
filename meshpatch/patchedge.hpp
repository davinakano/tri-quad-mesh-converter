/** 
 * \file patchedge.hpp
 *  
 * \brief Definition and implementation  of the class PatchEdge, which
 * represents  an   edge  of  a  triangulation  of   a  surface  patch
 * homeomorphic to a disk.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date October 2010
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef PATCHEDGE_HPP
#define PATCHEDGE_HPP


/** 
 * \defgroup MESHPATCHNameSpace Namespace meshpatch.
 * @{
 */

/**
 * \namespace meshpatch
 *
 * \brief   The  namespace  meshpatch   contains  the   definition  and
 * implementation of  all classes of the data  structure MeshPatch data
 * structure,  which can  be used  to  represent a  triangulation of  a
 * surface patch homeomorphic  to a disk.
 */


// --------------------------------------------------------------------
//
// The MeshPatch data  structure is based on the  "directed edge" data
// structure. See  paper "Directed Edges --  A Scalable Representation
// for  Triangular Meshes"  by  S.  Campagna,  L.  Kobbelt, and  H.-P.
// Seidel at the  Journal of graphics Tools, Vol. 3,  No. 4, pp. 1-12,
// 1998.
//
// --------------------------------------------------------------------

namespace meshpatch {

  /**
   * \class PatchEdge
   *
   * \brief This  class represents an  edge of the triangulation  of a
   * surface patch homeomorphic to a closed disk in \f$R^2\f$.
   */
  template < 
             typename VAttrib = int ,
             typename EAttrib = int ,
             typename PAttrib = int 
           >
  class PatchEdge {
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn PatchEdge( unsigned origin , unsigned mate , unsigned prev , bool status = false )
     *
     * \brief Creates an instance of this class.
     *
     * \param origin Index of the origin vertex of this edge.
     * \param mate Index of the mate edge of this edge or its index in
     * the boundary edge array if this edge is a boundary edge.
     * \param prev Index  of the edge preceding this  edge in the same
     * face cycle.
     * \param status A flag to indicate if this edge is a boundary edge.
     */
    PatchEdge(
	      unsigned origin ,
	      unsigned mate ,
              unsigned prev ,
	      bool status = false 
	     )
    {
      set_origin( origin ) ;
      set_mate( mate ) ;
      set_prev( prev ) ;
      set_boundary_status( status ) ;
    
      return ;
    }


    /**
     * \fn PatchEdge( const PatchEdge< VAttrib , EAttrib , PAttrib >& patch )
     *
     * \brief Creates an instance of this class.
     *
     * \param patch A reference to an object to be cloned.
     */
    PatchEdge( const PatchEdge< VAttrib , EAttrib , PAttrib >& patch )
    {
      set_origin( patch.get_origin() ) ;
      set_mate( patch.get_mate() ) ;
      set_prev( patch.get_prev() ) ;
      set_boundary_status( patch.get_boundary_status() ) ;
      this->_attributes = patch._attributes ;
    
      return ;
    }


    /**
     * \fn ~PatchEdge()
     *
     * \brief Releases the memory held by an instance of this class.
     */
    virtual ~PatchEdge() 
    {}


    /**
     * \fn inline unsigned get_origin() const
     *
     * \brief Returns the index of the origin vertex of this edge.
     *
     * \return The index of the origin vertex of this edge.
     */
    inline unsigned get_origin() const
    {
      return _origin ;
    }


    /**
     * \fn inline unsigned get_mate() const
     *
     * \brief Returns the  index of the mate of this edge.
     *
     * \return The index of the mate of this edge.
     */
    inline unsigned get_mate() const
    {
      return _mate ;
    }


    /**
     * \fn inline unsigned get_prev() const
     *
     * \brief Returns the index of the edge preceding this edge in the
     * same face cycle.
     *
     * \return The index  of the edge preceding this  edge in the same
     * face cycle.
     */
    inline unsigned get_prev() const
    {
      return _prev ;  
    }


    /**
     * \fn inline bool get_boundary_status() const
     *
     * \brief Returns the boundary status of this edge.
     *
     * \return The logic  value true if this edge  is a boundary edge,
     * and the logic value false otherwise.
     */
    inline bool get_boundary_status() const
    {
      return _status ;  
    }


    /**
     * \fn inline EAttrib& get_attributes()
     *
     * \brief Returns the set of attributes associated with this edge.
     *
     * \return The set of attributes associated with this edge.
     */
    inline EAttrib& get_attributes()
    {
      return _attributes ;  
    }


    /**
     * \fn inline bool is_boundary_edge() const
     *
     * \brief Determines whether this edge is a boundary edge.
     *
     * \return The logic  value true if this edge  is a boundary edge,
     * and the logic value false otherwise.
     */
    inline bool is_boundary_edge() const
    {
      return _status ;  
    }


    /**
     * \fn inline void set_origin( unsigned i )
     *
     * \brief Sets the origin vertex index of this edge.
     *
     * \param i The index of a surface patch vertex.
     */
    inline void set_origin( unsigned i )
    {
      _origin = i ;

      return ;
    }


    /**
     * \fn inline void set_mate( unsigned i )
     *
     * \brief Sets the mate index of this edge.
     *
     * \param i The index of a surface patch edge.
     */
    inline void set_mate( unsigned i )
    {
      _mate = i ;

      return ;
    }


    /**
     * \fn inline void set_prev( unsigned i )
     *
     * \brief Sets the previous edge index of this edge.
     *
     * \param i The index of a surface patch edge.
     */
    inline void set_prev( unsigned i )
    {
      _prev = i ;

      return ;
    }


    /**
     * \fn inline void set_boundary_status( bool status )
     *
     * \brief Sets the status of this edge.
     *
     * \param status A logic value to indicate the status (true if the
     * edge is a boundary edge, and false otherwise.)
     */
    inline void set_boundary_status( bool status )
    {
      _status = status ;

      return ;
    }


    /**
     * \fn inline void set_attributes( const EAttrib& attrib )
     *
     * \brief Assigns a set of attributes with this edge.
     *
     * \param attrib A set of edge attributes.
     */
    inline void set_attributes( const EAttrib& attrib )
    {
      _attributes = attrib ;

      return ;
    }


  protected:
    // ---------------------------------------------------------------
    //
    // Protected member functions
    //
    // ---------------------------------------------------------------

    unsigned _origin ;  ///< Index of the origin vertex of this edge.

    unsigned _mate ;  ///< Index of the mate of this edge.

    unsigned _prev ;  ///< Index of the edge preceding this one in the same face cycle.

    bool _status ;   ///< Status of this edge.

    EAttrib _attributes ;  ///< Set of attributes associated with this edge. 

  } ;

}

/** @} */ //end of group class.

#endif   // PATCHEDGE_HPP
