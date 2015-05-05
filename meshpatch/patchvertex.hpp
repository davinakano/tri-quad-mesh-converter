/** 
 * \file patchvertex.hpp
 *  
 * \brief  Definition  and implementation  of  the class  PatchVertex,
 * which represents  a vertex  of a triangulation  of a  surface patch
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

#ifndef PATCHVERTEX_HPP
#define PATCHVERTEX_HPP

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
   * \class PatchVertex
   *
   * \brief This class  represents a vertex of the  triangulation of a
   * surface patch homeomorphic to a closed disk in \f$R^2\f$.
   */
  template < 
             typename VAttrib = int ,
             typename EAttrib = int ,
             typename PAttrib = int 
           >
  class PatchVertex {
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn PatchVertex()
     *
     * \brief Creates an instance of this class.
     */
    PatchVertex()
    {}


    /**
     * \fn PatchVertex( unsigned e , double x , double y , double z )
     *
     * \brief Creates an instance of this class.
     *
     * \param  e Index  of a  surface  patch edge  incident with  this
     * vertex.
     * \param x The first Cartesian coordinate of the vertex location.
     * \param  y  The  second   Cartesian  coordinate  of  the  vertex
     * location.
     * \param z The third Cartesian coordinate of the vertex location.
     */
    PatchVertex(
		unsigned e ,
		double x ,
		double y ,
		double z
	       )
    {
      set_edge( e ) ;
      set_x_coord( x ) ;
      set_y_coord( y ) ;
      set_z_coord( z ) ;

      return ;
    }


    /**
     * \fn PatchVertex( const PatchVertex< VAttrib , EAttrib , PAttrib >& patch )
     *
     * \brief Creates an instance of this class.
     *
     * \param patch A reference to an object to be cloned.
     */
    PatchVertex(  const PatchVertex< VAttrib , EAttrib , PAttrib >& patch ) 
    {
      set_edge( patch.get_edge() ) ;
      set_x_coord( patch.x() ) ;
      set_y_coord( patch.y() ) ;
      set_z_coord( patch.z() ) ;
      this->_attributes = patch._attributes ;

      return ;
    }


    /**
     * \fn ~PatchVertex()
     *
     * \brief Releases the memory held by an instance of this class.
     */
    virtual ~PatchVertex() 
    {}


    /**
     * \fn inline unsigned get_edge() const
     *
     * \brief Returns the index of an edge with origin at this vertex.
     *
     * \return The index of an edge with origin at this vertex.
     */
    inline unsigned get_edge() const
    {
      return _edge ;
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
     * \fn inline VAttrib& get_attributes()
     *
     * \brief  Returns  the set  of  attributes  associated with  this
     * vertex.
     *
     * \return The set of attributes associated with this vertex.
     */
    inline VAttrib& get_attributes()
    {
      return _attributes ;  
    }


    /**
     * \fn inline void set_edge( unsigned e )
     *
     * \brief Sets origin edge index of this vertex .
     *
     * \param e The index of a surface patch edge.
     */
    inline void set_edge( unsigned e )
    {
      _edge = e ;

      return ;
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
     * \fn inline void set_attributes( const VAttrib& attrib )
     *
     * \brief Assigns a set of attributes with this vertex.
     *
     * \param attrib A set of vertex attributes.
     */
    inline void set_attributes( const VAttrib& attrib )
    {
      _attributes = attrib ;

      return ;
    }


    // ---------------------------------------------------------------
    //
    // Protected data members
    //
    // ---------------------------------------------------------------

    unsigned _edge ;  ///< Index of an edge with origin at this vertex.

    double _x ;  ///< First Cartesian coordinate of this vertex.

    double _y ;  ///< Second Cartesian coordinate of this vertex.

    double _z ;  ///< Third Cartesian coordinate of this vertex.

    VAttrib _attributes ;  ///< Set of attributes associated with this vertex. 

  } ;

}

/** @} */ //end of group class.

#endif   // PATCHVERTEX_HPP

