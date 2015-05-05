/** 
 * \file mypair.hpp
 *  
 * \brief  Definition and  implementation of  the class  MyPair, which
 * represents an ordered pair of generic type objects.
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

#ifndef MYPAIR_HPP
#define MYPAIR_HPP


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
   * \class MyPair
   *
   * \brief  This class  represents  an ordered  pair of  parametrized
   * types.
   */
  template < 
             typename X = unsigned ,
             typename Y = unsigned  
           >
  class MyPair {
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn MyPair( const X& first , const Y& second )
     *
     * \brief Creates an instance of this class.
     *
     * \param first The first element of the pair.
     * \param second The second element of the pair.
     */
    MyPair( const X& first , const Y& second ) 
    {
      _first = first ;
      _second = second ;
    }


    /**
     * \fn MyPair( const MyPair< X , Y >& pair )
     *
     * \brief Creates an instance  of this class from another instance
     * of this class.
     *
     * \param pair An instance of this class.
     */
    MyPair( const MyPair< X , Y >& pair  )
    {
      _first = pair._first ;
      _second = pair._second ;
    }


    /** 
     * \fn inline X first() const
     *
     * \brief Returns the first element of the pair.
     *
     * \return The first element of the pair.
     */
    inline X first() const
    {
      return _first ;
    }


    /** 
     * \fn inline Y second() const
     *
     * \brief Returns the second element of the pair.
     *
     * \return The second element of the pair.
     */
    inline Y second() const
    {
      return _second ;
    }


    /** 
     * \fn inline MyPair< X , Y >& operator=( const MyPair< X , Y >& pair )
     *
     * \brief Overloads assignment operator.
     *
     * \param pair An instance of this class.
     *
     * \return A reference to this instance.
     */
    MyPair< X , Y >& operator=( const MyPair< X , Y >& pair )
    {
      _first = pair._first ;
      _second = pair._second ;

      return *this ;
    }


    /** 
     * \fn inline bool operator==( const MyPair< X , Y >& pair ) const 
     *
     * \brief  Compares  a given  instance  of  this  class with  this
     * instance, and  returns true if  and only if both  instances are
     * the same.
     *
     * \param pair An instance of this class.
     *
     * \return The  logic value  true if this  instance and  the given
     * instance are the same, and the logic value false otherwise.
     */
    bool operator==( const MyPair< X , Y >& pair ) const 
    {
      return ( _first == pair._first ) && ( _second == pair._second ) ;
    }


    /** 
     * \fn inline bool operator<( const MyPair< X , Y >& pair ) const 
     *
     * \brief  Compares  a given  instance  of  this  class with  this
     * instance, and returns true if and only this instance is smaller
     * than the latter.
     *
     * \param pair An instance of this class.
     *
     * \return The logic  value true if this instance  is smaller than
     * the given instance, and the logic value false otherwise.
     */
    bool operator<( const MyPair< X , Y >& pair ) const 
    {
      if ( _first < pair._first ) {
	return true ;
      }
      else if ( _first == pair._first ) {
	return _second < pair._second ;
      }
      
      return false ;
    }


  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    X _first ;  ///< First element of this ordered pair.
    Y _second ;  ///< Second element of this ordered pair.

  } ;


}

/** @} */ //end of group class.

#endif   // MYPAIR_HPP
