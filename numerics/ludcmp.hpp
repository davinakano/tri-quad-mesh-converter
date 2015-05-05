/** 
 * \file ludcmp.hpp
 *  
 * \brief Definition  of the class  LUdcmp, which represents  a linear
 * system solver based on LU  decomposition. This class is basically a
 * wrapper for  the LU decomposition  code described in  the following
 * book:
 *
 * William H. Press; Saul.  A. Teukolsky; William T. Vetterling; Brian
 * P.   Flannery.   Numerical Recipes  in  C:  The  Art of  Scientific
 * Computing. Second Edition, Cambridge University Press, 1992.
 *
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
 * \date July 2009
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef LUDCMP_HPP
#define LUDCMP_HPP

#include <vector>    // std::vector


/** 
 * \defgroup NumericsNameSpace Namespace numerics.
 * @{
 */

/**
 * \namespace numerics
 *
 * \brief   The  namespace  numerics   contains  the   definition  and
 * implementation  of  several  classes  for  carrying  out  numerical
 * computations.
 *
 */

namespace numerics {


  /**
   * \class LUdcmp
   *
   * \brief This class  represents a linear system solver  based on LU
   * decomposition.
   *
   */
  class LUdcmp {
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn LUdcmp( double** mat , unsigned n )
     *
     * \brief Creates an instance of this class.
     *
     * \param mat A ( n x n ) matrix.
     * \param n The number of rows and columns of matrix mat.
     *
     */
    LUdcmp( double** mat , unsigned n ) ;


    /**
     * \fn LUdcmp( const LUdcmp& lu )
     *
     * \brief Creates an instance of this class from another instance.
     *
     * \param lu An instance of class LUdcmp.
     *
     */
    LUdcmp(const LUdcmp& lu ) ;


    /**
     * \fn ~LUdcmp() 
     *
     * \brief Destroys an instance of this class.
     *
     */
    ~LUdcmp() ;


    /**
     * \fn void solve( double* b , unsigned n ) const 
     *
     * \brief Solves a linear system using LU decomposition.
     *
     * \param b The column (and solution) vector.
     * \param n The number of elements of the vector 
     */
    void solve( double* b , unsigned n ) const;
    

  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn double** allocate( unsigned n ) const
     *
     * \brief Allocates memory for a ( n x n ) matrix.
     *
     * \param n The number of rows and columns of the matrix.
     *
     * \return The address of the matrix.
     */
    double** allocate( unsigned n ) const ;


    /**
     * \fn void deallocate( double** mat , unsigned n ) const
     *
     * \brief Releases the memory held by a ( n x n ) matrix.
     *
     * \param mat The address of the matrix.
     * \param n The number of rows and columns of the matrix.
     */
    void deallocate( double** mat , unsigned n ) const ;


    /**
     * \fn void decomp()
     *
     * \brief Carries out  a LU decomposition of the  matrix of stored
     * in \var{_matA}.
     *
     */
    void decomp() ;


    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------

    double** _matA;   ///< The coefficient matrix of the linear system.

    unsigned _nele;   ///< The number of rows and columns of the matrix.

    std::vector<unsigned> _perm;   ///< The permutation vector used by the LU decomposition algorithm

    int _sign;        ///< The sign flag used by the LU decomposition algorithm

  } ;

}

/** @} */ //end of group class.

#endif   // LUDCMP_HPP
