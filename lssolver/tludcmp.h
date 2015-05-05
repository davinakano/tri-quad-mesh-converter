/** 
 * \file ludcmp.h
 *  
 * \brief Definition  of the class tLUdcmp, which  represents a linear
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

#ifndef tLUDCMP_H
#define tLUDCMP_H

#include <vector>    // std::vector


/** 
* \defgroup LSSOLVERNameSpace Namespace lssolver.
* @{
*/

/**
* \namespace lssolver
*
* \brief   The  namespace   lssolver  contains   the   definition  and
* implementation of classes representing numerical solvers for systems
* of linear equations.
*/

namespace lssolver {


  /**
   * \class tLUdcmp
   *
   * \brief This class  represents a linear system solver  based on LU
   * decomposition.
   *
   */
  class tLUdcmp {
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn tLUdcmp( double** mat , unsigned n )
     *
     * \brief Creates an instance of this class.
     *
     * \param mat A ( n x n ) matrix.
     * \param n The number of rows and columns of matrix mat.
     *
     */
    tLUdcmp( double** mat , unsigned n ) ;


    /**
     * \fn tLUdcmp( const tLUdcmp& lu )
     *
     * \brief Creates an instance of this class from another instance.
     *
     * \param lu An instance of class LUdcmp.
     *
     */
    tLUdcmp(const tLUdcmp& lu ) ;


    /**
     * \fn ~tLUdcmp() 
     *
     * \brief Destroys an instance of this class.
     *
     */
    ~tLUdcmp() ;


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

    double** m_matA;   ///< The coefficient matrix of the linear system.

    unsigned m_nele;   ///< The number of rows and columns of the matrix.

    std::vector<unsigned> m_perm;   ///< The permutation vector used by the LU decomposition algorithm

    int m_sign;        ///< The sign flag used by the LU decomposition algorithm

  } ;

}

/** @} */ //end of group class.

#endif   // LUDCMP_H
