/** 
 * \file ludcmp.h
 *  
 * \brief  Implementation of  the  class tLUdcmp,  which represents  a
 * linear  system solver  based on  LU decomposition.   This  class is
 * basically a wrapper for the  LU decomposition code described in the
 * following book:
 *
 * William H. Press; Saul.  A. Teukolsky; William T. Vetterling; Brian
 * P.   Flannery.   Numerical Recipes  in  C:  The  Art of  Scientific
 * Computing. Second Edition, Cambridge University Press, 1992.
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

#include "tludcmp.h"

#include <cassert>    // assert
#include <cmath>      // fabs


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
   * \fn tLUdcmp( double** mat , unsigned n )
   *
   * \brief Creates an instance of this object.
   *
   * \param mat A ( n x n ) matrix.
   * \param n The number of rows and columns of matrix mat.
   *
   */
  tLUdcmp::tLUdcmp(
		   double** mat,
		   unsigned n
		  )
  {
    m_nele = n;
    m_matA = allocate( m_nele ) ;

    for ( unsigned i = 0 ; i < m_nele ; i++ ) {
      for ( unsigned j = 0 ; j < m_nele ; j++ ) {
        m_matA[ i ][ j ] = mat[ i ][ j ] ;       
      }
    }

    m_perm.resize( m_nele ) ;
    m_sign = 1 ;

    decomp() ;

    return ;
  }


  /**
   * \fn tLUdcmp( const tLUdcmp& lu )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param lu An instance of class LUdcmp.
   *
   */
  tLUdcmp::tLUdcmp( const tLUdcmp& lu ) : m_nele( lu.m_nele )
  {
    m_matA = allocate( m_nele ) ;

    m_perm.resize( m_nele ) ;
    m_sign = lu.m_sign ;

    for ( unsigned i = 0 ; i < m_nele; i++ ) {
      m_perm[ i ] = lu.m_perm[ i ] ;

      for ( unsigned j = 0 ; j < m_nele ; j++ ) {
        m_matA[ i ][ j ] = lu.m_matA[ i ][ j ] ;        
      }
    }

    return ;
  }


  /**
   * \fn ~tLUdcmp() 
   *
   * \brief Destroys an instance of this class.
   *
   */
  tLUdcmp::~tLUdcmp()
  {
    deallocate( m_matA , m_nele ) ;

    return ;
  }


  /**
   * \fn void solve( double* b , unsigned n ) const 
   *
   * \brief Solves a linear system using LU decomposition.
   *
   * \param b The column (and solution) vector.
   * \param n The number of elements of the vector 
   */
  void
  tLUdcmp::solve(
		 double* b,
		 unsigned n
		)
    const
  {
    assert( n == m_nele ) ;

    int i, ii = 0, ip, j;
    double sum;
    for ( i = 0 ; i < (int) n ; i++ ) {
      ip = m_perm[ i ] ;
      sum = b[ ip ] ;
      b[ ip ] = b[ i ] ;
      if ( ii != 0 ) {
        for ( j = ii - 1 ; j < i ; j++ ) {
          sum -= m_matA[ i ][ j ] * b[ j ] ;
        }
      }
      else if ( sum != 0.0 ) {
        ii = i + 1 ;
      }

      b[ i ] = sum ;
    }

    for ( i = n - 1 ; i >= 0 ; i-- ) {
      sum = b[ i ] ;
      for ( j = i + 1 ; j < (int) n ; j++ ) {
        sum -= m_matA[ i ][ j ] * b[ j ] ;
      }

      b[ i ] = sum / m_matA[ i ][ i ] ;
    }

    return ;
  }


  /**
   * \fn double** allocate( unsigned n ) const
   *
   * \brief Allocates memory for a ( n x n ) matrix.
   *
   * \param n The number of rows and columns of the matrix.
   *
   * \return The address of the matrix.
   */
  double**
  tLUdcmp::allocate( unsigned n ) const
  {
    double** mat = (double **) new double*[n];

    for ( unsigned i = 0 ; i < n ; i++ ) {
      mat[ i ] = ( double* ) new double[ n ] ;
    }

    return mat ;
  }


  /**
   * \fn void deallocate( double** mat , unsigned n ) const
   *
   * \brief Releases the memory held by a ( n x n ) matrix.
   *
   * \param mat The address of the matrix.
   * \param n The number of rows and columns of the matrix.
   */
  void
  tLUdcmp::deallocate(
		      double** mat,
		      unsigned n
		     )
    const
  {
    for ( unsigned i = 0 ; i < n ; i++ ) {
      delete[] mat[ i ] ;
    }

    delete mat ;

    return ;
  }


  /**
   * \fn void decomp()
   *
   * \brief Carries out a LU  decomposition of the matrix of stored in
   * \var{m_matA}.
   *
   */
  void
  tLUdcmp::decomp()
  {
    const double TINY = 1.0e-20 ;
    int i , imax=0 , j , k ; 
    double big , dum , sum , temp ;

    int n = ( int ) m_nele ;
    std::vector< double > vv( n ) ;
    m_sign = 1 ;
    for ( i = 0 ; i < n ; i++ ) {
      big = 0.0 ;
      for ( j = 0 ; j < n ; j++ ) {
        if ( ( temp = fabs( m_matA[ i ][ j ] ) ) > big ) {
          big = temp ;
        }
      }

      assert( big != 0.0 ) ;

      vv[ i ] = 1.0 / big ;
    }

    for ( j = 0 ; j < n ; j++ ) {
      for ( i = 0 ; i < j ; i++ ) {
        sum = m_matA[ i ][ j ] ;

        for ( k = 0 ; k < i ; k++ ) {
          sum -= m_matA[ i ][ k ] * m_matA[ k ][ j ] ;
        }

        m_matA[ i ][ j ] = sum ;
      }

      big = 0.0 ;

      for ( i = j ; i < n ; i++ ) {
        sum = m_matA[ i ][ j ] ;
        for ( k = 0 ; k < j ; k++ ) {
          sum -= m_matA[ i ][ k ] * m_matA[ k ][ j ] ;
        }

        m_matA[ i ][ j ] = sum ;

        if ( ( dum = vv[ i ] * fabs( sum ) ) >= big ) {
          big = dum ;
          imax = i ;
	}
      }

      if ( j != imax ) {
        for ( k = 0 ; k < n ; k++ ) {
          dum = m_matA[ imax ][ k ] ;
          m_matA[ imax ][ k ] = m_matA[ j ][ k ] ;
          m_matA[ j ][ k ] = dum ;
	}

        m_sign = -m_sign ;
        vv[ imax ] = vv[ j ] ;
      }

      m_perm[ j ] = imax ;

      if ( m_matA[ j ][ j ] == 0.0 ) {
        m_matA[ j ][ j ] = TINY ;
      }

      if ( j != ( n - 1 ) ) {
        dum = 1.0 / ( m_matA[ j ][ j ] ) ;
        for ( i = j + 1 ; i < n ; i++ ) {
          m_matA[ i ][ j ] *= dum ;
        }
      }
    }

    return ;
  }

}

/** @} */ //end of group class.
