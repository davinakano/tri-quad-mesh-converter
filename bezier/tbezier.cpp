/** 
 * \file bezier.cpp
 *  
 * \brief  Implementation  of  the  class tBezier  that  represents  a
 * rectangular B&eacute;zier patch in 3D.
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

#include "tbezier.h"      // tBezier

#include "lssolver/tludcmp.h"      // lssolver::tLUdcmp

#include <cassert>       // assert
#include <cmath>         // fabs

/** 
* \defgroup BEZIERNameSpace Namespace bezier.
* @{
*/

/**
* \namespace bezier
*
* \brief   The   namespace   bezier   contains  the   definition   and
* implementation of  a class representing  a rectangular B&eacute;zier
* patch in 3D.
*
*/

namespace bezier {


  /**
   * \fn tBezier( double* param_pts , double* image_pts , unsigned np , unsigned m , unsigned n , double rx , double sx , double ry , double sy )
   *
   * \brief  Creates an  instance of  this Bezier  class by  solving a
   * least squares fitting problem.  In particular, the control points
   * of the patch  are computed as the solution of  a linear system of
   * normal equations.  To set up  this system, the method  receives a
   * set of  parameter points and their corresponding  image points in
   * 3D. The patch is built so that it approximates the latter points.
   *
   * \param param_pts  an array  with the 2D  coordinates of a  set of
   * parameter points.
   * \param image_pts An array with the  3D coordinates of a set of 3D
   * points,  which are  the image  of the  paramenter points  under a
   * possibly unknown function.
   * \param np The number of parameter points.
   * \param n The first degree of this rectangular Bezier patch.
   * \param m The second degree of this rectangular Bezier patch.
   * \param rx  The first Cartesian  coordinate of the first  point of
   * the affine frame w.r.t. which this surface patch is defined.
   * \param ry The second coordinate  of the first point of the affine
   * frame w.r.t. which this surface patch is defined.
   * \param sx The  first Cartesian coordinate of the  second point of
   * the affine frame w.r.t. which this surface patch is defined.
   * \param sy The second Cartesian  coordinate of the second point of
   * the affine frame w.r.t. which this surface patch is defined.
   *
   */
  tBezier::tBezier(
		   double* param_pts ,
		   double* image_pts ,
		   unsigned np ,
		   unsigned m ,
		   unsigned n ,
                   double rx ,
                   double sx ,
                   double ry ,
                   double sy
		  ) 
  {
    /**
     * Check consistency of the input parameters.
     */
    assert( m > 0 ) ;
    assert( n > 0 ) ;
    assert( rx < sx ) ;
    assert( ry < sy ) ;
    assert( ( 2 * ( n + 1 ) * ( m + 1 ) ) <= np ) ;

    /**
     * Initialize private attributes.
     */

    set_bidegree_1st_index( m ) ;
    set_bidegree_2nd_index( n ) ;

    set_aff_1st_point_x_coord( rx ) ;
    set_aff_1st_point_y_coord( ry ) ;

    set_aff_2nd_point_x_coord( sx ) ;
    set_aff_2nd_point_y_coord( sy ) ;

    /**
     * Allocate memory for the control points.
     */
    unsigned ncp = ( get_bidegree_1st_index() + 1 ) * 
                   ( get_bidegree_2nd_index() + 1 ) ;

    m_ctrl_pts = new double*[ 
			     ( get_bidegree_1st_index() + 1 ) * 
			     ( get_bidegree_2nd_index() + 1 ) 
			    ] ;

    for ( unsigned i = 0 ; i < ncp ; i++ ) {
      m_ctrl_pts[ i ] = new double[ 3 ] ;
    }

    /**
     * Compute a (np x ncp)  matrix "A" with all Bernstein polynomials
     * of order "m_m + m_n" at the parameter points given in the input
     * vector "param_pts".
     */

    double** A = new double*[ np ] ;

    for ( unsigned i = 0 ; i < np ; i++ ) {
      A[ i ] = new double[ ncp ] ;
    }

    comp_bpoly_matrix( A , param_pts , np ) ;

    /**
     * Compute the  product "At x A",  where "At" is  the transpose of
     * "A".
     */
    double** AtA = new double*[ ncp ] ;

    for ( unsigned i = 0 ; i < ncp ; i++ ) {
      AtA[ i ] = new double[ ncp ] ;
    }

    comp_matrix_ata( A , np , ncp , AtA ) ;

    /**
     * Compute the  product "At x B",  where "At" is  the transpose of
     * "A".
     */
    double** AtB = new double*[ ncp ] ;

    for ( unsigned i = 0 ; i < ncp ; i++ ) {
      AtB[ i ] = new double[ 3 ] ;
    }

    comp_matrix_atb( A , image_pts , np , ncp , AtB ) ;

    /**
     * Solve the linear system "AtA Y = AtB[i]", where "AtB[i]" is the
     * i-th column  of "AtB",  for i=0,1,2. This  is done by  using LU
     * decomposition.
     */
    
    /**
     * Compute the LU decomposition of AtA.
     */
    lssolver::tLUdcmp LU( AtA , ncp );

    /**
     * Allocate  memory for the  solution vector,  "Y", of  the linear
     * system "AtA Y = AtB[i].
     *
     */

    double* Y = new double[ ncp ] ;

    for ( unsigned i = 0 ; i < 3 ; i++ ) {
      for ( unsigned j = 0 ; j < ncp ; j++ ) {
        Y[ j ] = AtB[ j ][ i ] ;
      }

      /**
       * Solve the linear system AtA Y = AtB[i].
       */
      LU.solve( Y , ncp ) ; 

      /**
       * Store  the  solution  in   the  i-th  column  of  the  matrix
       * "_ctrl_pts[j]".
       */
      for ( unsigned j = 0 ; j < ncp ; j++ ) {
        m_ctrl_pts[ j ][ i ] = Y[ j ] ;
      }
    }

    /**
     * Release memory
     */    
    if (A != 0) {
      for ( unsigned i = 0 ; i < np ; i++ ) {
        if ( A[ i ] != 0 ) {
          delete[] A[ i ] ;
        }
      }

      delete A;
    }

    if ( AtA != 0 ) {
      for ( unsigned i = 0 ; i < ncp ; i++ ) {
        if ( AtA[ i ] != 0 ) {
          delete[] AtA[ i ] ;
        }
      }

      delete AtA ;
    }

    if ( AtB != 0 ) {
      for ( unsigned i = 0 ; i < ncp ; i++ ) {
        if ( AtB[ i ] != 0 ) {
          delete[] AtB[ i ] ;
        }
      }

      delete AtB ;
    }

    delete[] Y ;

    return ;
  }


  /**
   * \fn tBezier( const tBezier& bz )
   *
   * \brief Clones  an instance of  this Bezier class (through  a deep
   * copy).
   *
   * \param bz The instance of the tBezier class to be cloned.
   */
  tBezier::tBezier( const tBezier& bz ) 
  {
    /**
     * Initialize private attributes.
     */

    set_bidegree_1st_index( bz.get_bidegree_1st_index() ) ;
    set_bidegree_2nd_index( bz.get_bidegree_2nd_index() ) ;

    set_aff_1st_point_x_coord( bz.get_aff_1st_point_x_coord() ) ;
    set_aff_1st_point_y_coord( bz.get_aff_1st_point_y_coord() ) ;

    set_aff_2nd_point_x_coord( bz.get_aff_2nd_point_x_coord() ) ;
    set_aff_2nd_point_y_coord( bz.get_aff_2nd_point_y_coord() ) ;

    /**
     * Allocate memory for the control points.
     */
    unsigned ncp = ( get_bidegree_1st_index() + 1 ) * 
                   ( get_bidegree_2nd_index() + 1 ) ;

    m_ctrl_pts = new double*[ 
			     ( get_bidegree_1st_index() + 1 ) * 
			     ( get_bidegree_2nd_index() + 1 ) 
			    ] ;

    for ( unsigned i = 0 ; i < ncp ; i++ ) {
      m_ctrl_pts[ i ] = new double[ 3 ] ;

      for ( unsigned j = 0 ; j < 3 ; j++ ) {
	m_ctrl_pts[ i ][ j ] = bz.m_ctrl_pts[ i ][ j ] ;
      }
    }

    return ;
  }


  /**
   * \fn ~tBezier()
   *
   * \brief Destroys an instance of this class.
   */
  tBezier::~tBezier()
  {
    if ( m_ctrl_pts != 0 ) {

      unsigned ncp = ( get_bidegree_1st_index() + 1 ) * 
	             ( get_bidegree_2nd_index() + 1 ) ;

      for ( unsigned i = 0 ; i < ncp ; i++ ) {
        if ( m_ctrl_pts[ i ] != 0 ) {
          delete[] m_ctrl_pts[ i ] ;
        }
      }

      delete[] m_ctrl_pts ;
    }

    return ;
  }


  /**
   * \fn void b( unsigned , unsigned , double& , double& , double& ) const
   *
   * \brief Get a specific control point defining this Bezier patch.
   *
   * \param i First index of the control point.
   * \param j Second index of the control point.
   * \param x First Cartesian coordinate of the control point.
   * \param y Second Caresian coordinate of the control point.
   * \param z Third Cartesian coordinate of the control point.
   * 
   */
  void
  tBezier::b(
             unsigned i ,
             unsigned j ,
	     double& x ,
	     double& y ,
	     double& z
            ) 
    const
  {
    assert( 
	   ( i <= get_bidegree_1st_index() ) && 
	   ( j <= get_bidegree_2nd_index() ) 
	  ) ;

    unsigned l = index( j , i ) ;

    x = m_ctrl_pts[ l ][ 0 ] ;
    y = m_ctrl_pts[ l ][ 1 ] ;
    z = m_ctrl_pts[ l ][ 2 ] ;

    return ;
  }


  /**
   * \fn void point( double u , double v , double& x , double& y , double& z ) const
   *
   * \brief Compute a point on this rectangular Bezier patch.
   *
   * \param u First Cartesian coordinate of the parameter point.
   * \param v Second Cartesian coordinate of the parameter point.
   * \param x First Cartesian coordinate of the point on the patch.
   * \param y Second Caresian coordinate of the point on the patch.
   * \param z Third Cartesian coordinate of the point on the patch.
   */
  void
  tBezier::point(
		 double u , 
		 double v ,
		 double& x ,
		 double& y ,
		 double& z 
		) 
    const
  {
    /**
     * Map the point to the affine frame [0,1].
     */
    double rx = get_aff_1st_point_x_coord() ;
    double ry = get_aff_1st_point_y_coord() ;
    double sx = get_aff_2nd_point_x_coord() ;
    double sy = get_aff_2nd_point_y_coord() ;
    
    /**
     * Compute all Bernstein polynomials of degree \var{_m}.
     */

    double uu = ( u - rx ) / ( sx - rx ) ;
    double vv = ( v - ry ) / ( sy - ry ) ;

    if ( fabs( uu ) <= 1e-15 ) {
      uu = 0 ;
    }
    else if ( fabs( 1 - uu ) <= 1e-15 ) {
      uu = 1 ;
    }

    if ( fabs( vv ) <= 1e-15 ) {
      vv = 0 ;
    }
    else if ( fabs( 1 - vv ) <= 1e-15 ) {
      vv = 1 ;
    }

    unsigned m = get_bidegree_1st_index() ;
    unsigned n = get_bidegree_2nd_index() ;

    std::vector< double > bu( m + 1 ) ;
    all_bernstein( m , uu , bu );

    /**
     * Compute all Bernstein polynomials of degree \var{_n}.
     */
    std::vector< double > bv( n + 1 ) ;
    all_bernstein( n , vv , bv );

    /**
     * Compute the image point, (x,y,z), of (u,v) on this patch.
     */
    x = 0 ;
    y = 0 ;
    z = 0 ;
    for ( unsigned j = 0 ; j <= n ; j++ ) {
      for ( unsigned i = 0 ; i <= m ; i++ ) {
        double buv = bu[ i ] * bv[ j ] ;
	double xaux ;
	double yaux ;
	double zaux ;
	b( i , j , xaux , yaux , zaux ) ;

	x += buv * xaux ;
	y += buv * yaux ;
	z += buv * zaux ;
      }
    }

    return ;
  }


  /**
   * \fn double bernstein( unsigned n , unsigned i , double u ) const
   *
   * \brief Computes the value  of i-th Bernstein polynomial od degree
   * n (B_(n,i)) at a given parameter value with respect to the affine
   * frame [0,1].
   *
   * \param n The degree of the Bernstein polynomial.
   * \param i The index identifying the Bernstein polynomial.
   * \param u A parameter value.
   *
   * \return The value of the i-th Bernstein polynomial of degree n at
   * the given parameter values.
   */
  double
  tBezier::bernstein(
		     unsigned n ,
		     unsigned i ,
		     double u
		    )
    const 
  {
    assert( i <= n ) ;

    std::vector< double > temp( n + 1 ) ;

    for ( unsigned j = 0 ; j <= n ; j++ ) {
      temp[ j ] = 0 ;
    }

    temp[ n - i ] = 1 ;

    double u1 = 1 - u ;

    for ( unsigned k = 1 ; k <= n ; k++ ) {
      for ( unsigned j = n ; j >= k ; j-- ) {
        temp[ j ] = ( u1 * temp[ j ] ) + ( u * temp[ j - 1 ] ) ;
      }
    }

    return temp[ n ] ;
  }


  /**
   * This function computes the  value of all Bernstein polynomials at
   * a given parameter value with respect to the affine frame [0,1].
   *
   * \param n The degree of the Bernstein polynomial.
   * \param u A parameter value.
   * \param b An array with  the value of all Bernstein polynomials at
   * the given parameter value.
   */
  void
  tBezier::all_bernstein(
			 unsigned n ,
			 double u ,
			 std::vector<double>& b
			)
    const
  {
    b[ 0 ] = 1 ;
    double u1 = 1 - u ;

    for ( unsigned j = 1 ; j <= n ; j++ ) {
      double saved = 0 ;
      for ( unsigned k = 0 ; k < j ;  k++ ) {
        double temp = b[ k ] ;
        b[ k ] = saved + ( u1 * temp ) ;
        saved = u * temp ;
      }

      b[ j ] = saved ;
    }

    return ;
  }

  /**
   * \fn void comp_bpoly_matrix( double**& a , double* param_pts , unsigned np ) const
   *
   * \brief Computes the pairwise product of all Bernstein polynomials
   * of  two degrees  at a  given set  of parameters  points.  The two
   * degrees are  the the first and  second index of  the bi-degree of
   * this patch.
   *
   * \param a A (np x (( m + 1) * ( n + 1)) matrix such that ( m , n )
   * is the bi-degree of this patch  and each element of the matrix is
   * a  product of  two  Bernstein polynomials  at  a given  parameter
   * point.
   * \param patch_pts An array of 2D coordinates of parameter points.
   * \param np The number of parameter points.
   */
  void
  tBezier::comp_bpoly_matrix(
			     double**& a ,
			     double* param_pts , 
			     unsigned np 
			    ) 
      const
  {
    /**
     * Compute the value of the Bernstein polynomials at the parameter
     * points.
     */
    double rx = get_aff_1st_point_x_coord() ;
    double ry = get_aff_1st_point_y_coord() ;
    double sx = get_aff_2nd_point_x_coord() ;
    double sy = get_aff_2nd_point_y_coord() ;

    for ( unsigned i = 0 ; i < np ; i++ ) {
      double x = ( param_pts[ ( i << 1 )     ] - rx ) / 
	( sx - rx ) ;

      double y = ( param_pts[ ( i << 1 ) + 1 ] - ry ) / 
	( sy - ry ) ;

      unsigned m = get_bidegree_1st_index() ;
      unsigned n = get_bidegree_2nd_index() ;

      std::vector< double > bu( m + 1 ) ;
      std::vector< double > bv( n + 1 ) ;

      all_bernstein( m , x , bu ) ;
      all_bernstein( n , y , bv ) ;

      for ( unsigned k = 0 ; k <= n ; ++k ) {
        for ( unsigned j = 0 ; j <= m ; ++j ) {
          unsigned l = index( k , j ) ;
          a[ i ][ l ] = bu[ j ] * bv[ k ] ;
        }
      }
    }

    return ;
  }


  /**
   * \fn void comp_matrix_ata( double** a , unsigned n , unsigned p , double**& ata ) const
   *
   * \brief Computes the product of  the transpose of a given matrix
   * with the matrix itself.
   *
   * \param a A (n x p) matrix.
   * \param n The number of rows of \var{a}.
   * \param p The number of columns of \var{a}.
   * \param ata A  (\var{p} x \var{p}) matrix equal  to the product of
   * the transpose of a given matrix with the matrix itself.
   */
  void
  tBezier::comp_matrix_ata(
			   double** a ,
			   unsigned n ,
			   unsigned p ,
			   double**& ata
			  )
    const
  { 
    for ( unsigned i = 0 ; i < p ; i++ ) {
      for ( unsigned k = 0 ; k < p ; k++ ) {
        ata[ i ][ k ] = 0 ;
        for ( unsigned j = 0 ; j < n ; j++ ) {
          ata[ i ][ k ] += a[ j ][ i ] * a[ j ][ k ] ;
        }
      }
    }

    return ;
  }


  /**
   * \fn void comp_matrix_atb( double** a , double* b , unsigned n , unsigned p , double**& atb ) const ;
   *
   * \brief Computes the product of the transpose of a given ( n x p )
   *  matrix "a" and a ( n x 3 ) matrix "b". The result is a ( p x 3 )
   *  matrix.
   *
   * \param a The first matrix.
   * \param b The second matrix.
   * \param n The number of rows of the given matrices.
   * \param p The number of columns of the first matrix.
   * \param atb  The matrix resulting  from the multiplication  of the
   * transpose of the first matrix by the second matrix.
   */
  void
  tBezier::comp_matrix_atb(
			   double** a ,
			   double* b ,
			   unsigned n ,
			   unsigned p ,
			   double**& atb
			  )
    const
  {
    /**
     * Compute the  product (at x b),  where "at" is  the transpose of
     * "a".
     */
    for ( unsigned i = 0 ; i < p ; i++ ) {
      atb[ i ] = new double[ 3 ] ;
      for ( unsigned k = 0 ; k < 3 ; k++ ) {
        atb[ i ][ k ] = 0 ;
        for ( unsigned j = 0 ; j < n ; j++ ) {
          atb[ i ][ k ] += a[ j ][ i ] * b[ ( 3 * j ) + k ] ;
        }
      }
    }

    return ;
  }

}

/** @} */ //end of group class.
