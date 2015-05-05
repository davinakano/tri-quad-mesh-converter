/** 
 * \file bezier.h
 *  
 * \brief  Definition   of  the   class  tBezier  that   represents  a
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

#ifndef tBEZIER_H
#define tBEZIER_H

#include <vector>           // std::vector
#include <cassert>          // assert

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
   * \class tBezier
   *
   * \brief This class represents a rectangular B&eacute;zier patch in
   * 3D.
   *
   */
  class tBezier {
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn tBezier( double* param_pts , double* image_pts , unsigned np , unsigned m , unsigned n , double rx , double sx , double ry , double sy )
     *
     * \brief Creates  an instance of  this Bezier class by  solving a
     * least  squares  fitting problem.   In  particular, the  control
     * points of  the patch are computed  as the solution  of a linear
     * system of normal  equations. To set up this  system, the method
     * receives  a set  of  parameter points  and their  corresponding
     * image points in 3D. The  patch is built so that it approximates
     * the latter points.
     *
     * \param param_pts an  array with the 2D coordinates  of a set of
     * parameter points.
     * \param image_pts An  array with the 3D coordinates  of a set of
     * 3D points, which are the image of the paramenter points under a
     * possibly unknown function.
     * \param np The number of parameter points.
     * \param n The first degree of this rectangular Bezier patch.
     * \param m The second degree of this rectangular Bezier patch.
     * \param rx The first Cartesian  coordinate of the first point of
     * the affine frame w.r.t. which this surface patch is defined.
     * \param  sx The  second coordinate  of the  first point  of the
     * affine frame w.r.t. which this surface patch is defined.
     * \param ry The first Cartesian coordinate of the second point of
     * the affine frame w.r.t. which this surface patch is defined.
     * \param sy  The second Cartesian coordinate of  the second point
     * of the affine frame w.r.t. which this surface patch is defined.
     *
     */
    tBezier( 
	    double* param_pts ,
	    double* image_pts ,
	    unsigned np ,
	    unsigned m ,
	    unsigned n ,
	    double rx ,
	    double sx ,
	    double ry ,
	    double sy
	   ) ;


    /**
     * \fn tBezier( const tBezier& bz )
     *
     * \brief Clones an instance of  this Bezier class (through a deep
     * copy).
     *
     * \param bz The instance of the tBezier class to be cloned.
     */
    tBezier( const tBezier& bz ) ;


    /**
     * \fn ~tBezier()
     *
     * \brief Destroys an instance of this class.
     */
    ~tBezier() ;


    /**
     * \fn inline unsigned get_bidegree_1st_index() const
     *
     * \brief Returns the first index of the bi-degree of this patch.
     *
     * \return The first index of the bi-degree of this patch.
     */
    inline unsigned get_bidegree_1st_index() const
    {
      return m_m ;
    }


    /**
     * \fn inline unsigned get_bidegree_2nd_index() const
     *
     * \brief Returns the second index of the bi-degree of this patch.
     *
     * \return The second index of the bi-degree of this patch.
     */
    inline unsigned get_bidegree_2nd_index() const
    {
      return m_n ;
    }


    /**
     * \fn inline double get_aff_1st_point_x_coord() const
     *
     * \brief  Returns the  first  Cartesian coordinate  of the  first
     * point of the  affine frame with respect to  which this patch is
     * defined.
     *
     * \return The  first Cartesian coordinate  of the first  point of
     * the affine frame with respect to which this patch is defined.
     */
    inline double get_aff_1st_point_x_coord() const
    {
      return m_rx ;
    }


    /**
     * \fn inline unsigned get_aff_1st_point_y_coord() const
     *
     * \brief  Returns the  second Cartesian  coordinate of  the first
     * point of the  affine frame with respect to  which this patch is
     * defined.
     *
     * \return The  second Cartesian coordinate of the  first point of
     * the affine frame with respect to which this patch is defined.
     */
    inline double get_aff_1st_point_y_coord() const
    {
      return m_ry ;
    }


    /**
     * \fn inline unsigned get_aff_2nd_point_x_coord() const
     *
     * \brief  Returns the  first Cartesian  coordinate of  the second
     * point of the  affine frame with respect to  which this patch is
     * defined.
     *
     * \return The  first Cartesian coordinate of the  second point of
     * the affine frame with respect to which this patch is defined.
     */
    inline double get_aff_2nd_point_x_coord() const
    {
      return m_sx ;
    }


    /**
     * \fn inline double get_aff_2nd_point_y_coord() const
     *
     * \brief Returns  the second  Cartesian coordinate of  the second
     * point of the  affine frame with respect to  which this patch is
     * defined.
     *
     * \return The second Cartesian  coordinate of the second point of
     * the affine frame with respect to which this patch is defined.
     */
    inline double get_aff_2nd_point_y_coord() const
    {
      return m_sy ;
    }


    /**
     * \fn void b( unsigned i , unsigned j , double& x , double& y , double& z ) const
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
    void b( unsigned i , unsigned j , double& x , double& y , double& z )
      const ;


    /**
     * \fn void point( double , double , double& x , double& y , double& z ) const
     *
     * \brief Compute a point on this rectangular Bezier patch.
     *
     * \param u First Cartesian coordinate of the parameter point.
     * \param v Second Cartesian coordinate of the parameter point.
     * \param x First Cartesian coordinate of the point on the patch.
     * \param y Second Caresian coordinate of the point on the patch.
     * \param z Third Cartesian coordinate of the point on the patch.
     */
    void point( double , double , double& x , double& y , double& z ) const ;


  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------


    /**
     * \fn inline void set_bidegree_1st_index( unsigned m )
     *
     * \brief Assigns a  value to the first index  of the bi-degree of
     * this patch.
     *
     * \param m  The value to  be assigned to  the first index  of the
     * bi-degree of this patch.
     */
    inline void set_bidegree_1st_index( unsigned m )
    {
      assert( m > 0 ) ;

      m_m = m ;
    }


    /**
     * \fn inline void set_bidegree_2nd_index( unsigned n )
     *
     * \brief Assigns a value to  the second index of the bi-degree of
     * this patch.
     *
     * \param n  The value to be  assigned to the second  index of the
     * bi-degree of this patch.
     */
    inline void set_bidegree_2nd_index( unsigned n )
    {
      assert( n > 0 ) ;

      m_n = n ;
    }


   /**
     * \fn inline void set_aff_1st_point_x_coord( double rx )
     *
     * \brief Assigns a value to the first Cartesian coordinate of the
     * first  point of  the affine  frame with  respect to  which this
     * patch is defined.
     *
     * \param  rx  Value  to   be  assigned  to  the  first  Cartesian
     * coordinate of the first point  of the affine frame with respect
     * to which this patch is defined.
     */
    inline void set_aff_1st_point_x_coord( double rx )
    {
      m_rx = rx ;
    }


    /**
     * \fn inline void set_aff_1st_point_y_coord( double ry )
     *
     * \brief Assigns  a value to  the second Cartesian  coordinate of
     * the first point of the  affine frame with respect to which this
     * patch is defined.
     *
     * \param  ry  Value  to  be  assigned  to  the  second  Cartesian
     * coordinate of the first point  of the affine frame with respect
     * to which this patch is defined.
     */
    inline void set_aff_1st_point_y_coord( double ry )
    {
      m_ry = ry ;
    }


    /**
     * \fn inline void set_aff_2nd_point_x_coord( double sx )
     *
     * \brief Assigns a value to the first Cartesian coordinate of the
     * second point  of the  affine frame with  respect to  which this
     * patch is defined.
     *
     * \return Assigns  a value to  the first Cartesian  coordinate of
     * the second point of the affine frame with respect to which this
     * patch is defined.
     */
    inline void set_aff_2nd_point_x_coord( double sx )
    {
      m_sx = sx ;
    }


    /**
     * \fn inline void set_aff_2nd_point_y_coord( double sy )
     *
     * \brief Assigns  a value to  the second Cartesian  coordinate of
     * the second point of the affine frame with respect to which this
     * patch is defined.
     *
     * \param  sy  Value  to  be  assigned  to  the  second  Cartesian
     * coordinate of the second point of the affine frame with respect
     * to which this patch is defined.
     */
    inline void set_aff_2nd_point_y_coord( double sy )
    {
      m_sy = sy ;
    }


    /**
     * \fn inline unsigned index( unsigned i , unsigned j ) const
     *
     * \brief Computes the "linear" index of a control point.
     *
     * \param i First index of the control point.
     * \param j Second index of the control point.
     *
     * \return The linear index of a control point.
     */
    inline unsigned index( unsigned i , unsigned j ) const
    {
      return ( i * ( m_m + 1 ) ) + j ;
    }


    /**
     * \fn double bernstein( unsigned n , unsigned i , double u ) const
     *
     * \brief  Computes  the value  of  i-th  Bernstein polynomial  of
     * degree n (B_(n,i))  at a given parameter value  with respect to
     * the affine frame [0,1].
     *
     * \param n The degree of the Bernstein polynomial.
     * \param i The index identifying the Bernstein polynomial.
     * \param u A parameter value.
     *
     * \return The value of the  i-th Bernstein polynomial of degree n
     * at the given parameter values.
     */
    double bernstein( unsigned , unsigned , double ) const ;


    /**
     * \fn void all_bernstein( unsigned n , double u , std::vector<double>& b ) const
     *
     * \brief  Computes  the value  of  all  Bernstein polynomials  of
     * degree n at a given  parameter value with respect to the affine
     * frame [0,1].
     *
     * \param n The degree of the Bernstein polynomials.
     * \param u A parameter value.
     * \param b An  array with the value of  all Bernstein polynomials
     * at the given parameter value.
     */
    void all_bernstein( unsigned n , double u , 
			std::vector<double>& b ) const ;


    /**
     * \fn void comp_bpoly_matrix( double**& a , double* param_pts , unsigned np ) const
     *
     * \brief   Computes  the  pairwise   product  of   all  Bernstein
     * polynomials  of  two  degrees  at  a given  set  of  parameters
     * points. The two  degrees are the the first  and second index of
     * the bi-degree of this patch.
     *
     * \param a A (np x (( m + 1) * ( n + 1)) matrix such that ( m , n
     * ) is the bi-degree of this patch and each element of the matrix
     * is a product of two  Bernstein polynomials at a given parameter
     * point.
     * \param  patch_pts  An  array  of 2D  coordinates  of  parameter
     * points.
     * \param np The number of parameter points.
     */
    void comp_bpoly_matrix( double**& a , double* param_pts , unsigned np ) 
      const ;


    /**
     * \fn void comp_matrix_ata( double** a , unsigned n , unsigned p , double**& ata ) const
     *
     * \brief Computes the product of  the transpose of a given matrix
     * with the matrix itself.
     *
     * \param a A (n x p) matrix.
     * \param n The number of rows of \var{a}.
     * \param p The number of columns of \var{a}.
     * \param ata A (\var{p} x \var{p}) matrix equal to the product of
     * the transpose of a given matrix with the matrix itself.
     */
    void comp_matrix_ata( double** a , unsigned n , unsigned p , double**& ata )
      const ;


    /**
     * \fn void comp_matrix_atb( double** a , double* b , unsigned n , unsigned p , double**& atb ) const ;
     *
     * \brief Computes the product of the transpose of a given ( n x p
     *  ) matrix "a" and a ( n x 3 ) matrix "b". The result is a ( p x
     *  3 ) matrix.
     *
     * \param a The first matrix.
     * \param b The second matrix.
     * \param n The number of rows of the given matrices.
     * \param p The number of columns of the first matrix.
     * \param atb The matrix  resulting from the multiplication of the
     * transpose of the first matrix by the second matrix.
     */
    void comp_matrix_atb( double** a , double* b , unsigned n , unsigned p , double**& atb ) const ;


    // ---------------------------------------------------------------
    //
    // Private data members
    //
    // ---------------------------------------------------------------    

    /**
     * Bi-degree of the Bezier patch.
     */
    unsigned m_m ;   ///< First bi-degree index.
    unsigned m_n ;   ///< Second bi-degree index.


    /**
     * Affine  frame  with respect  to  which  this  surface patch  is
     * defined.
     */
    double m_rx ;    ///< First Cartesian coordinate of the first point of the affine frame.
    double m_sx ;    ///< Second Cartesian coordinate of the first point of the affine frame.
    double m_ry ;    ///< First Cartesian coordinate of the second point of the affine frame.
    double m_sy ;    ///< Second Cartesian coordinate of the second point of the affine frame.

    /**     
     * Control points.
     */
    double** m_ctrl_pts ;   ///< Matrix of control points. 

  } ;

}

/** @} */ //end of group class.

#endif   // BEZIER_H
