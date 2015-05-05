/** 
 * \file floaterpar.cpp
 *  
 * \brief  Implementation  of  class  FloaterPar, which  implements  a
 * method for parametrizing a  disk-like surface patch in \f$R^3\f$ to
 * \f$R^2\f$.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * mfsiqueira at gmail (dot) com
 *
 * \version 1.0
 * \date March 2011
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include "floaterpar.hpp"     // FloaterPar

#include "geodesics/geometric.hpp"      // common::Geometric::distance_3d()
#include "numerics/ludcmp.hpp"         // numerics::LUdcmp

#include <suitesparse/umfpack.h>

#include <cmath>              // HUGE_VAL, acos, fabs
#include <climits>            // UNIT_MAX
#include <cassert>            // assert
#include <list>               // std::list

#include <iostream>
#include <algorithm>

/** 
 * \defgroup FLOATERPARNameSpace Namespace floaterpar.
 * @{
 */

/**
 * \namespace floaterpar
 *
 * \brief The namespace floaterpar  contains all classes related to an
 * implementation  of  the mesh  parametrization  method developed  by
 * Floater.
 */

namespace floaterpar {


/**
   * \fn void FloaterPar::initialize_vertex_attributes( double* tv ) 
   *
   * \brief A method for creating vertex attributes.
   *
   * \param tv  An array with the  t values associated  with the patch
   * vertices.
   */
void
FloaterPar::initialize_vertex_attributes( double* tv )
{
    /*
     * For each  patch vertex, assign  a t value and  dummy parametric
     * coordinates.
     */
    for ( unsigned i = 0 ; i < patch()->get_number_of_vertices() ; ++i ) {
        /*
       * Set the t value using the given t value.
       */
        patch()->get_vertex_attributes( i ).set_tval( tv[ i ] ) ;

        /*
       * Set dummy coordinate values.
       */
        patch()->get_vertex_attributes( i ).set_ucoord( HUGE_VAL ) ;
        patch()->get_vertex_attributes( i ).set_vcoord( HUGE_VAL ) ;
    }

    return ;
}


/**
   * \fn void FloaterPar::parametrize() 
   *
   * \brief computes a parametrization of the associated surface patch
   * using  Floater's  method  as  described  in  the  paper  "M.   S.
   * Floater,  Parametrization  and  smooth approximation  of  surface
   * triangulations,  Computer  Aided   Geometric  Design  14  (1997),
   * 231-250."
   */
void
FloaterPar::parametrize()
{
    /*
     * Separate inner from boundary vertices.
     */
    std::vector< unsigned > iv ;
    std::vector< unsigned > bv ;

    separate_vertices( iv , bv ) ;

    /*
     * Compute the parameter coordinates of the boundary vertices.
     */
    comp_bv_pcoords( bv ) ;

    /*
     * Build  a table to  index the  vertices in  the arrays  "iv" and
     * "bv".
     */
    const unsigned K = iv.size() ;
    const unsigned N = bv.size() ;

    if( K==0 ) {
        _is_done = true ;
        return;
    }

    std::vector< unsigned > vtab( K + N ) ;

    for ( unsigned i = 0 ; i < K ; i++ ) {
        vtab[ iv[ i ] ] = i ;
    }

    for ( unsigned i = 0 ; i < N ; i++ ) {
        vtab[ bv[ i ] ] = i ;
    }

    /*
     * Compute the  mean value coordinates  of each inner  vertex with
     * respect to the  vertices in its star. We  use those coordinates
     * to create two systems of  K linear equations, AX=U and AX=V, in
     * K unknowns each.
     */

    //    double** A = new double*[ K ] ;
    //    for ( unsigned i = 0 ; i < K ; i++ ) {
    //        A[ i ] = new double[ K ] ;

    //        for ( unsigned j = 0 ; j < K ; j++ ) {
    //            A[ i ][ j ] = 0 ;
    //        }
    //    }
    //boost::numeric::ublas::matrix<double> A( K, K );
    //A *= 0;
    std::vector<int> Ap(K+1);
    std::vector<int> Ai;
    std::vector<double> Ax;
    Ap[0] = 0;


    std::vector<double> U(K);
    //    double* U = new double[ K ] ;
    for ( unsigned i = 0 ; i < K ; i++ ) {
        U[ i ] = 0 ;
    }
    //boost::numeric::ublas::vector<double> U( K );
    //U *= 0;


    std::vector<double> V(K);
    //    double* V = new double[ K ] ;
    for ( unsigned i = 0 ; i < K ; i++ ) {
        V[ i ] = 0 ;
    }
    //boost::numeric::ublas::vector<double> V( K );

    std::vector< std::pair< std::pair<int,int>, double > > A;

    for ( unsigned i = 0 ; i < K ; i++ ) {
        /*
       * Get the star of the i-th inner vertex.
       */
        std::vector< unsigned > star ;
        patch()->get_vertex_star( iv[ i ] , star ) ;

        comp_mv_coord(
                    iv[ i ] ,
                    star ,
                    vtab ,
                    A, i,
                    U[i] ,
                    V[i]
                    ) ;

        //A( i, i ) = 1 ;
        A.push_back( std::pair< std::pair<int,int>, double >(std::pair<int,int>(i,i), 1.0 ) );
    }


    std::sort( A.begin(), A.end() );
    int nz = A.size();
    int col_anterior = 0;
    int col, lin;
    for(int i=0; i<nz; i++) {

        col = A[i].first.first;
        lin = A[i].first.second;

        Ai.push_back(lin);

        if( col != col_anterior ) {
            for( int j=col_anterior; j<col; j++){
                Ap[j+1] = i;
            }
            col_anterior = col;
        }

        Ax.push_back( A[i].second );
    }
    for( int j=col_anterior; j<=col; j++){
        Ap[j+1] = nz;
    }
    A.reserve(0);


    /*
     * Solve the linear systems using LU decomposition.
     *
     */
    //numerics::LUdcmp lu( A , K ) ;
    //lu.solve( U , K ) ;
    //lu.solve( V , K ) ;
    //boost::numeric::ublas::permutation_matrix<std::size_t> pm(A.size1());
    //int ret = boost::numeric::ublas::lu_factorize(A, pm);
    //assert( ret == 0 );
    //boost::numeric::ublas::lu_substitute(A,pm, U);
    //boost::numeric::ublas::lu_substitute(A,pm, V);

    double *null = (double *) NULL ;
    int i ;
    void *Symbolic, *Numeric ;
    (void) umfpack_di_symbolic (K, K, &(Ap[0]), &(Ai[0]), &(Ax[0]), &Symbolic, null, null) ;
    (void) umfpack_di_numeric (&(Ap[0]), &(Ai[0]), &(Ax[0]), Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;

    double *Ux = new double[K];
    double *Vx = new double[K];

    (void) umfpack_di_solve (UMFPACK_A, &(Ap[0]), &(Ai[0]), &(Ax[0]), Ux, &(U[0]), Numeric, null, null) ;
    (void) umfpack_di_solve (UMFPACK_A, &(Ap[0]), &(Ai[0]), &(Ax[0]), Vx, &(V[0]), Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;


    /*
             * Update the parameter coordinates of the inner vertices.
             */
    for ( unsigned i = 0 ; i < K ; i++ ) {
        patch()->get_vertex_attributes( iv[ i ] ).set_ucoord( Ux[ i ] ) ;
        patch()->get_vertex_attributes( iv[ i ] ).set_vcoord( Vx[ i ] ) ;
    }

    delete [] Ux;
    delete [] Vx;

    /*
     * Release dynamically allocated memory.
     *
     */
    //    if ( A != 0 ) {
    //        for ( unsigned i = 0 ; i < K ; i++ ) {
    //            if ( A[ i ] != 0 ) {
    //                delete[] A[ i ] ;
    //            }
    //        }

    //        delete[] A ;
    //    }

    //    if ( U != 0 ) {
    //        delete[] U ;
    //    }

    //    if ( V != 0 ) {
    //        delete[] V ;
    //    }

    /*
     * If we reach this point, a parametrization has been computed.
     */
    _is_done = true ;

    return ;
}


/**
   * \fn void FloaterPar::sample( unsigned nv , double* param_pts , double*& image_pts )
   *
   * \brief Sample the parametrization  associated with this object to
   * create a  set of  3d points,  which are the  image of  the sample
   * points.
   *
   * \param nv The number of sampled points.
   * \param  param_pts   The  Cartesian  coordinates   of  the  sample
   * (parameter) points.
   * \param image_pts  The Cartesian coordinates  of the image  of the
   * parameter points.
   */
void
FloaterPar::sample(
    unsigned nv ,
    double* param_pts ,
    double*& image_pts
    )
{
    /*
     * Make sure  we have  a parametrization. If  we don't,  thrown an
     * exception.
     */
    assert( _is_done ) ;

    /*
     * Compute a bounding box for each parametrization triangle.
     */
    std::vector< BoundingBox > bbset( get_number_of_faces() ) ;
    compute_bounding_box( bbset ) ;

    /*
     * Compute the image of the sample points.
     */
    std::vector< double > pset ;
    compute_image_points( nv , param_pts , bbset , pset ) ;

    /*
     * Select vertices of the star parametrization.
     */
    // select_parametrization_points( pset1 , pset2 ) ;

#ifdef DEBUGMODE
    assert( ( pset.size() % 3 ) ==  0 ) ;
    assert( nv == ( pset.size() / 3 ) ) ;
#endif

    /*
     * Copy the sampled points to the parameter variables.
     */
    image_pts = new double[ pset.size() ] ;

    for ( unsigned i = 0 ; i < pset.size() ; i++ ) {
        image_pts[ i ] = pset[ i ] ;
    }

    return ;
}


/**
   * \fn void FloaterPar::separate_vertices( std::vector< unsigned >& iv , std::vector< unsigned >& bv )
   *
   * \brief  Separates  inner  from  boundary  vertices.   During  the
   * separation process, the consistency of the t-values are checked.
   *
   * \param iv A reference for an array with the inner vertices.
   * \param bv A reference for an array with the boundary vertices.
   */
void
FloaterPar::separate_vertices(
    std::vector<unsigned>& iv ,
    std::vector<unsigned>& bv
    )
{
    //
    // Separates vertices using their t-values.
    //

    /* Create a counter for vertices with t-value equals to zero. */
    unsigned nbv = 0 ;

    /* Create a marker for the vertex with t-value equals to zero. */
    unsigned fbv = patch()->get_number_of_vertices() ;

    for ( unsigned i = 0 ; i < patch()->get_number_of_vertices() ; ++i ) {
        /*
       * If the i-th vertex is an  inner one, then its t-value must be
       * -1. If  not, we  have an inconsistency,  and an  exception is
       * thrown.
       */
        if ( is_inner_vertex( i ) ) {

            assert( patch()->get_vertex_attributes( i ).get_tval() == -1 ) ;

            /*
         * Store "i" into the array of inner vertices.
         */
            iv.push_back( i );
        }
        else {
            /*
         * Make sure the t-value is valid.
         */
            assert( patch()->get_vertex_attributes( i ).get_tval() >= 0 ) ;

            /* Increment the counter of boundary vertices. */
            ++nbv ;

            /*
         * If the t-value of the i-th vertex is 0, keep track of it. 
         */
            if ( patch()->get_vertex_attributes( i ).get_tval() == 0 ) {
                fbv = i ;
            }
        }
    }

    /*
     * There must be at least three boundary vertices. 
     */
    assert( nbv >= 3 ) ;

    /*
     * There must exist at least one vertex with t-value equals to 0.
     */
    assert( fbv < patch()->get_number_of_vertices() ) ;

#ifdef DEBUGMODE
    assert( patch()->get_vertex_attributes( fbv ).get_tval() == 0 ) ;
#endif

    /*
     * Store "fbv" into the array of boundary vertices.
     */
    bv.push_back( fbv ) ;

    /* Create an auxiliary variable to store "fbv". */
    unsigned vp = fbv ;

    /** Get the t-value of "vp". */
    double tv = patch()->get_vertex_attributes( vp ).get_tval() ;

    /* Initialize a boundary vertex counter */
    unsigned nbvaux = 0 ;

    /* Get a half-edge incident with "vp". */
    unsigned he1 = patch()->get_vertex_edge( vp ) ;

    /* Counter for the number of distinct integer t-values. */
    unsigned h = 0 ;

    /*
     * We assume that "fbv" is a boundary vertex. Then, we walk on the
     * boundary from "fbv", in a counterclockwise manner, visiting all
     * of its vertices.
     */
    do {
        /*
       * Find the  (unique) boundary half-edge whose  origin vertex is
       * "vp".
       */
        if ( !patch()->is_boundary_edge( he1 ) ) {
            /*
         * Get the next half-edge in a clockwise traverse of the cycle
         * of half-edges whose origin is "fbv".
         */
            unsigned he2 = patch()->get_next( patch()->get_mate( he1 ) ) ;

            /*
         * If  "he2"  is  not  a  boundary  half-edge,  then  we  keep
         * traversing the cycle.  If we get back to "he1" then "vp" is
         * not  a  boundary  vertex,  which  means that  the  data  is
         * inconsistent.
         */
            while ( !patch()->is_boundary_edge( he2 ) && ( he2 != he1 ) ) {
                /*
           * Get  the next half-edge  in a  clockwise traverse  of the
           * cycle of half-edges whose origin is "vp".
           */
                he2 = patch()->get_next( patch()->get_mate( he2 ) ) ;

#ifdef DEBUGMODE
                assert( patch()->get_origin_vertex( he2 ) == vp ) ;
#endif
            }

            /* If "he1 == he2" we have got an inconsistency! */
            assert( he2 != he1 ) ;

            he1 = he2 ;
        }

        /*
       * Get the vertex that follows "vp" along the boundary.
       */
        he1 = patch()->get_next( he1 ) ;

        vp = patch()->get_origin_vertex( he1 ) ;

        /*
       * Check if  the t-value of the  next vertex is  larger than the
       * previous one.   If not, we have an  inconsistency, except for
       * the case  in which  the current vertex  ( "vp" )  has t-value
       * equals to 0.
       */
        double tvaux = patch()->get_vertex_attributes( vp ).get_tval() ;

        assert( ( tv < tvaux ) || ( ( tv != 0 ) && ( tvaux == 0 ) ) ) ;

        /*
       * Check if the floor of both t-values changed. If so, we are on
       * another boundary curve of the  patch, and thus "h" is updated
       * accordingly.
       */
        if ( unsigned( tvaux ) == ( unsigned( tv ) + 1 ) ) {
            ++h;
        }

        tv = tvaux ;

        /**
       * Store "vp" into the array of boundary vertices.
       */
        if ( vp != fbv ) {
            bv.push_back( vp ) ;
        }

        /** Increment the boundary vertex counter. */
        ++nbvaux ;
    }
    while ( vp != fbv ) ;

    /*
     * The boundary  vertex counter must equal the  amount of boundary
     * vertices.
     */
    assert( nbv == nbvaux ) ;

    /**
     * Check if  t-values are all  distinct and contains  all integers
     * from 0 to N, where N  is the floor of the maximum t-value among
     * all t-values.
     */
    assert( h == unsigned( patch()->get_vertex_attributes(
                               bv[ bv.size() - 1 ] ).get_tval() ) ) ;

    return ;
}


/**
   * \fn void FloaterPar::comp_bv_pcoords( const std::vector< unsigned >& bv ) const
   *
   * \brief  Computes  the   parameter  coordinates  of  the  boundary
   * vertices.
   *
   * \param bv An array with the boundary vertices.
   */  
void
FloaterPar::comp_bv_pcoords( const std::vector< unsigned >& bv ) const
{
#ifdef DEBUGMODE
    assert ( bv.size() >= 3 ) ;
#endif

    /**
     * Get  the number  of  boundary vertices  of the  parametrization
     * domain,  i.e.,   the  n-gon  over   which  the  mesh   will  be
     * parametrized.
     */
    unsigned vp = bv[ bv.size() - 1 ] ;
    unsigned ns = unsigned( patch()->get_vertex_attributes( vp ).get_tval() ) + 1 ;

#ifdef DEBUGMODE
    assert ( ns >= 3 ) ;
#endif

    /*
     * Compute  the  length  of  the  polygonal line  defined  by  the
     * (single) boundary  of the simplicial surface,  and also compute
     * the  length  between  each  two  consecutive  vertices  of  the
     * polygonal line  itself. We assume  that this polygonal  line is
     * simple and closed (i.e., a Jordan curve). Otherwise, what we do
     * here does not make any sense.
     *
     * From now on, we refer to this boundary as "polygonal line".
     */

    /*
     * The element  sm[ i ] will  receive the length  of the polygonal
     * line to  be mapped  to the i-th  side of the  parameter domain,
     * i.e., a n-gon.
     */
    std::vector< double > sm( ns ) ;

    /*
     * The element df[  i ][ j ] will receive the  length of the first
     * "j" segments  of the  polygonal line to  be mapped to  the i-th
     * side of the parameter domain, i.e., a n-gon.
     */
    std::vector< std::vector< double > > df( ns ) ;

    /* Auxiliary variables to store vertex coordinates. */

    /*
     * They are initialized to avoid annoying compiler warnings.
     */
    double x1 = 0 ;
    double y1 = 0 ;
    double z1 = 0 ;

    /*
     * The  variable "h" will  store the  index in  "bv" of  the first
     * vertex of  the polygonal line  whose length is  being currently
     * computed.
     */
    unsigned h = 0 ;

    /*
     * The variable  "i" is a counter  for the number of  sides of the
     * parameter domain, i.e., a n-gon.
     */
    for ( unsigned i = 0 ; i < ns ; i++ ) {
        /*
       * Get  the 3D  coordinates  of  the first  vertex  of the  i-th
       * polygonal line.
       */
        patch()->get_vertex_coords( bv[ h ] , x1 , y1 , z1 ) ;

        /*
       * Since "bv[  h ]"  is the first  vertex, initialize "sm[  i ]"
       * with "0".
       */
        sm[ i ] = 0 ;

        /*
       * The distance between the first vertex and itself is 0.
       */
        df[ i ].push_back( sm[ i ] ) ;

        /*
       * Compute the  length between each two  consecutive vertices of
       * the i-th segment.  The t-values of those vertices  must be in
       * the range [ double(i) , double(i+1) ].
       */

        /*
       * Get the index "j" of  the second vertex of the i-th polygonal
       * line.
       */
        unsigned j = h + 1 ;

        bool endloop = false ;
        while ( j < bv.size() && !endloop ) {
            /*
         * Get the t-value of the j-th vertex.
         */
            double tv = patch()->get_vertex_attributes( bv[ j ] ).get_tval() ;

            if ( tv <= double( i + 1 ) ) {
                /*
           * Get the  3D coordinates of  the first vertex of  the j-th
           * polygonal line.
           *
           */
                double x2 = 0 ;
                double y2 = 0 ;
                double z2 = 0 ;
                patch()->get_vertex_coords( bv[ j ] , x2 , y2 , z2 ) ;

                /*
           * Compute  the coordinates  of  a vector  defined from  the
           * previous (the "j-1"-vertex) to  the current vertex of the
           * i-th polygonal line.
           */
                sm[ i ] += common::Geometric::distance_3d(
                            x1 ,
                            y1 ,
                            z1 ,
                            x2 ,
                            y2 ,
                            z2
                            ) ;

                /*
           * Also, store  the length of  the polygonal line up  to its
           * "j"-th vertex (the current one).
           */
                df[ i ].push_back( sm[ i ] ) ;

                /*
           * Update the  variables that  store the coordinates  of the
           * "previous" vertex,  so that we  can use them in  the next
           * loop iteration.
           */
                x1 = x2 ;
                y1 = y2 ;
                z1 = z2 ;

                /*
           * Increment the counter of  vertices of the current polygonal
           * line.
           */
                ++j ;
            }
            else {
                endloop = true ;
            }
        }

        /*
       * The "j"-th vertex is the  second vertex of the next polygonal
       * line (if  any), so we must  assign "j-1" to "h",  so that "h"
       * keeps the  index of  the first vertex  of the  next polygonal
       * line.
       */
        h = j - 1 ;
    }

    /*
     * Compute  the  length of  the  very  last  segment of  the  last
     * polygonal line.
     */
    double x2 = 0 ;
    double y2 = 0 ;
    double z2 = 0 ;
    patch()->get_vertex_coords( bv[ 0 ] , x2 , y2 , z2 ) ;

    sm[ ns - 1 ] += common::Geometric::distance_3d(
                x1 ,
                y1 ,
                z1 ,
                x2 ,
                y2 ,
                z2
                ) ;

    df[ ns - 1 ].push_back( sm[ ns - 1 ] ) ;

    /*
     * Finally, we  compute the coordinates of  the boundary vertices.
     * The  vertices will  be mapped  to  the edges  of the  parameter
     * domain (the n-gon). The j-th  vertex of the i-th polygonal line
     * is mapped to  the i-th side of the n-gon. If  "p0" and "p1" are
     * the end points of the n-gon  i-th side, then the j-th vertex is
     * mapped to the point
     *
     * (1 - t) * p0 + t * p1
     *
     * where "t" is the length of the first "j-1" segments of the i-th
     * polygonal  line  divided  by  the  total  length  of  the  i-th
     * polygonal line.
     *
     */

    /* Compute the vertices of the n-gon. */

    std::vector< double > u( ns + 1 ) ;
    std::vector< double > v( ns + 1 ) ;

    /*
     * If you would like to have [ ( 1 , 0 ) , ( 0 , 1 ) , ( 0 , 0 ) ]
     * as the affine frame for the parametrization domain, just change
     * the if-else below:  uncomment the if below and  comment the one
     * following it.
     *
     *
     * if ( ns == 3 ) {
     *
     *   We are using [ ( 0 , 0 ) , (  1 , 0 ) , ( 0 , 1 ) ] as affine
     *   frame.
     *
     *   u[ 0 ] = u[ 2 ] = u[ 3 ] = v[ 0 ] = v[ 1 ] = v[ 3 ] = 0 ;
     *   u[ 1 ] = v[ 2 ] = 1 ;
     * }
     *
     */

    if ( ns == 3 ) {
        u[ 0 ] = u[ 3 ] = v[ 1 ] = 1 ;
        u[ 1 ] = u[ 2 ] = v[ 0 ] = v[ 2 ] = v[ 3 ] = 0 ;
    } else if ( ns == 4 ) {
        u[ 0 ] = 1;
        v[ 0 ] = 0;
        u[ 1 ] = 1;
        v[ 1 ] = 1;
        u[ 2 ] = 0;
        v[ 2 ] = 1;
        u[ 3 ] = 0;
        v[ 3 ] = 0;
        u[ 4 ] = 1;
        v[ 4 ] = 0;
    } else {
        const double ang = ( 2 * acos( -1 ) ) / ns ;

        u[ 0 ] = u[ ns ] = 1 ;
        v[ 0 ] = v[ ns ] = 0 ;

        for ( unsigned i = 1 ; i < ns ; i++ ) {
            u[ i ] = cos( i * ang ) ;
            v[ i ] = sin( i * ang ) ;
        }
    }

    /*
     * Map  the vertices of  the polygonal  line to  the sides  of the
     * n-gon.
     */
    h = 0 ;
    for ( unsigned i = 0 ; i < ns ; i++ ) {
        /* Get the number of vertices to be mapped to the i-th side. */
        unsigned max = df[ i ].size() ;

        /*
       * Map the j-th vertex of the i-th segment of the polygonal line
       * to  the  i-th   side  of  the  n-gon.   We   use  the  affine
       * interpolation described before.
       */
        for ( unsigned j = 0 ; j < ( max - 1 ) ; j++ ) {
            /*
         * Compute the  ratio: length of  the first "j-1"  segments of
         * the i-th  polygonal line divided by its  total length. This
         * ratio is always a number between  0 and 1 (so, we perform a
         * convex combination).
         */
            double t = df[ i ][ j ] / sm[ i ] ;

            /**
         * Compute the  coordinates of  the point on  the side  of the
         * n-gon. 
         */
            double uu = ( u[ i ] * ( 1 - t ) ) + ( u[ i + 1 ] * t ) ;
            double vv = ( v[ i ] * ( 1 - t ) ) + ( v[ i + 1 ] * t ) ;

            patch()->get_vertex_attributes( bv[ h ] ).set_ucoord( uu ) ;
            patch()->get_vertex_attributes( bv[ h ] ).set_vcoord( vv ) ;

            ++h ;
        }
    }

    return ;
}


/**
   * \fn void FloaterPar::comp_mv_coord( unsigned vi , const std::vector< unsigned >& star , const std::vector< unsigned >& vtab ,  double*& A , double& U , double& V )
   *
   * \brief Computes the mean value coordinates of a given vertex with
   * respect to the vertices in its star. This is an implementation of
   * the method in
   *
   * M. S.  Floater, Mean  value coordinates, Computer Aided Geometric
   * Design 20 (2003), 19-27.
   *
   * \param vi The index of an inner patch vertex.
   * \param star  A list with  the vertices in  the star of  the given
   * inner vertex.
   * \param vtab A lookup table used to speed up access to vertex ID.
   * \param  A An array to store  the mean  value coordinates  of the
   * inner vertices of the star.
   * \param U  A variable to  store the convex  sum of the  mean value
   * coordinates times  the U coordinate  of the boundary  vertices of
   * the star.
   * \param V  A variable to  store the convex  sum of the  mean value
   * coordinates times  the U coordinate  of the boundary  vertices of
   * the star.
   */
void
FloaterPar::comp_mv_coord(
    unsigned vi ,
    const std::vector< unsigned >& star ,
    const std::vector< unsigned >& vtab ,
    //boost::numeric::ublas::matrix<double>& A , int Ai,
    std::vector< std::pair< std::pair<int,int>, double > >& A,
    int Ai,
    double& U ,
    double& V
    )
{
    /*
     * Make sure  there are  at least  3 vertices in  the link  of the
     * given inner vertex.
     */
    assert( star.size() >= 3 ) ;

    /*
     * Create an array to store the vertex weights.
     */
    std::vector< double > w( star.size() ) ;

    /*
     * Create an array to store the coordinates of the vertices in the
     * link  of the  given  inner vertex.  This  is just  to speed  up
     * computations.
     */
    std::vector< double > link_x( star.size() ) ;
    std::vector< double > link_y( star.size() ) ;
    std::vector< double > link_z( star.size() ) ;

    /*
     * Create an array to store the length of the edges connecting the
     * given inner vertex to the vertices in its link.
     */
    std::vector< double > length( star.size() ) ;

    /*
     * Get the coordinates of the given inner vertex.
     */
    double x0 = 0 ;
    double y0 = 0 ;
    double z0 = 0 ;
    patch()->get_vertex_coords( vi , x0 , y0 , z0 ) ;

    for ( unsigned i = 0 ; i < star.size() ; i++ ) {
        unsigned e = patch()->get_next( star[ i ] ) ;
        unsigned v = patch()->get_origin_vertex( e ) ;

        double x1 = 0 ;
        double y1 = 0 ;
        double z1 = 0 ;
        patch()->get_vertex_coords( v , x1 , y1 , z1 ) ;


        length[ i ] = common::Geometric::distance_3d(
                    x0 ,
                    y0 ,
                    z0 ,
                    x1 ,
                    y1 ,
                    z1
                    ) ;

        w[ i ] = 1 / length[ i ] ;

        link_x[ i ] = ( x1 - x0 ) / length[ i ] ;
        link_y[ i ] = ( y1 - y0 ) / length[ i ] ;
        link_z[ i ] = ( z1 - z0 ) / length[ i ] ;
    }

    std::vector< double > a( star.size() ) ;

    for ( unsigned i = 0 ; i < star.size() ; i++ ) {
        double temp =
                ( link_x[ i ] * link_x[ ( i + 1 ) % star.size() ] ) +
                ( link_y[ i ] * link_y[ ( i + 1 ) % star.size() ] ) +
                ( link_z[ i ] * link_z[ ( i + 1 ) % star.size() ] ) ;

        a[ i ] = acos( temp ) ;
    }

    double sum = 0 ;
    for ( unsigned i = 0 ; i < star.size() ; i++ ) {
        sum += a[ i ] ;
    }

    for ( unsigned i = 0 ; i < star.size() ; i++ ) {
        a[ i ] = ( _MYPI * a[ i ] ) / sum ;
    }

    w[ 0 ] *= ( tan( a[ 0 ] ) + tan( a[ star.size() - 1 ] ) ) ;

    assert( w[ 0 ] >= 0 ) ;

    sum = w[ 0 ] ;

    for ( unsigned i = 1 ; i < star.size() ; i++ ) {
        w[ i ] *= ( tan( a[ i - 1 ] ) + tan( a[ i ] ) ) ;

        assert( w[ i ] >= 0 ) ;

        sum  += w[ i ] ;
    }

    for ( unsigned i = 0 ; i < star.size() ; i++ ) {
        unsigned i1 = patch()->get_origin_vertex(
                    patch()->get_next( star[ i ] )
                    ) ;
        unsigned i2 = vtab[ i1 ] ;

        if ( patch()->get_vertex_attributes( i1 ).get_tval() == -1 ) {
            A.push_back( std::pair< std::pair<int,int>, double >(std::pair<int,int>(i2,Ai),-w[ i ] / sum) ) ;
        }
        else {
            double t = w[ i ] / sum ;
            U += ( t * patch()->get_vertex_attributes( i1 ).get_ucoord() ) ;
            V += ( t * patch()->get_vertex_attributes( i1 ).get_vcoord() ) ;
        }
    }

    return ;
}


/**
   * \fn bool FloaterPar::is_inner_vertex( unsigned i ) const
   *
   * \brief Decides whether a given patch vertex is an inner vertex.
   *
   * \param i The index of a patch vertex.
   *
   * \return The logic value true if the the vertex is an inner vertex
   * and the logic value false otherwise.
   */
bool
FloaterPar::is_inner_vertex( unsigned i ) const
{
    /* Create a list to store the edges of the i-th vertex star. */
    std::vector< unsigned > star ;

    /* Get the edges of the i-th vertex star. */
    patch()->get_vertex_star( i , star ) ;

    /*
     * If any edge  of the vertex star is a  boundary edge, the vertex
     * is  not an  inner vertex.  Otherwise,  the vertex  is an  inner
     * vertex.
     */
    for ( unsigned k = 0 ; k < star.size() ; k++ ) {
        if ( patch()->is_boundary_edge( star[ k ] ) ) {
            return false ;
        }
    }

    return true ;
}


/**
   * \fn void FloaterPar::compute_bounding_box( std::vector< BoundingBox >& bbset )
   *
   * \brief Computes a bounding box for each parametrization triangle.
   *
   * \param bbset An array of bounding boxes (one for each triangle).
   */
void
FloaterPar::compute_bounding_box( std::vector< BoundingBox >& bbset )
{
    /**
     * Loop   over  the   set   of  triangles   of  the   parametrized
     * triangulation. For  each triangle,  computes a bounding  box to
     * it.
     */

#ifdef DEBUGMODE
    assert( bbset.size() == get_number_of_faces() ) ;
#endif

    for( unsigned i = 0 ; i < get_number_of_faces() ; i++ ) {
        /*
       * Get the ID of the vertices of the i-th face.
       */
        unsigned v1 = patch()->get_origin_vertex( 3 * i     ) ;
        unsigned v2 = patch()->get_origin_vertex( 3 * i + 1 ) ;
        unsigned v3 = patch()->get_origin_vertex( 3 * i + 2 ) ;

        double x1 = patch()->get_vertex_attributes( v1 ).get_ucoord() ;
        double y1 = patch()->get_vertex_attributes( v1 ).get_vcoord() ;

        double x2 = patch()->get_vertex_attributes( v2 ).get_ucoord() ;
        double y2 = patch()->get_vertex_attributes( v2 ).get_vcoord() ;

        double x3 = patch()->get_vertex_attributes( v3 ).get_ucoord() ;
        double y3 = patch()->get_vertex_attributes( v3 ).get_vcoord() ;

        bbset[ i ]._xmin = x1 ;
        bbset[ i ]._xmax = x1 ;
        bbset[ i ]._ymin = y1 ;
        bbset[ i ]._ymax = y1 ;

        if ( bbset[ i ]._xmin > x2 ) bbset[ i ]._xmin = x2 ;
        if ( bbset[ i ]._xmin > x3 ) bbset[ i ]._xmin = x3 ;
        if ( bbset[ i ]._xmax < x2 ) bbset[ i ]._xmax = x2 ;
        if ( bbset[ i ]._xmax < x3 ) bbset[ i ]._xmax = x3 ;

        if ( bbset[ i ]._ymin > y2 ) bbset[ i ]._ymin = y2 ;
        if ( bbset[ i ]._ymin > y3 ) bbset[ i ]._ymin = y3 ;
        if ( bbset[ i ]._ymax < y2 ) bbset[ i ]._ymax = y2 ;
        if ( bbset[ i ]._ymax < y3 ) bbset[ i ]._ymax = y3 ;
    }

    return ;
}


/**
   * \fn void FloaterPar::compute_image_points( unsigned nv , double* param_pts , const std::vector< BoundingBox >& bbset , std::vector< double >& pset ) 
   *
   * \brief  Compute the  image points  of  a given  set of  parameter
   * points.
   *
   * \param nv The number of parameter points.
   * \param param_pts Coordinates of the parameter points.
   * \param bbset An array of bounding boxes (one for each triangle).
   * \param pset A reference to an array of image points.
   */
void
FloaterPar::compute_image_points(
    unsigned nv ,
    double* param_pts ,
    const std::vector< BoundingBox >& bbset ,
    std::vector< double >& pset
    )
{
    /*
     * For each parameter point, find out the parametrization triangle
     * that contains it. If such  a triangle exists, use a barycentric
     * mapping to compute the 3D coordinates of the image point on the
     * patch.
     */
    for ( unsigned i = 0 ; i < nv ; i++ ) {
        /**
       * Get the coordinates of the i-th parameter point.
       */
        double x = param_pts[ 2 * i     ] ;
        double y = param_pts[ 2 * i + 1 ] ;

        /**
       * Find a triangle containing the point, if any.
       */
        const unsigned nb = bbset.size() ;

        bool found = false ;
        unsigned j = 0 ;

        while ( !found && ( j < nb ) ) {
            /*
         * Perform the point-in-bounding box test.
         */
            if ( is_point_in_bounding_box( x , y , bbset[ j ] ) ) {
                /*
           * Get the ID of the vertices of the j-th triangle.
           */
                unsigned v1 = patch()->get_origin_vertex( 3 * j     ) ;
                unsigned v2 = patch()->get_origin_vertex( 3 * j + 1 ) ;
                unsigned v3 = patch()->get_origin_vertex( 3 * j + 2 ) ;

                /*
           * Get the parameter coordinates of the triangle vertices.
           */
                double x1 =
                        patch()->get_vertex_attributes( v1 ).get_ucoord() ;
                double y1 =
                        patch()->get_vertex_attributes( v1 ).get_vcoord() ;

                double x2 =
                        patch()->get_vertex_attributes( v2 ).get_ucoord() ;
                double y2 =
                        patch()->get_vertex_attributes( v2 ).get_vcoord() ;

                double x3 =
                        patch()->get_vertex_attributes( v3 ).get_ucoord() ;
                double y3 =
                        patch()->get_vertex_attributes( v3 ).get_vcoord() ;

                /*
           * Compute the  barycentric coordinates  of ( x  , y  ) with
           * respect to the triangle with vertices  ( x1 , y1 ) , ( x2
           * , y2 ) and ( x3 , y3 ).
           */
                double u , v , w ;
                compute_barycentric_coordinates(
                            x ,
                            y ,
                            x1 ,
                            y1 ,
                            x2 ,
                            y2 ,
                            x3 ,
                            y3 ,
                            u ,
                            v ,
                            w
                            ) ;

                /*
           * The point ( x , y ) is in the interior of the triangle or
           * in  the   interior  of  its  edges   iff  no  barycentric
           * coordinate is equal to 1.
           */
                if (
                        ( u >= 0 ) && ( u <= 1 ) && ( v >= 0 ) &&
                        ( v <= 1 ) && ( w >= 0 ) && ( w <= 1 )
                        )
                {
                    /*
             * Compute the image of the point ( x , y ) on the patch.
             */
                    double z1 , z2 , z3 ;
                    patch()->get_vertex_coords( v1 , x1 , y1 , z1 ) ;
                    patch()->get_vertex_coords( v2 , x2 , y2 , z2 ) ;
                    patch()->get_vertex_coords( v3 , x3 , y3 , z3 ) ;

                    /*
             * Compute the Cartesian coordinates of the point.
             */
                    x1 = ( u * x1 ) + ( v * x2 ) + ( w * x3 ) ;
                    y1 = ( u * y1 ) + ( v * y2 ) + ( w * y3 ) ;
                    z1 = ( u * z1 ) + ( v * z2 ) + ( w * z3 ) ;

                    pset.push_back( x1 ) ;
                    pset.push_back( y1 ) ;
                    pset.push_back( z1 ) ;

                    /*
             * Get out of the inner for loop.
             */
                    found = true ;
                }
                else {
                    /*
             * Try another bounding box.
             */

                    ++j ;
                }
            }
            else {
                /*
           * Try another bounding box.
           */

                ++j ;
            }
        }

        /*
       * Make sure all parameter points are inside the parametrization
       * domain.
       */
        assert( found ) ;
    }

    return ;
}


/**
   * \fn  void FloaterPar::compute_barycentric_coordinates( double xp , double yp , double x0 , double y0 , double x1 , double y1 , double x2 , double y2 , double& u , double& v , double& w ) const
   * 
   * \brief Computes the barycentric  coordinates of a given point (in
   * Cartesian  coordinates)   with  respect  to   a  given  reference
   * triangle.
   *
   * \param xp First Cartesian coordinate of the point.
   * \param yp Second Cartesian coordinate of the point.
   * \param x0 First  Cartesian coordinate of the first  vertex of the
   * reference triangle.
   * \param y0 Second Cartesian coordinate  of the first vertex of the
   * reference triangle.
   * \param x1 First Cartesian coordinate  of the second vertex of the
   * reference triangle.
   * \param y1 Second Cartesian coordinate of the second vertex of the
   * reference triangle.
   * \param x2 First  Cartesian coordinate of the third  vertex of the
   * reference triangle.
   * \param y2 Second Cartesian coordinate  of the third vertex of the
   * reference triangle.
   * \param u First barycentric coordinate of the point.
   * \param v Second barycentric coordinate of the point.
   * \param w Third barycentric coordinate of the point.
   *
   */
void
FloaterPar::compute_barycentric_coordinates(
    double xp ,
    double yp ,
    double x0 ,
    double y0 ,
    double x1 ,
    double y1 ,
    double x2 ,
    double y2 ,
    double& u ,
    double& v ,
    double& w
    )
const
{
    /*
     * Compute the determinant.
     */
    double dd = ( x1 * y0 ) - ( x2 * y0 ) - ( x0 * y1 )
            + ( x2 * y1 ) + ( x0 * y2 ) - ( x1 * y2 ) ;

    /*
     * The determinant cannot be zero.
     */
    assert( fabs( dd ) > 1e-16 ) ;

    /*
     * Compute the barycentric coordinates. 
     */
    u = ( x2 * y1 ) - ( xp * y1 ) - ( x1 * y2 )
            + ( xp * y2 ) + ( x1 * yp ) - ( x2 * yp ) ;

    u /= dd ;

    v = ( xp * y0 ) - ( x2 * y0 ) + ( x0 * y2 )
            - ( xp * y2 ) - ( x0 * yp ) + ( x2 * yp ) ;

    v /= dd ;

    if ( fabs( u ) < 1e-10 ) {
        u = 0 ;
    }
    else if ( fabs( 1 - u ) < 1e-10 ) {
        u = 1 ;
    }

    if ( fabs( v ) < 1e-10 ) {
        v = 0 ;
    }
    else if ( fabs( 1 - v ) < 1e-10 ) {
        v = 1 ;
    }

    w = 1 - u - v ;

    if ( fabs( w ) < 1e-10 ) {
        w = 0 ;
    }
    else if ( fabs( 1 - w ) < 1e-10 ) {
        w = 1 ;
    }

    return ;
}

}

/** @} */ //end of group class.
