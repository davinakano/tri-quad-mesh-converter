/** 
 * \file geometric.hpp
 *  
 * \brief Definition of a class with a set of geometric methods.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * mfsiqueira at gmail (dot) com
 *
 * \version 1.0
 * \date January 2010
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef GEOMETRIC_HPP
#define	GEOMETRIC_HPP

#include <cmath>               // sqrt, acos


/** 
 * \defgroup COMMONNameSpace Namespace common.
 * @{
 */

/**
 * \namespace common
 *
 * \brief   The   namespace  common   contains   the  definition   and
 * implementation  of several  general classes,  which  contains basic
 * algorithms and/or  basic data structures  that can be used  by more
 * specialized classes.
 */


namespace common {

  /**
   * \class Geometric
   *
   * \brief This class contains a  set of geometric procedures used by
   * the  algorithm  for  computing  discrete  geodesics  on  triangle
   * surface meshes.
   */
  class Geometric {
  public:

    /**
     * \typedef INTERSECTION_TYPE
     *
     * \brief Definition of a type name for all possible return values
     * of  the algorithm that  computes line  segment to  line segment
     * intersection.
     */
    typedef enum {
      NOINTERSECTION = 0 ,
      POINT = 1 ,
      INTERVAL = 2 ,
      PARALLEL_NOINTERSECTION = 3
    }
    INTERSECTION_TYPE ;
    
    // ---------------------------------------------------------------
    //
    // Public methods.
    //
    // ---------------------------------------------------------------

    /**
     * \fn Geometric()
     *
     * \brief Creates an  instance of this class.
     */    
   Geometric()
   {}


    /**
     * \fn ~Geometric()
     *
     * \brief Destroys an instance of this class.
     */    
    ~Geometric() 
    {}


    /**
     * \fn static INTERSECTION_TYPE compute_segment_intersection( double x1 , double y1 , double x2 , double y2 , double x3 , double y3 , double x4 , double y4 , double& sa , double& sb )
     *
     * \brief  Computes the  intersection point  or interval  (if any)
     * between  two given  line segments  in 2D.   This  procedure was
     * taken  from the book  Geomtric Tools  for Computer  Graphics by
     * Philip  J.   Schneider and  David  H. Eberly,  Morgan-Kaufmann,
     * 2003.
     *
     *
     * \param x1 The  X coordinate of the first  endpoint of the first
     * segment.
     * \param y1 The  Y coordinate of the first  endpoint of the first
     * segment.
     * \param x2 The X coordinate  of the second endpoint of the first
     * segment.
     * \param y2 The Y coordinate  of the second endpoint of the first
     * segment.
     * \param x3 The X coordinate  of the first endpoint of the second
     * segment.
     * \param y3 The Y coordinate  of the first endpoint of the second
     * segment.
     * \param x4 The X coordinate of the second endpoint of the second
     * segment.
     * \param y4 The Y coordinate of the second endpoint of the second
     * segment.
     * \param  sa A  reference to  the  parameter value  of the  first
     * segment  that defines  the intersection  point between  the two
     * segments. If there is  no intersection, the value is undefined.
     * If the  intersection is an interval, then  the reference refers
     * to the  first segment parameter value that  defines the closest
     * intersection point between the two segments with respect to the
     * first end point of the first segment.
     * \param  sb A reference  to the  first segment  parameter value,
     * which defines the intersection  point between the two segments.
     * If  there is  no intersection,  the  value referred  to by  the
     * reference is  undefined.  If  the intersection is  an interval,
     * then the reference refers  to the first segment parameter value
     * that defines  the farthest  intersection point between  the two
     * segments  with respect  to the  first  end point  of the  first
     * segment.
     */
    static INTERSECTION_TYPE compute_segment_intersection(
							  double x1 ,
							  double y1 ,
							  double x2 ,
							  double y2 ,
							  double x3 ,
							  double y3 ,
							  double x4 , 
							  double y4 , 
							  double& sa ,
							  double& sb 
							 )
    {
      /*
       * Compute the direction vector of the first segment.
       */
      double v1x = x2 - x1 ;
      double v1y = y2 - y1 ;
      
      /*
       * Compute the direction vector of the second segment.
       */
      double v2x = x4 - x3 ;
      double v2y = y4 - y3 ;
      
      /*
       * Compute  one coordinate of  the cross  product vector  of the
       * direction  vectors. They  can  be used  to  determine if  the
       * segment  supporting  lines  are  parallel (and  the  same  or
       * distinct), or concurrent.
       */
      double kross = ( v1x * v2y ) - ( v1y * v2x ) ;
      
      /*
       * To  reduce  the effects  of  finite precision  floating-point
       * arithmetic, we perform a test  using the square of kross, and
       * the squares of the length of the direction vectors (they need
       * not be unit vectors).
       */
      double sqr_kross = kross * kross ;
      double sqr_dirv1 = ( v1x * v1x ) + ( v1y * v1y ) ;
      double sqr_dirv2 = ( v2x * v2x ) + ( v2y * v2y ) ;
      
      double ex = x3 - x1 ;
      double ey = y3 - y1 ;
      
      /*
       * Get  the square  of the  tolerance to  consider  two segments
       * parallel.
       */
      const double _TOL_T_SEGP =  1e-16 ;

      double sqr_toler = _TOL_T_SEGP * _TOL_T_SEGP ;
      
      if ( sqr_kross > ( sqr_dirv1 * sqr_dirv2 * sqr_toler ) ) {
	/*
	 * The segment supporting lines are not parallel.
	 */
	
	double s = ( ( ex * v2y ) - ( ey * v2x ) ) / kross ;
	double t = ( ( ex * v1y ) - ( ey * v1x ) ) / kross ;
	
	sa = s ;
	sb = t ;
	
	if ( ( s < 0 ) || ( s > 1 ) ) {
	  return NOINTERSECTION ;
	}
	
	if ( ( t < 0 ) || ( t > 1 ) ) {	
	  return NOINTERSECTION ;
	}
	
	/*
	 * Segments intersect at exactly one point.
	 *
	 * xa = xb = x1 + s * v1x ;
	 * ya = yb = y1 + s * v1y ;
	 */
	return POINT ;
      }
      
      /*
       * The segment supporting lines are parallel
       */
      
      /*
       * Find out whether the segments are on the same supporting line.
       */
      double sqr_vece = ( ex * ex ) + ( ey * ey ) ;
      kross = ( ex * v1y ) - ( ey * v1x ) ;
      sqr_kross = kross * kross ;
      if ( sqr_kross > ( sqr_dirv1 * sqr_vece * sqr_toler ) ) {
	/*
	 * The segments belong to distinct supporting lines.
	 */
	
	return PARALLEL_NOINTERSECTION ;
      }
      
      /*
       * Find out  whether the segments  overlap.  If so,  compute the
       * overlapping interval.
       */
      double s1 = ( ( v1x * ex ) + ( v1y * ey ) ) / sqr_dirv1 ;
      double s2 = s1 + ( ( v1x * v2x ) + ( v1y * v2y ) ) / sqr_dirv1 ;

      double smin = ( s1 <= s2 ) ? s1 : s2 ;
      double smax = ( s1 >= s2 ) ? s1 : s2 ;
      
      if ( ( smin > 1.0 ) || ( smax < 0.0 ) ) {
	/*
	 * Segments do not overlap each other.
	 */
	return PARALLEL_NOINTERSECTION ;
      }
      
      if ( smin < 1.0 ) {
	/*
	 * There is an overlap at one point or an interval.
	 */
	if ( smax > 0.0 ) {
	  /*
	   * The overlap is an interval.
	   */
	  if ( smin > 0.0 ) {
	    // xa = x1 + smin * v1x ;
	    // ya = y1 + smin * v1y ;
	    
	    sa = smin ;  
	  } 
	  else {
	    // xa = x1 ;
	    // ya = y1 ;
	    
	    sa = 0 ;
	  }
	  
	  if ( smax < 1.0 ) {
	    // xb = x1 + smax * v1x ;
	    // yb = y1 + smax * v1y ;
	    
	    sb = smax ;
	    
	  } 
	  else {
	    // xb = x2 ;
	    // yb = y2 ;
	    
	    sb = 1 ;
	  }
	  
	  return INTERVAL ;
	}
	else {
	  /*
	   * The overlap is an endpoint.
	   *
	   * xa = xb = x1 ;
	   * ya = yb = y1 ;
	   */
	  
	  sa = 0 ;
	  sb = 1 ;
	  
	  return POINT ;
	}
      } 
      else {
	/*
	 * If we get here, it is because smin is equal to 1.
	 *
	 * xa = xb = x2 ;
	 * ya = yb = y2 ;
	 */
	
	sa = 1 ;
	sb = 0 ;
	
	return POINT ;
      }
      
      /*
       * This point should never be reached.     
       */
      
      assert( false ) ;
      
      return NOINTERSECTION ;
    }
    
    
    /**
     * \fn static double get_angle( double xa , double ya , double za , double xb , double yb , double zb , double xc , double yc , double zc )
     *
     * \brief Computes  the angle defined at  the vertex A  = (xa, ya,
     * za) of a given triangle ABC, where B = ( xb, yb, zb ) and C = (
     * xc, yc, zc ).
     *
     * \param xa The X coordinate of vertex A.
     * \param ya The Y coordinate of vertex A.
     * \param za The Z coordinate of vertex A.
     * \param xb The X coordinate of vertex B.
     * \param yb The Y coordinate of vertex B.
     * \param zb The Z coordinate of vertex B.
     * \param xc The X coordinate of vertex C.
     * \param yc The Y coordinate of vertex C.
     * \param zc The Z coordinate of vertex C.
     *
     * \return The angle value.
     */
    static double get_angle(
			    double xa ,
			    double ya ,
			    double za ,				     
			    double xb ,
			    double yb ,
			    double zb ,
			    double xc ,
			    double yc ,
			    double zc 
			   )
    {
      /*
       * Compute vector (B - A).
       */
      double abx = xb - xa ;
      double aby = yb - ya ;
      double abz = zb - za ;
      
      /*
       * Compute vector (C - A).
       */
      double acx = xc - xa ;
      double acy = yc - ya ;
      double acz = zc - za ;
      
      /*
       * Compute the length of vector (B - A).
       */
      double lab = sqrt( 
			( abx * abx ) + 
			( aby * aby ) + 
			( abz * abz ) 
		       ) ;
      
      /*
       * Compute the length of vector (C - A).
       */
      double lac = sqrt( 
			( acx * acx ) + 
			( acy * acy ) + 
			( acz * acz ) 
		       ) ;
      
      /*
       * If one  of (B - A)  and (C - A)  is almost a  null vector, we
       * return the angle PI. This return values always makes sense in
       * the context of the computation of geodesics using the present
       * algorithm.
       */
      if ( ( lab < 1e-14 ) && ( lac < 1e-14 ) ) {
	return acos( -1 ) ;
      }
      
      /*
       * Normalize the vectors (B - A) and (C - A).
       */
      abx /= lab ;
      aby /= lab ;
      abz /= lab ;
      
      acx /= lac ;
      acy /= lac ;
      acz /= lac ;
      
      /*
       * Compute the cosine of the angle we wish to compute .
       */
      double cost = ( abx * acx ) + ( aby * acy ) + ( abz * acz ) ;
      
      /*
       * Make sure the cosine is in the interval [-1,1].
       */
      if ( cost > 1 ) {
	cost = 1 ;
      } 
      else if ( cost < -1 ) {
	cost = -1 ;
      }
      
      /*
       * Return the angle defined by (B - A) and (C - A).
       */
      return acos( cost ) ;
    }
    

    /**
     * \fn static void compute_third_vertex_position( double xa , double ya , double xb , double yb , double lac , double lbc , double& xc , double& yc )
     *
     * \brief  Computes the Cartesian  coordinates of  the image  of a
     * triangle  vertex  in  R3  under  an isometric  mapping  of  the
     * triangle  to  the plane.   The  planar  coordinates (under  the
     * isometric mapping)  of two vertices of the  triangle are given,
     * as well as the distances from these two vertices to the (third)
     * vertex we wish to compute the planar coordinates.
     *
     * \param xa The X planar coordinate of the first triangle vertex.
     * \param ya The X planar coordinate of the first triangle vertex.
     * \param  xb  The X  planar  coordinate  of  the second  triangle
     * vertex.
     * \param  yb  The Y  planar  coordinate  of  the second  triangle
     * vertex.
     * \param  lac Euclidean  distance  from the  first  to the  third
     * vertex.
     * \param  lbc Euclidean  distance from  the second  to  the third
     * vertex.
     * \param xc A  reference to the X planar  coordinate of the third
     * vertex.
     * \param yc A  reference to the Y planar  coordinate of the third
     * vertex.
     */
    static void compute_third_vertex_position(
					      double xa ,
					      double ya ,
					      double xb ,
					      double yb ,
					      double lac ,
					      double lbc ,
					      double& xc ,
					      double& yc 
					     ) 
    {
      /*
       * The planar coordinates of the third vertex can be computed by
       * solving the two circle  intersection problem. In other words,
       * the point C = (xc, yc) we  wish to compute can be seen as one
       * of the two intersection points of the circles centered at A =
       * (xa,  ya)   and  B  =   (xb,yb)  with  radii  lac   and  lbc,
       * respectively.  We start by computing the point P1 = ( x1, y1)
       * corresponding  to the projection  of (xc,  yc) onto  the edge
       * connecting A to B.  To do that, we first compute the distance
       * from A to P1, and from P1 to B. Next, we compute the height h
       * (the distance from  P1 to C) of triangle  ABC with respect to
       * edge AB.  Finally, we use triangle similarity  to compute the
       * offsets  of  the  coordinates   of  C  with  respect  to  the
       * coordinates of P1. Once we have these offsets, we have C.
       */

      /*
       * Compute the distance from A to B.
       */
      double lab = sqrt( 
			( xa - xb ) * ( xa - xb ) + 
			( ya - yb ) * ( ya - yb ) 
		       ) ;
    
      /*
       * Compute the unit vector in the direction from A to B.
       */
      double vbax = ( xb - xa ) ;
      double vbay = ( yb - ya ) ;
      double norm = sqrt( ( vbax * vbax ) + ( vbay * vbay ) ) ;
      
      vbax = ( 1 / norm ) * vbax ;
      vbay = ( 1 / norm ) * vbay ;
      
      /*
       * Compute the (signed) distance from B to P1.
       */
      double dbp1 = ( ( lab * lab ) + ( lbc * lbc ) - ( lac * lac ) ) 
	/ ( 2 * lab ) ;
      
      /*
       * Compute the coordinates of P1.
       */
      const double length = lab - dbp1 ;
      double x1 = xa + length * vbax ;
      double y1 = ya + length * vbay ;
      
      /*
       * Compute the distance from C to P1.
       */
      double dcp1 = sqrt( ( lbc * lbc ) - ( dbp1 * dbp1 ) ) ;
      
      /*
       * Compute the offsets  of the coordinates of C  with respect to
       * P1.
       */
      const double r2 = dcp1 / lab ;
      
      double offset_x = r2 * ( yb - ya ) ;
      double offset_y = r2 * ( xb - xa ) ;
      
      /*
       * Compute the coordinates of C.
       */
      xc = x1 - offset_x ;
      yc = y1 + offset_y ;
      
      return ;
    }
    

    /**
     * \fn static double norm_3d( double x , double y , double z )
     *
     * \brief Computes the length of a given 3D vector.
     *
     * \param x The first Cartesian coordinate of the vector.
     * \param y The second Cartesian coordinate of the vector.
     * \param z The third Cartesian coordinates of the vector.
     *
     * \return The length of the given 3D vector.
     */
    static double norm_3d(
			  double x , 
			  double y , 
			  double z
			 )
    {
      return sqrt( ( x * x ) + ( y * y ) + ( z * z ) ) ;
    }


    /**
     * \fn static double distance_3d( double xa , double ya , double za , double xb , double yb , double zb )
     *
     * \brief  Computes   the  Euclidean  distance   between  two  given
     * points.
     *
     * \param xa The first Cartesian coordinate of the first point.
     * \param ya The second Cartesian coordinate of the first point.
     * \param za The third Cartesian coordinates of the first point.
     * \param xb The first Cartesian coordinate of the second point.
     * \param yb The second Cartesian coordinate of the second point.
     * \param zb The third Cartesian coordinates of the second point.
     *
     * \return The Euclidean distance between two given vertices.
     */
    static double distance_3d( 
			      double xa ,
			      double ya ,
			      double za , 
			      double xb , 
			      double yb ,
			      double zb 
			     )
    {
      double dx = xa - xb ;
      double dy = ya - yb ;
      double dz = za - zb ;
      
      return norm_3d( dx , dy , dz ) ;
    } 


    /**
     * \fn static double dot_3d( double xa , double ya , double za , double xb , double yb , double zb )
     *
     * \brief Computes the dot product of two 3D vectors.
     *
     * \param xa The first Cartesian coordinate of the first vector.
     * \param ya The second Cartesian coordinate of the first vector.
     * \param za The third Cartesian coordinates of the first vector.
     * \param xb The first Cartesian coordinate of the second vector.
     * \param yb The second Cartesian coordinate of the second vector.
     * \param zb The third Cartesian coordinates of the second vector.
     *
     * \return The dot product of the two given 3D vectors.
     */
    static double dot_3d(
			 double xa ,
			 double ya ,
			 double za , 
			 double xb , 
			 double yb , 
			 double zb
			)
    {
      return ( xa * xb ) + ( ya * yb ) + ( za * zb ) ;
    }


    /**
     * \fn static void cross_3d( double xa , double ya , double za , double xb , double yb , double zb , double& xc , double& yc , double& zc )
     *
     * \brief Computes the cross product of two given 3D vectors.
     *
     * \param xa The first Cartesian coordinate of the first vector.
     * \param ya The second Cartesian coordinate of the first vector.
     * \param za The third Cartesian coordinates of the first vector.
     * \param xb The first Cartesian coordinate of the second vector.
     * \param yb The second Cartesian coordinate of the second vector.
     * \param zb The third Cartesian coordinates of the second vector.
     * \param xc The first Cartesian coordinate of the resulting vector.
     * \param yc The second Cartesian coordinate of the resulting vector.
     * \param zc The third Cartesian coordinates of the resulting vector.
     */
    static void cross_3d(
			 double xa ,
			 double ya ,
			 double za , 
			 double xb , 
			 double yb , 
			 double zb , 
			 double& xc , 
			 double& yc , 
			 double& zc 
			)
    {
      xc = ( ya * zb ) - ( za * yb ) ; 
      yc = ( xa * zb ) - ( za * xb ) ; 
      zc = ( xa * yb ) - ( ya * xb ) ;
      
      return ;
    }


    /**
     * \fn static double compute_triangle_area( double l1 , double l2 , double l3 )
     *
     * \brief Computes the area of a triangle using the lengths of its
     * sides.
     *
     * \param l1 The length of the first side of the triangle.
     * \param l2 The length of the second side of the triangle.
     * \param l3 The length of the third side of the triangle.
     *
     * \return The area of the given triangle.
     */
    static double compute_triangle_area(
					double l1,
					double l2,
					double l3
				       )
    {
      return 0.25 * sqrt (
			    ( l1 + l2 + l3 ) 
			  * ( l2 + l3 - l1 ) 
			  * ( l3 + l1 - l2 ) 
			  * ( l1 + l2 - l3 ) 
			 ) ;
    }


    /**
     * \fn static void compute_triangle_normal( double xa , double ya , double za , double xb , double yb , double zb , double xc , double yc , double zc , double& nx , double& ny , double& nz )
     *
     * \brief Computes the unit normal  of a triangle. The vertices of
     * the triangle  are given in counterclockwise order,  so that the
     * normal is computed according to the order in which the vertices
     * are given.
     *
     * \param xa The first  Cartesian coordinate of the first triangle
     * vertex.  
     * \param ya The second Cartesian coordinate of the first triangle
     * vertex.
     * \param za The third Cartesian coordinates of the first triangle
     * vertex.
     * \param xb The first Cartesian coordinate of the second triangle
     * vertex.
     * \param  yb  The  second  Cartesian  coordinate  of  the  second
     * triangle vertex.
     * \param  zb  The  third  Cartesian  coordinates  of  the  second
     * triangle vertex.
     * \param xc The first  Cartesian coordinate of the third triangle
     * vertex.
     * \param yc The second Cartesian coordinate of the third triangle
     * vertex.
     * \param zc The third Cartesian coordinates of the third triangle
     * vertex.
     * \param  nx  The  first  Cartesian coordinate  of  the  triangle
     * normal.
     * \param  ny  The second  Cartesian  coordinate  of the  triangle
     * normal.
     * \param  nz  The third  Cartesian  coordinates  of the  triangle
     * normal.
     *
     * \return The unit normal to the given triangle.
     */
    static void compute_triangle_normal(
					double xa ,
					double ya , 
					double za ,
					double xb , 
					double yb , 
					double zb , 
					double xc , 
					double yc , 
					double zc ,
					double& nx , 
					double& ny , 
					double& nz 
				       )
    {
      double a1 = xb - xa ;
      double a2 = yb - ya ;
      double a3 = zb - za ;
      
      double b1 = xc - xa ;
      double b2 = yc - ya ;
      double b3 = zc - za ;
      
      cross_3d( a1 , a2 , a3 , b1 , b2 , b3 , nx , ny , nz ) ;
      
      double norm = norm_3d( nx , ny , nz ) ;
      
      nx /= norm ;
      ny /= norm ;
      nz /= norm ;
      
      return ;
    }


    /**
     * \fn static double circumradius( double a , double b , double c )
     *
     * \brief Computes  the circumradius of a triangle  from the lengths
     * of the sides of the triangle.
     *
     * \param a The length of the first side.
     * \param b The length of the second side
     * \param c The length of the third side.
     *
     * \return The circumradius of a triangle.
     */
    static double circumradius(
			       double a ,
			       double b ,
			       double c 
			      )
    {
      double abc = a * b * c ;
      double den = ( a + b + c ) * ( b + c - a ) 
	         * ( c + a - b ) * ( a + b - c ) ;

      return abc / sqrt( den ) ;
    }


  } ;

}

/** @} */ //end of group class.

#endif	/* GEOMETRIC_HPP */

