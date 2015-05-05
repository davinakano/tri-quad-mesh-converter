#include "types.h"
#include "myAdgGeodesics.h"



bool myApproxGeodesics::inRegion( Face* f, QFace* qf ) const {

#ifdef DEBUGMODE
  assert( f != 0 ) ;
#endif
  std::set<QFace*>::iterator it = f->get_attributes().region.find( qf );

  return it != f->get_attributes().region.end();
}




/**
 * \fn void myApproxGeodesics::get_restricted_adj_vertices( Vertex* v , std::vector< Vertex* >& a ) const
 *
 * \brief Finds all  mesh vertices connected to a  given mesh vertex
 * by an  edge that is  incident to a  face belonging to  the active
 * face list.
 *
 * \param v A pointer to a vertex.
 * \param a A reference to a list of vertices.
 */
void
myApproxGeodesics::get_restricted_adj_vertices(
                                                      Vertex* v ,
                                                      std::vector< Vertex* >& a
                                                     )
  const
{
  Halfedge* h1 = _mesh->get_halfedge( v ) ;
  Halfedge* h2 = h1 ;

  do {
    Face* f1 = _mesh->get_face( h2 ) ;
    Face* f2 = _mesh->get_face( _mesh->get_mate( h2 ) ) ;

    /*
     * Insert the origin vertex of  "he2" in the list of vertices if
     * and only if the edge  connecting it to vertex "v" is incident
     * with an active face.
     */
    if ( _mesh->is_active( f1 ) || _mesh->is_active( f2 ) ) {

        if( (!_region)||(inRegion(f1,_region))||(inRegion(f2,_region) ) )
            a.push_back( _mesh->get_org( _mesh->get_next( h2 ) ) ) ;
    }

    h2 = _mesh->get_next( _mesh->get_mate( h2 ) ) ;
  }
  while ( h2 != h1 ) ;

  return ;
}
