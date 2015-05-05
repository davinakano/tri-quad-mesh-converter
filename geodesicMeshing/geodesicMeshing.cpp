#include "geodesicMeshing/geodesicMeshing.h"
#include "meshref.h"

geodesicMeshing::geodesicMeshing( TMesh* refinMesh ) {

    _refinMesh = refinMesh;

    _adgMesh = new adgFDCel( _refinMesh );

    _geo = new geodesics( _adgMesh, 1e-10, 10000) ;
}

geodesicMeshing::~geodesicMeshing() {
    delete _geo;
    delete _adgMesh;
}

bool geodesicMeshing::getGeodesic( QFace* region, Vertex* vi, Vertex* ve, std::list< Vertex* >& vertices ) {

    geodesics::GP* gvi = new geodesics::GPV( vi );
    geodesics::GP* gve = new geodesics::GPV( ve );

    std::list< geodesics::GP* > pts_first ;
    std::list< geodesics::GP* > pts ;

    _geo->setRegion(region);

    _geo->compute_first_approximation( gvi , gve , pts_first ) ;
    _geo->shorten_geodesic( pts_first , pts ) ;

    _todelete.insert( pts_first.begin(), pts_first.end() );
    _todelete.insert( pts.begin(), pts.end() );

    remeshingGeodesics( pts, vertices );
    //remeshingGeodesics( pts_first, vertices );

    deleteGPs(); // memory leaking?

    return true;
}


bool geodesicMeshing::deleteGPs() {

    for( std::set<geodesics::GP*>::iterator it = _todelete.begin(); it != _todelete.end(); ++it ) {
        geodesics::GP* pt = *it;
        delete pt;
    }

    _todelete.clear();

    return true;
}

bool geodesicMeshing::remeshingGeodesics( std::list< geodesics::GP* >& pts, std::list< Vertex* >& vertices ) {

    std::list< geodesics::GP* >::iterator it1 = pts.begin();
    std::list< geodesics::GP* >::iterator it2 = it1;

    ASSERT( pts.front()->on_vertex() );
    ASSERT( pts.back()->on_vertex() );

    Halfedge* he_left = 0;
    Halfedge* he_right = 0;

    for( ++it2; it2 != pts.end(); ++it1, ++it2) {

        geodesics::GP* gp_prev = *it1;
        geodesics::GP* gp_next = *it2;

        if( gp_prev->on_vertex() && gp_next->on_edge() ) {

            geodesics::GPV* gpv_prev = dynamic_cast<geodesics::GPV*>(gp_prev);
            geodesics::GPE* gpe_next = dynamic_cast<geodesics::GPE*>(gp_next);

            split2( gpv_prev, gpe_next, he_left, he_right );

            vertices.push_back( he_left->get_origin() );

        } else if( gp_prev->on_edge() && gp_next->on_edge() ) {

            geodesics::GPE* gpe_prev = dynamic_cast<geodesics::GPE*>(gp_prev);
            geodesics::GPE* gpe_next = dynamic_cast<geodesics::GPE*>(gp_next);

            split3( gpe_prev, gpe_next, he_left, he_right );

            vertices.push_back( he_left->get_origin() );

        } else if( gp_prev->on_edge() && gp_next->on_vertex() ) {

            geodesics::GPE* gpe_prev = dynamic_cast<geodesics::GPE*>(gp_prev);
            geodesics::GPV* gpv_next = dynamic_cast<geodesics::GPV*>(gp_next);

            split2( gpe_prev, gpv_next, he_left, he_right );

            vertices.push_back( gpv_next->get_vertex() );

        } else if( gp_prev->on_vertex() && gp_next->on_vertex() ) {

            geodesics::GPV* gpv_next = dynamic_cast<geodesics::GPV*>(gp_next);

            he_left = he_right = 0;

            vertices.push_back( gpv_next->get_vertex() );

        } else {
            ASSERT( false );
        }
    }

    ASSERT( (*it1)->on_vertex() );

    vertices.pop_back();

    return true;
}

bool geodesicMeshing::split2( geodesics::GPV* gpv_prev, geodesics::GPE* gpe_next, Halfedge*& he_left, Halfedge*& he_right) {

    ASSERT( gpv_prev );
    ASSERT( gpe_next );
    ASSERT( he_left == 0 );
    ASSERT( he_right == 0 );


    // defining reference halfedge
    Halfedge* h = 0;

    Halfedge* htemp[2] = { gpe_next->get_edge()->get_first_halfedge(), gpe_next->get_edge()->get_second_halfedge() };
    if( htemp[0]->get_prev()->get_origin() == gpv_prev->get_vertex() )
        h = htemp[0]->get_prev();
    else if( htemp[1]->get_prev()->get_origin() == gpv_prev->get_vertex() )
        h = htemp[1]->get_prev();
    else {
        ASSERT( false );
    }

    // initialization
    Halfedge* h_next = h->get_next();
    Halfedge* h_prev = h->get_prev();

    //Halfedge* hv = h->get_mate();
    Halfedge* hv_next = h_next->get_mate();
    Halfedge* hv_prev = h_prev->get_mate();

    //Edge* e = h->get_edge();
    //Edge* e_next = h_next->get_edge();
    Edge* e_prev = h_prev->get_edge();

    Vertex* v = h->get_origin();
    //Vertex* v_next = h_next->get_origin();
    Vertex* v_prev = h_prev->get_origin();

    // creating new structures
    Halfedge* new_h = new Halfedge(0,0,0,0,0);
    Halfedge* new_h_next = new Halfedge(0,0,0,0,0);
    Halfedge* new_h_prev = new Halfedge(0,0,0,0,0);

    double nx, ny, nz;
    gpe_next->get_coords( _adgMesh, nx, ny, nz );
    Vertex* new_vertex = new Vertex(nx,ny,nz, 0);

    Edge* new1_edge = new Edge( 0, 0 );
    Edge* new2_edge = new Edge( 0, 0 );

    Face* new_cell = new Face( 0 );

    // new vertices
    new_vertex->set_halfedge( new_h_next );

    // old vertices
    v_prev->set_halfedge( new_h_prev );

    // new halfedges
    new_h->set_next( new_h_next );
    new_h_next->set_next( new_h_prev );
    new_h_prev->set_next( new_h );

    new_h->set_prev( new_h_prev );
    new_h_next->set_prev( new_h );
    new_h_prev->set_prev( new_h_next );

    new_h->set_origin( v );
    new_h_next->set_origin( new_vertex );
    new_h_prev->set_origin( v_prev );

    new_h->set_edge( new1_edge );
    new_h_next->set_edge( new2_edge );
    new_h_prev->set_edge( e_prev );

    new_h->set_face( new_cell );
    new_h_next->set_face( new_cell );
    new_h_prev->set_face( new_cell );

    // old halfedges
    h_prev->set_origin( new_vertex );
    h_prev->set_edge( new1_edge );

    // new edges
    new1_edge->set_first_halfedge( h_prev );
    new1_edge->set_second_halfedge( new_h );

    new2_edge->set_first_halfedge( new_h_next );
    new2_edge->set_second_halfedge( hv_next );

    // old edges
    e_prev->set_first_halfedge( hv_prev );
    e_prev->set_second_halfedge( new_h_prev );

    // new cells
    new_cell->set_halfedge( new_h );

    // adding new structures
    _refinMesh->add_vertex( new_vertex );
    _refinMesh->add_edge( new1_edge );
    _refinMesh->add_edge( new2_edge );
    _refinMesh->add_face( new_cell );

    // returning halfedges
    he_left = new_h_next;
    he_right = h_next;

    return true;
}

bool geodesicMeshing::split2( geodesics::GPE* gpe_prev, geodesics::GPV* gpv_next, Halfedge*& he_left, Halfedge*& he_right) {

    ASSERT( gpe_prev );
    ASSERT( gpv_next );
    ASSERT( he_left );
    ASSERT( he_right );

    // defining reference halfedge
    Halfedge* h = 0;

    Halfedge* htemp[2] = { gpe_prev->get_edge()->get_first_halfedge(), gpe_prev->get_edge()->get_second_halfedge() };
    if( htemp[0]->get_prev()->get_origin() == gpv_next->get_vertex() )
        h = htemp[0];
    else if( htemp[1]->get_prev()->get_origin() == gpv_next->get_vertex() )
        h = htemp[1];
    else {
        ASSERT( false );
    }

    // initialization
    Halfedge* h_next = h->get_next();
    Halfedge* h_prev = h->get_prev();

    Halfedge* hv_next = h_next->get_mate();
    //Halfedge* hv_prev = h_prev->get_mate();

    Edge* e1 = he_left->get_edge();
    Edge* e2 = he_right->get_edge();
    Edge* e_next = h_next->get_edge();
    //Edge* e_prev = h_prev->get_edge();

    //Vertex* v = h->get_origin();
    Vertex* v_next = h_next->get_origin();
    Vertex* v_prev = h_prev->get_origin();
    Vertex* mid_v = he_left->get_origin();

    // creating new structures
    Halfedge* new_h = new Halfedge(0,0,0,0,0);
    Halfedge* new_h_next = new Halfedge(0,0,0,0,0);
    Halfedge* new_h_prev = new Halfedge(0,0,0,0,0);

    Edge* new_edge = new Edge( 0, 0 );

    Face* new_cell = new Face( 0 );

    // old vertices
    v_next->set_halfedge( new_h_next );
    mid_v->set_halfedge( new_h );

    // new halfedges
    new_h->set_next( new_h_next );
    new_h_next->set_next( new_h_prev );
    new_h_prev->set_next( new_h );

    new_h->set_prev( new_h_prev );
    new_h_next->set_prev( new_h );
    new_h_prev->set_prev( new_h_next );

    new_h->set_origin( mid_v );
    new_h_next->set_origin( v_next );
    new_h_prev->set_origin( v_prev );

    new_h->set_edge( e2 );
    new_h_next->set_edge( e_next );
    new_h_prev->set_edge( new_edge );

    new_h->set_face( new_cell );
    new_h_next->set_face( new_cell );
    new_h_prev->set_face( new_cell );

    // old halfedges
    h->set_edge( e1 );
    h_next->set_edge( new_edge );
    h_next->set_origin( mid_v );

    // new edges
    new_edge->set_first_halfedge( h_next );
    new_edge->set_second_halfedge( new_h_prev );

    // old edges
    e1->set_first_halfedge( h  );
    e1->set_second_halfedge( he_left );

    e2->set_first_halfedge( new_h  );
    e2->set_second_halfedge( he_right );

    e_next->set_first_halfedge( hv_next  );
    e_next->set_second_halfedge( new_h_next );

    // new cells
    new_cell->set_halfedge( new_h );

    // adding new structures
    _refinMesh->add_edge( new_edge );
    _refinMesh->add_face( new_cell );

    // returning halfedges
    he_left = 0;
    he_right = 0;

    return true;
}

bool geodesicMeshing::split3( geodesics::GPE* gpe_prev, geodesics::GPE* gpe_next, Halfedge*& he_left, Halfedge*& he_right) {

    ASSERT( gpe_prev );
    ASSERT( gpe_next );
    ASSERT( he_left );
    ASSERT( he_right );

    // defining reference halfedge
    Halfedge* h = 0;

    int type = -1;

    Halfedge* htemp[2] = { gpe_prev->get_edge()->get_first_halfedge(), gpe_prev->get_edge()->get_second_halfedge() };
    if( htemp[0]->get_prev()->get_edge() == gpe_next->get_edge() ) {
        h = htemp[0];
        type = 0;
    }
    else if( htemp[1]->get_prev()->get_edge() == gpe_next->get_edge() ) {
        h = htemp[1];
        type = 0;
    }
    else if( htemp[0]->get_next()->get_edge() == gpe_next->get_edge() ) {
        h = htemp[0];
        type = 1;
    }
    else if( htemp[1]->get_next()->get_edge() == gpe_next->get_edge() ) {
        h = htemp[1];
        type = 1;
    }
    else {
        ASSERT( false );
    }

    // initialization
    Halfedge* h_next = h->get_next();
    Halfedge* h_prev = h->get_prev();

    Halfedge* hv_next = h_next->get_mate();
    Halfedge* hv_prev = h_prev->get_mate();

    Edge* e1 = he_left->get_edge();
    Edge* e2 = he_right->get_edge();
    Edge* e_next = h_next->get_edge();
    Edge* e_prev = h_prev->get_edge();

    Vertex* v = h->get_origin();
    Vertex* v_next = h_next->get_origin();
    Vertex* v_prev = h_prev->get_origin();
    Vertex* mid_v = he_left->get_origin();

    // creating new structures
    Halfedge* new1_h = new Halfedge(0,0,0,0,0);
    Halfedge* new1_h_next = new Halfedge(0,0,0,0,0);
    Halfedge* new1_h_prev = new Halfedge(0,0,0,0,0);

    Halfedge* new2_h = new Halfedge(0,0,0,0,0);
    Halfedge* new2_h_next = new Halfedge(0,0,0,0,0);
    Halfedge* new2_h_prev = new Halfedge(0,0,0,0,0);

    double nx, ny, nz;
    gpe_next->get_coords( _adgMesh, nx, ny, nz );
    Vertex* new_vertex = new Vertex(nx,ny,nz, 0);

    Edge* new1_edge = new Edge( 0, 0 );
    Edge* new2_edge = new Edge( 0, 0 );
    Edge* new3_edge = new Edge( 0, 0 );

    Face* new1_cell = new Face( 0 );
    Face* new2_cell = new Face( 0 );

    if ( type == 0 ) {

        // new vertices
        new_vertex->set_halfedge( h_prev );

        // old vertices
        v_next->set_halfedge( new2_h_next );
        v_prev->set_halfedge( new2_h_prev );
        mid_v->set_halfedge(new1_h );

        // new halfedges
        new1_h->set_next( new1_h_next );
        new1_h_next->set_next( new1_h_prev );
        new1_h_prev->set_next( new1_h );
        new2_h->set_next( new2_h_next );
        new2_h_next->set_next( new2_h_prev );
        new2_h_prev->set_next( new2_h );

        new1_h->set_prev( new1_h_prev );
        new1_h_next->set_prev( new1_h );
        new1_h_prev->set_prev( new1_h_next );
        new2_h->set_prev( new2_h_prev );
        new2_h_next->set_prev( new2_h );
        new2_h_prev->set_prev( new2_h_next );

        new1_h->set_origin( mid_v );
        new1_h_next->set_origin( v_next );
        new1_h_prev->set_origin( new_vertex );
        new2_h->set_origin( new_vertex );
        new2_h_next->set_origin( v_next );
        new2_h_prev->set_origin( v_prev );

        new1_h->set_edge( e2 );
        new1_h_next->set_edge( new1_edge );
        new1_h_prev->set_edge( new3_edge );
        new2_h->set_edge( new1_edge );
        new2_h_next->set_edge( e_next );
        new2_h_prev->set_edge( e_prev );

        new1_h->set_face( new1_cell );
        new1_h_next->set_face( new1_cell );
        new1_h_prev->set_face( new1_cell );
        new2_h->set_face( new2_cell );
        new2_h_next->set_face( new2_cell );
        new2_h_prev->set_face( new2_cell );

        // old halfedges
        h_next->set_origin( mid_v );
        h_prev->set_origin( new_vertex );

        h_next->set_edge( new3_edge );
        h_prev->set_edge( new2_edge );
        h->set_edge( e1 );

        // new edges
        new1_edge->set_first_halfedge( new1_h_next );
        new1_edge->set_second_halfedge( new2_h );

        new2_edge->set_first_halfedge( h_prev );
        new2_edge->set_second_halfedge( hv_prev );

        new3_edge->set_first_halfedge( new1_h_prev );
        new3_edge->set_second_halfedge( h_next );

        // old edges
        e1->set_first_halfedge( he_left );
        e1->set_second_halfedge( h );

        e2->set_first_halfedge( he_right );
        e2->set_second_halfedge( new1_h );

        e_next->set_first_halfedge( hv_next );
        e_next->set_second_halfedge( new2_h_next );

        e_prev->set_first_halfedge( hv_prev );
        e_prev->set_second_halfedge( new2_h_prev );

    } else {

        // new vertices
        new_vertex->set_halfedge( h_prev );

        // old vertices
        v->set_halfedge( new1_h );
        v_prev->set_halfedge( new2_h_prev );
        mid_v->set_halfedge( h );

        // new halfedges
        new1_h->set_next( new1_h_next );
        new1_h_next->set_next( new1_h_prev );
        new1_h_prev->set_next( new1_h );
        new2_h->set_next( new2_h_next );
        new2_h_next->set_next( new2_h_prev );
        new2_h_prev->set_next( new2_h );

        new1_h->set_prev( new1_h_prev );
        new1_h_next->set_prev( new1_h );
        new1_h_prev->set_prev( new1_h_next );
        new2_h->set_prev( new2_h_prev );
        new2_h_next->set_prev( new2_h );
        new2_h_prev->set_prev( new2_h_next );

        new1_h->set_origin( v );
        new1_h_next->set_origin( mid_v );
        new1_h_prev->set_origin( new_vertex );
        new2_h->set_origin( v );
        new2_h_next->set_origin( new_vertex );
        new2_h_prev->set_origin( v_prev );

        new1_h->set_edge( e1 );
        new1_h_next->set_edge( new3_edge );
        new1_h_prev->set_edge( new1_edge );
        new2_h->set_edge( new1_edge );
        new2_h_next->set_edge( new2_edge );
        new2_h_prev->set_edge( e_prev );

        new1_h->set_face( new1_cell );
        new1_h_next->set_face( new1_cell );
        new1_h_prev->set_face( new1_cell );
        new2_h->set_face( new2_cell );
        new2_h_next->set_face( new2_cell );
        new2_h_prev->set_face( new2_cell );

        // old halfedges
        h->set_origin( mid_v );
        h_prev->set_origin( new_vertex );

        h->set_edge( e2 );
        h_prev->set_edge( new3_edge );

        // new edges
        new1_edge->set_first_halfedge( new1_h_prev );
        new1_edge->set_second_halfedge( new2_h );

        new2_edge->set_first_halfedge( new2_h_next );
        new2_edge->set_second_halfedge( hv_next );

        new3_edge->set_first_halfedge( new1_h_next );
        new3_edge->set_second_halfedge( h_prev );

        // old edges
        e1->set_first_halfedge( he_left );
        e1->set_second_halfedge( new1_h );

        e2->set_first_halfedge( he_right );
        e2->set_second_halfedge( h );

        e_prev->set_first_halfedge( hv_prev );
        e_prev->set_second_halfedge( new2_h_prev );

    }


    // new cells
    new1_cell->set_halfedge( new1_h );
    new2_cell->set_halfedge( new2_h );

    // adding new structures
    _refinMesh->add_vertex( new_vertex );
    _refinMesh->add_edge( new1_edge );
    _refinMesh->add_edge( new2_edge );
    _refinMesh->add_edge( new3_edge );
    _refinMesh->add_face( new1_cell );
    _refinMesh->add_face( new2_cell );

    // returning halfedges

    if ( type == 0 ) {
        he_left = h_prev;
        he_right = new2_h_prev;
    } else {
        he_left = new2_h_next;
        he_right = h_next;
    }

    return true;
}

