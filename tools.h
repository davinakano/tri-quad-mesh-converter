#ifndef TOOLS_H
#define TOOLS_H

#include "types.h"
#include "debug.h"
#include "pack.h"
#include "vtkWriter.h"
#include "stdlib.h"
#include "meshref.h"
#include "geodesicMeshing/geodesicMeshing.h"
#include "floater/floaterpar.hpp"
#include "bezier/tbezier.h"
#include <sstream>
#include "cgalmethods.h"

geodesics* geo = NULL;

#define P_THRESHOLD 1e-6
#define C_THRESHOLD 1e-10

class tools
{

private:
    int _threshold;
public :


    tools(int threshold) {
        _threshold = threshold;
    }

    // ****************************************
    //    FACE <-> QFACE POINTER VERIFICATION
    // ****************************************

    /**
      * \fn void faceQFacePointerVerification( TMesh* tm, QMesh* qm )
      *
      * \brief Verifies the amount of Faces (triangular) per QFace and
      *        compares it to the total of Faces in the original Mesh.
      *
      * \param TMesh Pointer to a triangular mesh.
      *
      * \param QMesh Pointer to a quadrilateral mesh.
      */

    void faceQFacePointerVerification( TMesh* tm, QMesh* qm ) {

        int cont = 0;

        for( QMesh::QFaceIterator it = qm->faces_begin(); it != qm->faces_end(); it++ ) {

            QFace* qf = *it;

            int qcont = 0;

            for( std::set< Face* >::iterator itf = qf->get_attributes().faces.begin(); itf != qf->get_attributes().faces.end(); itf++ ) {

                Face* f = *itf;
                ASSERT ( f->get_attributes().region.size() == 1 );
                ASSERT ( *(f->get_attributes().region.begin()) == qf );

                cont++;
                qcont++;
            }

            ASSERT( qcont == (int)qf->get_attributes().faces.size() );
        }

        for( TMesh::FaceIterator it = tm->faces_begin(); it != tm->faces_end(); it++ ) {

            Face* f = *it;
            ASSERT ( f->get_attributes().region.size() == 1 );
            QFace* qf = *(f->get_attributes().region.begin());

            ASSERT( qf->get_attributes().faces.find( f ) != qf->get_attributes().faces.end() );
        }

        ASSERT (cont == (int)tm->get_number_of_faces());
    }

    // ****************************************
    //              EDGE MARK
    // ****************************************

    // No need to call correctOrientation() for the Half-Edges

    /**
      * \fn void edgeMark( std::list< Vertex* >& lv, QEdge* qe )
      *
      * \brief Every vertex (triangular) in lv has its star checked and
      *        marked either for "qe first half-edge's QFace f1" or the
      *        "qe second half-edge's QFace f2".
      *        This method is recommended for consistent meshes.
      *
      * \param lv List of vertex (triangular) created by the geodesical cut.
      *
      * \param qe QEdge that separates two regions (QFaces).
      */

    void edgeMark( std::list< Vertex* >& lv, QEdge* qe ) {

        QFace* f1 = qe->get_first_halfedge()->get_face() ;
        QFace* f2 = qe->get_second_halfedge()->get_face() ;

        ASSERT( f1 != f2 );


        if( lv.size() == 2 )
        {
            std::list< Vertex* >::iterator iti = lv.begin();
            std::list< Vertex* >::iterator ite = lv.begin();
            ite++;

            Vertex* v = *iti;
            Halfedge* he = v->get_halfedge();

            while ( he->get_next()->get_origin() != *ite ) {
                he = he->get_prev()->get_mate();
                ASSERT(he != v->get_halfedge());
            }

            ASSERT( he->get_face()->get_attributes().region.empty() ||
                    (he->get_face()->get_attributes().region.size() == 1) );

            he->get_face()->get_attributes().region.insert( f1 );

            ASSERT(he->get_face()->get_attributes().region.size() == 1);

            f1->get_attributes().faces.insert( he->get_face() );

            he = he->get_mate();

            ASSERT( he->get_face()->get_attributes().region.empty() ||
                    (he->get_face()->get_attributes().region.size() == 1) );

            he->get_face()->get_attributes().region.insert( f2 );

            ASSERT(he->get_face()->get_attributes().region.size() == 1);

            f2->get_attributes().faces.insert( he->get_face() );


            return;
        }

        ASSERT( lv.size() >= 3 );

        std::list< Vertex* >::iterator it1 = lv.begin();
        std::list< Vertex* >::iterator it2 = it1;
        it2++;
        std::list< Vertex* >::iterator it3 = it2;
        it3++;

        int count = 0;
        while( it3 != lv.end() ) {

            Vertex* gpv = *it2;

            Halfedge* he = gpv->get_halfedge();

            while( he->get_next()->get_origin() != *it1 ){
                he = he->get_prev()->get_mate();
                ASSERT( he->get_origin() == gpv );
                ASSERT( he != gpv->get_halfedge() );
            }

            do {

                ASSERT( he->get_face()->get_attributes().region.empty() ||
                        (he->get_face()->get_attributes().region.size() == 1) );

                he->get_face()->get_attributes().region.insert( f2 );

                ASSERT(he->get_face()->get_attributes().region.size() == 1);

                f2->get_attributes().faces.insert( he->get_face() );

                he = he->get_prev()->get_mate();
                ASSERT( he->get_origin() == gpv );
                ASSERT( he->get_next()->get_origin() != *it1 );
            } while( he->get_next()->get_origin() != *it3 );

            do {

                ASSERT( he->get_face()->get_attributes().region.empty() ||
                        (he->get_face()->get_attributes().region.size() == 1) );

                he->get_face()->get_attributes().region.insert( f1 );

                ASSERT(he->get_face()->get_attributes().region.size() == 1);

                f1->get_attributes().faces.insert( he->get_face() );

                he = he->get_prev()->get_mate();
                ASSERT( he->get_origin() == gpv );
                ASSERT( he->get_next()->get_origin() != *it3 );
            } while( he->get_next()->get_origin() != *it1 );

            he = gpv->get_halfedge();
            do {
                he = he->get_prev()->get_mate();
                ASSERT( ! he->get_face()->get_attributes().region.empty() );
                ASSERT( he->get_origin() == gpv );
            } while( he != gpv->get_halfedge() );

            it1++;
            it2++;
            it3++;
            count++;
        }

    }

    // ****************************************
    //    EDGE MARK (for unconsistent meshes)
    // ****************************************
    // f1 interna
    // f2 externa

    /**
      * \fn void edgeMark( std::list< Vertex* >& lv, QFace *f1, QFace* f2 )
      *
      * \brief Every vertex (triangular) in lv has its star checked and
      *        marked either for internal QFace f1 or the external QFace f2.
      *        This method is recommended for unconsistent meshes.
      *
      * \param lv List of vertex (triangular) created by the geodesical cut.
      *
      * \param f1 Pointer of the internal QFace.
      *
      * \param f2 Pointer of the external QFace.
      */

    void edgeMark( std::list< Vertex* >& lv, QFace *f1, QFace* f2 ) {

        ASSERT( f1 != f2 );


        if( lv.size() == 2 )
        {
            std::list< Vertex* >::iterator iti = lv.begin();
            std::list< Vertex* >::iterator ite = lv.begin();
            ite++;

            Vertex* v = *iti;
            Halfedge* he = v->get_halfedge();

            while ( he->get_next()->get_origin() != *ite ) {
                he = he->get_prev()->get_mate();
                ASSERT(he != v->get_halfedge());
            }

            ASSERT( he->get_face()->get_attributes().region.empty() ||
                    (he->get_face()->get_attributes().region.size() == 1) );

            he->get_face()->get_attributes().region.insert( f1 );

            ASSERT(he->get_face()->get_attributes().region.size() == 1);

            f1->get_attributes().faces.insert( he->get_face() );

            he = he->get_mate();

            ASSERT( he->get_face()->get_attributes().region.empty() ||
                    (he->get_face()->get_attributes().region.size() == 1) );

            he->get_face()->get_attributes().region.insert( f2 );

            ASSERT(he->get_face()->get_attributes().region.size() == 1);

            f2->get_attributes().faces.insert( he->get_face() );


            return;
        }



        ASSERT( lv.size() >= 3 );

        std::list< Vertex* >::iterator it1 = lv.begin();
        std::list< Vertex* >::iterator it2 = it1;
        it2++;
        std::list< Vertex* >::iterator it3 = it2;
        it3++;

        int cont = 0;

        while( it3 != lv.end() ) {

            cont++;
            Vertex* gpv = *it2;

            Halfedge* he = gpv->get_halfedge();

            while( he->get_next()->get_origin() != *it1 ){
                he = he->get_prev()->get_mate();
                ASSERT( he->get_origin() == gpv );
                ASSERT( he != gpv->get_halfedge() );
            }

            do {

                ASSERT( he->get_face()->get_attributes().region.empty() ||
                        (he->get_face()->get_attributes().region.size() == 1) );

                he->get_face()->get_attributes().region.insert( f2 );

                ASSERT(he->get_face()->get_attributes().region.size() == 1);

                f2->get_attributes().faces.insert( he->get_face() );

                he = he->get_prev()->get_mate();
                ASSERT( he->get_origin() == gpv );
                ASSERT( he->get_next()->get_origin() != *it1 );
            } while( he->get_next()->get_origin() != *it3 );

            do {

                ASSERT( he->get_face()->get_attributes().region.empty() ||
                        (he->get_face()->get_attributes().region.size() == 1) );

                he->get_face()->get_attributes().region.insert( f1 );

                ASSERT(he->get_face()->get_attributes().region.size() == 1);

                f1->get_attributes().faces.insert( he->get_face() );

                he = he->get_prev()->get_mate();
                ASSERT( he->get_origin() == gpv );
                ASSERT( he->get_next()->get_origin() != *it3 );
            } while( he->get_next()->get_origin() != *it1 );

            he = gpv->get_halfedge();
            do {
                he = he->get_prev()->get_mate();
                ASSERT( ! he->get_face()->get_attributes().region.empty() );
                ASSERT( he->get_origin() == gpv );
            } while( he != gpv->get_halfedge() );

            it1++;
            it2++;
            it3++;
        }

    }

    // ****************************************
    //              FLOOD REGION
    // ****************************************

    /**
      * \fn void flood_region( Face* tf )
      *
      * \brief Floods a region of triangular faces based on its borders'.
      *
      * \param tf Triangular face that serves as a start point for the flooding.
      */

    void flood_region( Face* tf ) {

        ASSERT( tf->get_attributes().region.empty() );

        std::set< QFace* > region;
        std::set< Face* > faces;
        std::list< Face* > toconsume;

        toconsume.push_back( tf );
        faces.insert( tf );

        while( !toconsume.empty() ) {

            Face* tf = toconsume.front();
            toconsume.erase( toconsume.begin() );

            Halfedge* he = tf->get_halfedge();
            do {
                Halfedge* hemate = he->get_mate();
                ASSERT( hemate->get_face()->get_attributes().region.empty() || hemate->get_face()->get_attributes().region.size() == 1 );

                if( hemate->get_face()->get_attributes().region.empty() ) {

                    if( faces.find( hemate->get_face() ) == faces.end() ) {
                        toconsume.push_back( hemate->get_face() );
                        faces.insert( hemate->get_face() );
                    }

                } else {

                    if(region.empty()) {
                        region = hemate->get_face()->get_attributes().region;
                    } else {

                        ASSERT( region.size() == 1 );

                        ASSERT( *(region.begin()) == *(hemate->get_face()->get_attributes().region.begin()) );
                    }
                }
                he = he->get_next();
            } while( he != tf->get_halfedge() );

        }

        ASSERT( region.size() == 1 );
        QFace* r = *region.begin();

        for( std::set< Face* >::iterator it = faces.begin(); it!= faces.end(); it++ ) {
            Face* tf = *it;
            tf->get_attributes().region.insert( r );
        }

        r->get_attributes().faces.insert( faces.begin(), faces.end() );

        return;
    }

    void getGeodesics( QFace* qf, Vertex* vi, Vertex*ve, std::list< Vertex* >& lv, std::set< pack3<Vertex*,double,double> >& vertices, geodesicMeshing* geo ) {

        geodesics::GP* gp_ini = new geodesics::GPV( vi );
        geodesics::GP* gp_end = new geodesics::GPV( ve );

        std::set< pack3<Vertex*,double,double> >::iterator it = vertices.find( pack3<Vertex*,double,double>( vi, 0.0, 0.0) );
        ASSERT(it != vertices.end() );
        pack3<Vertex*,double,double> pit = *it;
        double u_ini = pit.second;
        double v_ini = pit.third;

        it = vertices.find( pack3<Vertex*,double,double>( ve, 0.0, 0.0) );
        ASSERT(it != vertices.end() );
        pit = *it;
        double u_end = pit.second;
        double v_end = pit.third;

        if(debug::iteration==6 && debug::templatenumber==784)
            int a = 0;

        std::list<geodesics::GP*> gplist;
        geodesics2D(qf->get_attributes().faces, vertices,
                    gp_ini, u_ini, v_ini,
                    gp_end, u_end, v_end,
                    gplist);



        std::list< std::pair<double,double> > gplist_uv;
        for( std::list<geodesics::GP*>::iterator it = gplist.begin(); it != gplist.end(); it++) {
            std::pair<double,double> p;
            geodesics::GP* gp = *it;
            if( gp->on_vertex()){
                geodesics::GPV* gpv = dynamic_cast< geodesics::GPV* >(gp) ;
                std::set< pack3<Vertex*,double,double> >::iterator itv = vertices.find( pack3<Vertex*,double,double>( gpv->get_vertex(), 0.0, 0.0) );
                ASSERT(itv != vertices.end() );
                pack3<Vertex*,double,double> pit = *itv;
                p.first = pit.second;
                p.second = pit.third;
            } else {
                geodesics::GPE* gpe = dynamic_cast< geodesics::GPE* >(gp) ;
                std::set< pack3<Vertex*,double,double> >::iterator itv = vertices.find( pack3<Vertex*,double,double>( gpe->get_edge()->get_first_halfedge()->get_origin(), 0.0, 0.0) );
                ASSERT(itv != vertices.end() );
                pack3<Vertex*,double,double> pit0 = *itv;
                itv = vertices.find( pack3<Vertex*,double,double>( gpe->get_edge()->get_second_halfedge()->get_origin(), 0.0, 0.0) );
                ASSERT(itv != vertices.end() );
                pack3<Vertex*,double,double> pit1 = *itv;

                p.first = (1.0-gpe->get_t())*pit0.second + gpe->get_t()*pit1.second;
                p.second = (1.0-gpe->get_t())*pit0.third + gpe->get_t()*pit1.third;
            }
            gplist_uv.push_back( p );
        }

        geo->remeshingGeodesics(gplist, lv);

        ASSERT( gplist_uv.size() == lv.size()+2);

        std::list< std::pair<double,double> >::iterator it_gp = gplist_uv.begin();
        it_gp++;
        for( std::list< Vertex* >::iterator it = lv.begin(); it != lv.end(); it++) {
            //ASSERT( vertices.count( pack3<Vertex*,double,double>(*it,0,0) ) == 0 );
            vertices.insert( pack3<Vertex*,double,double>( *it, it_gp->first, it_gp->second ));
            it_gp++;
        }
        it_gp++;
        ASSERT( it_gp == gplist_uv.end());

        for( std::list<geodesics::GP*>::iterator it = gplist.begin(); it != gplist.end(); it++) {
            geodesics::GP* gp = *it;
            delete gp;
        }



        lv.push_front(vi);
        lv.push_back(ve);

        for( std::list< Vertex* >::iterator it = lv.begin(); it != lv.end(); it++) {

            std::list< Vertex* >::iterator it2 = it;
            it2++;

            if( it2 != lv.end()) {

                Vertex* v1 = *it;
                Vertex* v2 = *it2;
                Halfedge* he = v1->get_halfedge();
                while( he->get_next()->get_origin() != v2 ) {
                    he = he->get_prev()->get_mate();
                    ASSERT( he->get_origin() == v1);
                    ASSERT( he != v1->get_halfedge() );
                }
            }
        }

        lv.pop_back();
        lv.pop_front();
    }

    double intersect( double v0x, double v0y, double v1x, double v1y, double p0x, double p0y, double p1x, double p1y ) {

        double tv, tp;
        common::Geometric::INTERSECTION_TYPE ret =
                common::Geometric::compute_segment_intersection(
                    v0x, v0y,
                    v1x, v1y,
                    p0x, p0y,
                    p1x, p1y,
                    tv ,
                    tp) ;

        ASSERT( ( ret == common::Geometric::POINT ) || ( ret == common::Geometric::NOINTERSECTION ) ) ;

        return tp;

    }

    double dist2( double v0x, double v0y, double v1x, double v1y ) {

        return (v1x - v0x)*(v1x - v0x) + (v1y - v0y)*(v1y - v0y);
    }

    int inLeft(double v0x, double v0y, double v1x, double v1y, double px, double py){

        // Retorna 1 caso esteja a esquerda
        // Retorna 0 caso esteja a direita
        // Retorna 2 caso esteja sobre

        double i = ( (v1x - v0x)*(py - v0y) - (px - v0x)*(v1y - v0y) );

        if (i > P_THRESHOLD)
            return 1;
        else
        {
            if(i < - P_THRESHOLD)
                return 0;
            else
                return 2;
        }
    }

    double get_coord( double ini, double end, double t)
    {
        return ini + t*(end-ini);
    }


    void geodesics2D( std::set<Face*>& faces, std::set< pack3<Vertex*,double,double> >& vertices,
                      geodesics::GP* gp_ini, double u_ini, double v_ini,
                      geodesics::GP* gp_end, double u_end, double v_end,
                      std::list< geodesics::GP* >& gplist)
    {
        ASSERT( !faces.empty() );
        ASSERT( vertices.size() > 2 );
        ASSERT( gp_ini );
        ASSERT( gp_end );
        ASSERT( gp_ini != gp_end );
        ASSERT( gplist.empty() );

        geodesics::GP* gp = gp_ini;
        Halfedge* he;

        if( gp->on_vertex() ) {
            he = 0;
        } else {
            ASSERT( gp->on_edge() );
            geodesics::GPE* gpe = dynamic_cast< geodesics::GPE* >(gp) ;
            Edge* e = gpe->get_edge();
            std::set< pack3<Vertex*,double,double> >::iterator it = vertices.find( pack3<Vertex*,double,double>( e->get_first_halfedge()->get_origin(), 0.0, 0.0) );
            ASSERT(it != vertices.end() );
            double u0 = it->second;
            double v0 = it->third;

            it = vertices.find( pack3<Vertex*,double,double>(e->get_second_halfedge()->get_origin(),0,0) );
            ASSERT(it != vertices.end() );
            double u1 = it->second;
            double v1 = it->third;

            if(inLeft( u0,v0,u1,v1, u_end, v_end) == 1 )
                he = e->get_first_halfedge();
            else
                he = e->get_second_halfedge();
        }

        double u_atu = u_ini;
        double v_atu = v_ini;

        gplist.push_back( gp_ini );

        while( (fabs(u_atu-u_end)>C_THRESHOLD )||
               (fabs(v_atu-v_end)>C_THRESHOLD )
               ) {


            if( gp->on_vertex() ) {
                ASSERT( !he );
                geodesics::GPV* gpv = dynamic_cast< geodesics::GPV* >(gp) ;
                bool cont = true;
                Halfedge* he_temp = gpv->get_vertex()->get_halfedge();
                do {

                    if( faces.count(he_temp->get_face()) == 1 ) {

                        Vertex* vtx0 = he_temp->get_next()->get_origin();
                        Vertex* vtx1 = he_temp->get_prev()->get_origin();

                        std::set< pack3<Vertex*,double,double> >::iterator it = vertices.find( pack3<Vertex*,double,double>( vtx0, 0.0, 0.0) );
                        ASSERT(it != vertices.end() );
                        double u0 = it->second;
                        double v0 = it->third;

                        it = vertices.find( pack3<Vertex*,double,double>(vtx1,0,0) );
                        ASSERT(it != vertices.end() );
                        double u1 = it->second;
                        double v1 = it->third;

                        double t01;
                        int tipo;
                        cgalMethods::intersection(u0, v0, u1, v1, u_ini, v_ini, u_end, v_end, tipo, t01);
                        if( tipo == 0 ) { // t01 == 0
                            // passa por vtx0
                            if( dist2( u0, v0, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) ) {
                                geodesics::GPV* gpv_new = new geodesics::GPV( vtx0 );
                                gplist.push_back( gpv_new );
                                he = 0;
                                u_atu = u0;
                                v_atu = v0;
                                cont = false;
                                gp = gpv_new;
                            }
                        } else if( tipo == 1 ) {
                            // passa por vtx1
                            if( dist2( u1, v1, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) ) {
                                geodesics::GPV* gpv_new = new geodesics::GPV( vtx1 );
                                gplist.push_back( gpv_new );
                                he = 0;
                                u_atu = u1;
                                v_atu = v1;
                                cont = false;
                                gp = gpv_new;
                            }
                        } else if( tipo == 2) {
                            // passa entre vtx0 e vtx1
                            if( dist2( get_coord(u0, u1, t01), get_coord(v0, v1, t01), u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) ) {
                                geodesics::GPE* gpe_new;
                                if( he_temp->get_next()->get_edge()->get_first_halfedge() == he_temp->get_next())
                                    gpe_new = new geodesics::GPE( he_temp->get_next()->get_edge(), t01);
                                else
                                    gpe_new = new geodesics::GPE( he_temp->get_next()->get_edge(), 1.0 - t01);

                                gplist.push_back( gpe_new );
                                he = he_temp->get_next()->get_mate();
                                u_atu = get_coord(u0, u1, t01);
                                v_atu = get_coord(v0, v1, t01);
                                cont = false;
                                gp = gpe_new;
                            }
                        } else if( tipo == 4 ) {
                        } else {
                            ASSERT( false );
                        }



                        //                        if( fabs( t01 ) < P_THRESHOLD ) { // t01 == 0
                        //                            // passa por v0
                        //                            if( dist2( u0, v0, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) ) {
                        //                                geodesics::GPV* gpv_new = new geodesics::GPV( vtx0 );
                        //                                gplist.push_back( gpv_new );
                        //                                he = 0;
                        //                                u_atu = u0;
                        //                                v_atu = v0;
                        //                                cont = false;
                        //                                gp = gpv_new;
                        //                            }
                        //                        } else if( fabs( t01 - 1.0 ) < P_THRESHOLD ) {
                        //                            // passa por v2
                        //                            if( dist2( u1, v1, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) ) {
                        //                                geodesics::GPV* gpv_new = new geodesics::GPV( vtx1 );
                        //                                gplist.push_back( gpv_new );
                        //                                he = 0;
                        //                                u_atu = u1;
                        //                                v_atu = v1;
                        //                                cont = false;
                        //                                gp = gpv_new;
                        //                            }
                        //                        } else if( (t01>0)&&(t01<1.0)) {
                        //                            // passa entre v0 e v2
                        //                            if( dist2( get_coord(u0, u1, t01), get_coord(v0, v1, t01), u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) ) {
                        //                                geodesics::GPE* gpe_new;
                        //                                if( he_temp->get_next()->get_edge()->get_first_halfedge() == he_temp->get_next())
                        //                                    gpe_new = new geodesics::GPE( he_temp->get_next()->get_edge(), t01);
                        //                                else
                        //                                    gpe_new = new geodesics::GPE( he_temp->get_next()->get_edge(), 1.0 - t01);

                        //                                gplist.push_back( gpe_new );
                        //                                he = he_temp->get_next()->get_mate();
                        //                                u_atu = get_coord(u0, u1, t01);
                        //                                v_atu = get_coord(v0, v1, t01);
                        //                                cont = false;
                        //                                gp = gpe_new;
                        //                            }
                        //                        }
                    }
                    he_temp = he_temp->get_prev()->get_mate();
                    ASSERT( he_temp->get_origin() == gpv->get_vertex() );
                    ASSERT( (!cont)||(he_temp != gpv->get_vertex()->get_halfedge()) );
                } while( cont );


            } else {
                ASSERT( gp->on_edge() );
                ASSERT( he );


                Vertex* vtx0 = he->get_origin();
                Vertex* vtx1 = he->get_next()->get_origin();
                Vertex* vtx2 = he->get_prev()->get_origin();

                std::set< pack3<Vertex*,double,double> >::iterator it = vertices.find( pack3<Vertex*,double,double>( vtx0, 0.0, 0.0) );
                ASSERT(it != vertices.end() );
                double u0 = it->second;
                double v0 = it->third;

                it = vertices.find( pack3<Vertex*,double,double>(vtx1,0,0) );
                ASSERT(it != vertices.end() );
                double u1 = it->second;
                double v1 = it->third;

                it = vertices.find( pack3<Vertex*,double,double>(vtx2,0,0) );
                ASSERT(it != vertices.end() );
                double u2 = it->second;
                double v2 = it->third;

                //double t12 = intersect(u_ini, v_ini, u_end, v_end, u1, v1, u2, v2);
                double t12;
                int tipo;
                cgalMethods::intersection(u1, v1, u2, v2, u_ini, v_ini, u_end, v_end, tipo, t12);
                if( tipo == 0 ) { // t12 == 0
                    ASSERT( false);
                }
                else if( tipo == 1 )
                {
                    // passa por v2
                    geodesics::GPV* gpv_new = new geodesics::GPV( vtx2 );
                    gplist.push_back( gpv_new );
                    he = 0;
                    ASSERT( dist2( u2, v2, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                    u_atu = u2;
                    v_atu = v2;
                    gp = gpv_new;
                } else if( tipo == 2) {
                    // passa entre v1 e v2

                    geodesics::GPE* gpe_new;
                    if( he->get_next()->get_edge()->get_first_halfedge() == he->get_next() )
                        gpe_new = new geodesics::GPE( he->get_next()->get_edge(), t12);
                    else
                        gpe_new = new geodesics::GPE( he->get_next()->get_edge(), 1.0 - t12);

                    gplist.push_back( gpe_new );
                    he = he->get_next()->get_mate();
                    //ASSERT( dist2( get_coord(u1, u2, t12), get_coord(v1, v2, t12), u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                    u_atu = get_coord(u1, u2, t12);
                    v_atu = get_coord(v1, v2, t12);
                    gp = gpe_new;
                } else if( tipo == 4 ) {
                    // double t02 = intersect(u_ini, v_ini, u_end, v_end, u0, v0, u2, v2);
                    double t02;
                    int tipo;
                    cgalMethods::intersection(u0, v0, u2, v2, u_ini, v_ini, u_end, v_end, tipo, t02);
                    if( tipo == 0 ) {
                        ASSERT( false );
                    } else if( tipo == 1 ) {
                        // passa por v2
                        geodesics::GPV* gpv_new = new geodesics::GPV( vtx2 );
                        gplist.push_back( gpv_new );
                        he = 0;
                        ASSERT( dist2( u2, v2, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                        u_atu = u2;
                        v_atu = v2;
                        gp = gpv_new;
                    } else if ( tipo == 2 ){
                        ASSERT( (t02>0)&&(t02<1.0));
                        // passa entre v0 e v2
                        geodesics::GPE* gpe_new;
                        if( he->get_prev()->get_edge()->get_first_halfedge() == he->get_prev() )
                            gpe_new = new geodesics::GPE( he->get_prev()->get_edge(), 1.0 - t02);
                        else
                            gpe_new = new geodesics::GPE( he->get_prev()->get_edge(), t02);

                        gplist.push_back( gpe_new );
                        he = he->get_prev()->get_mate();
                        ASSERT( dist2( get_coord(u0, u2, t02), get_coord(v0, v2, t02), u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                        u_atu = get_coord(u0, u2, t02);
                        v_atu = get_coord(v0, v2, t02);
                        gp = gpe_new;
                    } else if( tipo == 4 ) {
                    } else {
                        ASSERT( false );
                    }
                } else {
                    ASSERT( false );
                }

                //                if( fabs( t12 ) < P_THRESHOLD ) // passa por v1 - nao permitido
                //                {
                //                    ASSERT( t12 > 0 );
                //                    t12 = P_THRESHOLD;
                //                }
                //                if( fabs( t12 - 1.0 ) < P_THRESHOLD )
                //                {
                //                    // passa por v2
                //                    geodesics::GPV* gpv_new = new geodesics::GPV( vtx2 );
                //                    gplist.push_back( gpv_new );
                //                    he = 0;
                //                    ASSERT( dist2( u2, v2, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                //                    u_atu = u2;
                //                    v_atu = v2;
                //                    gp = gpv_new;
                //                } else if( (t12>0)&&(t12<1.0)) {
                //                    // passa entre v1 e v2

                //                    geodesics::GPE* gpe_new;
                //                    if( he->get_next()->get_edge()->get_first_halfedge() == he->get_next() )
                //                        gpe_new = new geodesics::GPE( he->get_next()->get_edge(), t12);
                //                    else
                //                        gpe_new = new geodesics::GPE( he->get_next()->get_edge(), 1.0 - t12);

                //                    gplist.push_back( gpe_new );
                //                    he = he->get_next()->get_mate();
                //                    //ASSERT( dist2( get_coord(u1, u2, t12), get_coord(v1, v2, t12), u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                //                    u_atu = get_coord(u1, u2, t12);
                //                    v_atu = get_coord(v1, v2, t12);
                //                    gp = gpe_new;
                //                } else {
                //                    double t02 = intersect(u_ini, v_ini, u_end, v_end, u0, v0, u2, v2);

                //                    if( fabs( t02 ) < P_THRESHOLD ) // passa por v0 - nao permitido
                //                    {
                //                        ASSERT( t02 > 0 );
                //                        t02 = P_THRESHOLD;
                //                    }


                //                    if( fabs( t02 - 1.0 ) < P_THRESHOLD )
                //                    {
                //                        // passa por v2
                //                        geodesics::GPV* gpv_new = new geodesics::GPV( vtx2 );
                //                        gplist.push_back( gpv_new );
                //                        he = 0;
                //                        ASSERT( dist2( u2, v2, u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                //                        u_atu = u2;
                //                        v_atu = v2;
                //                        gp = gpv_new;
                //                    } else {
                //                        ASSERT( (t02>0)&&(t02<1.0));
                //                        // passa entre v0 e v2
                //                        geodesics::GPE* gpe_new;
                //                        if( he->get_prev()->get_edge()->get_first_halfedge() == he->get_prev() )
                //                            gpe_new = new geodesics::GPE( he->get_prev()->get_edge(), 1.0 - t02);
                //                        else
                //                            gpe_new = new geodesics::GPE( he->get_prev()->get_edge(), t02);

                //                        gplist.push_back( gpe_new );
                //                        he = he->get_prev()->get_mate();
                //                        ASSERT( dist2( get_coord(u0, u2, t02), get_coord(v0, v2, t02), u_end, v_end) < dist2( u_atu, v_atu, u_end, v_end) );
                //                        u_atu = get_coord(u0, u2, t02);
                //                        v_atu = get_coord(v0, v2, t02);
                //                        gp = gpe_new;
                //                    }
                //                }
            }

        }

        gplist.pop_back();
        gplist.push_back( gp_end );
    }


    std::list< Vertex* > getBoundary( QVertex* qv0, QVertex* qv1, QFace* qf) {

        ASSERT( qv0->get_attributes().getGP()->on_vertex() );
        ASSERT( qv1->get_attributes().getGP()->on_vertex() );
        geodesics::GPV* gpv0 = dynamic_cast< geodesics::GPV*>(qv0->get_attributes().getGP() );
        geodesics::GPV* gpv1 = dynamic_cast< geodesics::GPV*>(qv1->get_attributes().getGP() );
        Vertex* v0 = gpv0->get_vertex();
        Vertex* v1 = gpv1->get_vertex();

        std::list< Vertex* > lv;
        lv.push_back( v0 );

        do {
            Halfedge* he = v0->get_halfedge();

            ASSERT(he->get_face()->get_attributes().region.size() == 1);

            QFace* qfhe = *(he->get_face()->get_attributes().region.begin());
            QFace* qfmate = *(he->get_mate()->get_face()->get_attributes().region.begin());
            while( !((qfhe == qf)&&(qfmate != qf)) ) {
                he = he->get_prev()->get_mate();

                qfhe = *(he->get_face()->get_attributes().region.begin());
                qfmate = *(he->get_mate()->get_face()->get_attributes().region.begin());

                ASSERT( he->get_origin() == v0);
                ASSERT( he != v0->get_halfedge() );
            }

            he = he->get_next();
            v0 = he->get_origin();
            lv.push_back( v0 );

        } while ( v0 != v1 );

        return lv;
    }


    // ****************************************
    //              FLOATER
    // ****************************************

    /**
      * \fn void floater( QFace* qf )
      *
      * \brief Parameterizes a region given by a QFace. If it's on debug mode,
      *        it saves both 2-D and 3-D patches.
      *
      * \param qf Quadrilateral face that will be parameterized.
      */

    void floater( QFace* qf, std::set< pack3<Vertex*,double,double> >& vertices_uv ) {

        ASSERT( qf );
        ASSERT( ! qf->get_attributes().faces.empty() );
        ASSERT( vertices_uv.empty() );

        std::set< Vertex* > vertices;

        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {

            Face* f = *it;
            Vertex* v0 = f->get_halfedge()->get_origin();
            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

            vertices.insert( v0 );
            vertices.insert( v1 );
            vertices.insert( v2 );
        }

        double* vset = new double[3*vertices.size()];
        unsigned i = 0;
        for( std::set< Vertex* >::iterator it = vertices.begin(); it != vertices.end(); it++ ) {
            Vertex* v = *it;
            v->get_attributes().tempInt = i;
            v->get_attributes().tempDouble = -1;

            vset[3*i] = v->x();
            vset[3*i+1] = v->y();
            vset[3*i+2] = v->z();

            i++;
        }

        i = 0;
        unsigned* fset = new unsigned[ 3*qf->get_attributes().faces.size()];
        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {
            Face* f = *it;
            Vertex* v0 = f->get_halfedge()->get_origin();
            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

            fset[3*i] = v0->get_attributes().tempInt;
            fset[3*i+1] = v1->get_attributes().tempInt;
            fset[3*i+2] = v2->get_attributes().tempInt;
            i++;
        }

        double offset = 0.0;
        QHalfedge *he = qf->get_halfedge();
        do {

            std::list< Vertex* > geoline = getBoundary(he->get_origin(), he->get_next()->get_origin(), qf);

            ASSERT( geoline.size() >= 2 );

            std::list< double > dists;

            double disttotal = 0;
            std::list< Vertex* >::iterator it1 = geoline.begin();
            std::list< Vertex* >::iterator it2 = it1;
            it2++;

            while(it2 != geoline.end()) {

                Vertex* v1 = *it1;
                Vertex* v2 = *it2;

                ASSERT( vertices.find( v1 ) != vertices.end() );
                ASSERT( vertices.find( v2 ) != vertices.end() );
                double d = common::Geometric::distance_3d( v1->x(), v1->y(), v1->z(), v2->x(), v2->y(), v2->z() );
                dists.push_back( d );
                disttotal += d;

                it1++;
                it2++;
            }

            it1 = geoline.begin();
            Vertex* v  = *it1;
            v->get_attributes().tempDouble = offset;
            it1++;
            dists.pop_back();
            double distatual = 0;
            for( std::list< double >::iterator it = dists.begin(); it!= dists.end(); it++ ) {
                distatual += *it;
                Vertex* v  = *it1;
                v->get_attributes().tempDouble = offset + distatual/disttotal;
                it1++;
            }
            it1++;
            ASSERT( it1 == geoline.end() );

            offset += 1.0;

            he = he->get_next();
        }while( he != qf->get_halfedge() );

        i = 0;
        double* tvalue = new double[vertices.size()];
        for( std::set< Vertex* >::iterator it = vertices.begin(); it != vertices.end(); it++ ) {

            Vertex* v = *it;
            tvalue[i] = v->get_attributes().tempDouble;
            i++;
        }

        floaterpar::FloaterPar fp( vertices.size(), vset, qf->get_attributes().faces.size(), fset, tvalue);
        fp.parametrize();

        //#ifdef DEBUG_PATCH
        //        vtkReader reader;

        //        std::string filename = "debug-patch";
        //        std::stringstream str;
        //        str << qf->get_attributes().id;

        //        std::string filename1 = filename + std::string("-3D-") + str.str();
        //        reader.saveVTKTriPatch( filename1, vertices.size(), vset, qf->get_attributes().faces.size(), fset, tvalue);

        //        double *uv_ = new double[3*vertices.size()];
        //        for( i=0; i < vertices.size(); i++ ) {

        //            uv_[3*i] = fp.get_ucoord(i);
        //            uv_[3*i+1] = fp.get_vcoord(i);
        //            uv_[3*i+2] = 0.0;
        //        }

        //        std::string filename2 = filename + std::string("-2D-") + str.str();
        //        reader.saveVTKTriPatch( filename2, vertices.size(), uv_, qf->get_attributes().faces.size(), fset, tvalue);

        //        delete [] uv_;
        //#endif

        double *uv = new double[2*vertices.size()];
        for( i=0; i < vertices.size(); i++ ) {

            uv[2*i] = fp.get_ucoord(i);
            uv[2*i+1] = fp.get_vcoord(i);
        }

        //this->error_metric(uv, vset, 3*vertices.size(), qf, 0);

        //bezier::tBezier bezier( uv, vset, vertices.size(), 3, 3, 0, 1, 0, 1);


        i = 0;
        for( std::set< Vertex* >::iterator it = vertices.begin(); it != vertices.end(); it++ ) {
            Vertex* v = *it;
            vertices_uv.insert( pack3<Vertex*,double,double>( v, uv[2*i], uv[2*i+1]) );
            i++;
        }


        delete [] uv;
        delete [] vset;
        delete [] fset;
        delete [] tvalue;
    }

    // ****************************************
    //              INITIAL SUBDIVISION
    // ****************************************

    /**
      * \fn void initialDivision( TMesh* tm, QMesh* qm, geodesicMeshing* geo )
      *
      * \brief It makes the first subdivisions of a given mesh, generating
      *        6 quads (as a cube) that correspond to 6 regions of the triangle
      *        mesh.
      *
      * \param tm Pointer of a triangular mesh.
      *
      * \param qm Pointer of a quadrilateral mesh.
      *
      * \param geo Pointer of a geodesic object that have cut methods.
      */

    void initialDivision( TMesh* tm, QMesh* qm, geodesicMeshing* geo ) {

        ASSERT( tm );
        ASSERT( qm );
        ASSERT( geo );

        for( QMesh::QEdgeIterator it = qm->edges_begin() ;
             it != qm->edges_end(); it++ ) {

            QEdge* qe = *it;

            ASSERT( qe->get_attributes().geoline.empty() );

            /*
            geodesics::GP* v1 = qe->get_first_halfedge()->get_origin()->get_attributes().getGP();
            geodesics::GP* v2 = qe->get_second_halfedge()->get_origin()->get_attributes().getGP();

            std::list< geodesics::GP* > lv;

            geo->setRegion( 0 );
            geo->compute( v1, v2, lv );
            */

            geodesics::GP* gp1 = qe->get_first_halfedge()->get_origin()->get_attributes().getGP();
            geodesics::GP* gp2 = qe->get_second_halfedge()->get_origin()->get_attributes().getGP();

            ASSERT( gp1->on_vertex() );
            ASSERT( gp2->on_vertex() );

            geodesics::GPV* gpv1 = dynamic_cast< geodesics::GPV* >( gp1 ) ;
            geodesics::GPV* gpv2 = dynamic_cast< geodesics::GPV* >( gp2 ) ;

            Vertex* v1 = gpv1->get_vertex();
            Vertex* v2 = gpv2->get_vertex();

            std::list< Vertex* > lv;

            geo->getGeodesic(0,v1,v2,lv);

            ASSERT( qm->is_consistent() );

            ASSERT( lv.size() > 0 );

            lv.push_front( v1 );
            lv.push_back( v2 );

            qe->get_attributes().geoline = lv;
        }

        for( QMesh::QEdgeIterator it = qm->edges_begin() ;
             it != qm->edges_end(); it++ ) {

            QEdge* qe = *it;

            ASSERT( !qe->get_attributes().geoline.empty() );
            edgeMark( qe->get_attributes().geoline, qe );
        }

        for( TMesh::FaceIterator it = tm->faces_begin() ;
             it != tm->faces_end(); it++ ) {

            Face* tf = *it;

            if( tf->get_attributes().region.empty() ) {

                flood_region( tf);
            }

        }

        return;
    }

    // ****************************************
    //              UPDATE SUBDIVISION
    // ****************************************

    //    void updateDivision( QHalfedge* qhe ) {

    //        QEdge* qenext = qhe->get_next()->get_edge();

    //        ASSERT( qenext->get_attributes().geoline.empty() );

    //        geodesics::GP* v1 = qenext->get_first_halfedge()->get_origin()->get_attributes().getGP();
    //        geodesics::GP* v2 = qenext->get_second_halfedge()->get_origin()->get_attributes().getGP();

    //        std::list< geodesics::GP* > lv;

    //        geo->setRegion( 0 );
    //        geo->compute( v1, v2, lv );

    //        qenext->get_attributes().geoline = lv;


    //        edgeMark( qenext->get_attributes().geoline, qenext );


    //    }

    // ****************************************
    //              GET BARICENTRIC X,Y,Z
    // ****************************************

    /**
      * \fn void PPS< Mesh >::from_barycentric_to_Cartesian( double u ,
    double v , double w , double x0 , double y0 , double x1 , double y1 ,
    double x2 , double y2 , double& x , double& y ) const
      *
      * \brief  Converts the  barycentric  coordinates of  a point,  with
      * respect to a reference triangle, to Cartesian coordinates.
      *
      * \param u First barycentric coordinate of the point.
      * \param v Second barycentric coordinate of the point.
      * \param w Third barycentric coordinate of the point.
      * \param x0 First  Cartesian coordinate of the first  vertex of the
      * reference triangle.
      * \param y0 Second Cartesian coordinate  of the first vertex of the
      * reference triangle.
      * \param x1 First Cartesian coordinate  of the second vertex of the
      * reference triangle.
      * \param y1 Second Cartesian coordinate of the second vertex of the
      * reference triangle.
      * \param x2 First  coordinate of the third vertex  of the reference
      * triangle.
      * \param y2 Second coordinate of  the third vertex of the reference
      * triangle.
      * \param x First Cartesian coordinate of the resulting point.
      * \param y Second Cartesian coordinate of the resulting point.
      *
      */
    void from_barycentric_to_Cartesian(
            double u ,
            double v ,
            double w ,
            double x0 ,
            double y0 ,
            double z0 ,
            double x1 ,
            double y1 ,
            double z1 ,
            double x2 ,
            double y2 ,
            double z2 ,
            double& x ,
            double& y ,
            double& z
            )
    const
    {
        x = ( u * x0 ) + ( v * x1 ) + ( w * x2 ) ;
        y = ( u * y0 ) + ( v * y1 ) + ( w * y2 ) ;
        z = ( u * z0 ) + ( v * z1 ) + ( w * z2 ) ;
    }

    // ****************************************
    //              GET BARICENTRIC UV
    // ****************************************

    /**
       * \fn void PPS< Mesh >::get_barycentric_coordinates(double x0 ,
     double y0 , double x1 , double y1 , double x2 , double y2 , double xp
     , double yp , double& u , double& v , double& w ) const
       *
       * \brief Computes the barycentric  coordinates of a given point (in
       * Cartesian  coordinates)   with  respect  to   a  given  reference
       * triangle.
       *
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
       * \param xp First Cartesian coordinate of the point.
       * \param yp Second Cartesian coordinate of the point.
       * \param u First barycentric coordinate of the point.
       * \param v Second barycentric coordinate of the point.
       * \param w Third barycentric coordinate of the point.
       *
       */
    void get_barycentric_coordinates(
            double x0 ,
            double y0 ,
            double x1 ,
            double y1 ,
            double x2 ,
            double y2 ,
            double xp ,
            double yp ,
            double& u ,
            double& v ,
            double& w
            )
    const
    {
        /**
         * Compute the determinant.
         */
        double dd = ( x1 * y0 ) - ( x2 * y0 ) - ( x0 * y1 )
                + ( x2 * y1 ) + ( x0 * y2 ) - ( x1 * y2 ) ;

        /**
         * The determinant cannot be zero.
         */
        assert( fabs( dd ) > 1e-16 ) ;

        /**
         * Compute the barycentric coordinates.
         */
        u = ( x2 * y1 ) - ( xp * y1 ) - ( x1 * y2 )
                + ( xp * y2 ) + ( x1 * yp ) - ( x2 * yp ) ;

        u /= dd ;

        v = ( xp * y0 ) - ( x2 * y0 ) + ( x0 * y2 )
                - ( xp * y2 ) - ( x0 * yp ) + ( x2 * yp ) ;

        v /= dd ;

        if ( fabs( u ) < 1e-14 ) {
            u = 0 ;
        }
        else if ( fabs( 1 - u ) < 1e-14 ) {
            u = 1 ;
        }

        if ( fabs( v ) < 1e-14 ) {
            v = 0 ;
        }
        else if ( fabs( 1 - v ) < 1e-14 ) {
            v = 1 ;
        }

        w = 1 - u - v ;

        if ( fabs( w ) < 1e-14 ) {
            w = 0 ;
        }
        else if ( fabs( 1 - w ) < 1e-14 ) {
            w = 1 ;
        }

        return ;
    }

    // ****************************************
    //            ERROR METRIC
    // ****************************************

    /**
       * \fn void error_metric(double *uv, double *vset, int size, QFace* qf, double errorMetric)
       *
       * \brief Computes the barycentric  coordinates of a given point (in
       * Cartesian  coordinates)   with  respect  to   a  given  reference
       * triangle.
       *
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
       *
       */

    void error_metric(double *uv, double *vset, int size, QFace* qf, double errorMetric)
    {
        double newU, newV, newW,
                newX, newY, newZ;

        double *newCoordSet = new double[3*size];

        double *P0 = new double[3];

        P0[0] = qf->get_halfedge()->get_origin()->x(); // P0x
        P0[1] = qf->get_halfedge()->get_origin()->y(); // P0y
        P0[2] = qf->get_halfedge()->get_origin()->z(); // P0z

        double *P1 = new double[3];

        P1[0] = qf->get_halfedge()->get_next()->get_origin()->x(); // P1x
        P1[1] = qf->get_halfedge()->get_next()->get_origin()->y(); // P1y
        P1[2] = qf->get_halfedge()->get_next()->get_origin()->z(); // P1z

        double *P2 = new double[3];

        P2[0] = qf->get_halfedge()->get_next()->get_next()->get_origin()->x(); // P2x
        P2[1] = qf->get_halfedge()->get_next()->get_next()->get_origin()->y(); // P2y
        P2[2] = qf->get_halfedge()->get_next()->get_next()->get_origin()->z(); // P2z

        double *P3 = new double[3];

        P3[0] = qf->get_halfedge()->get_prev()->get_origin()->x(); // P3x
        P3[1] = qf->get_halfedge()->get_prev()->get_origin()->y(); // P3y
        P3[2] = qf->get_halfedge()->get_prev()->get_origin()->z(); // P3z

        for ( int i = 0 ; i < size ; i++ ){

            // if u coordinate is greater than v coordinate it belongs to the triangle P0, P1, P2
            if (uv[2*i] > uv[2*i+1]){

                this->get_barycentric_coordinates(0,0,1,0,1,1,uv[2*i],uv[2*i+1],newU,newV,newW);
                this->from_barycentric_to_Cartesian(newU,newV,newW,
                                                    P3[0], P3[1], P3[2],
                                                    P0[0], P0[1], P0[2],
                                                    P1[0], P1[1], P1[2],
                                                    newX,newY,newZ);

            } else {

                this->get_barycentric_coordinates(0,0,1,1,0,1,uv[2*i],uv[2*i+1],newU,newV,newW);
                this->from_barycentric_to_Cartesian(newU,newV,newW,
                                                    P3[0], P3[1], P3[2],
                                                    P1[0], P1[1], P1[2],
                                                    P2[0], P2[1], P2[2],
                                                    newX,newY,newZ);

            }

            newCoordSet[3*i] = newX;
            newCoordSet[3*i+1] = newY;
            newCoordSet[3*i+2] = newZ;

            // Calculating the error metric

            //if ( (vset[3*i] - newX >= errorMetric) || (vset[3*i+1] - newY >= errorMetric) ||
            //     (vset[3*i+2] - newZ >= errorMetric) ) {
            //    qf->get_attributes().subdivide = true;
            //    break;
            //}
        }

        int i = 0;
        unsigned* fset = new unsigned[ 3*qf->get_attributes().faces.size()];
        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {
            Face* f = *it;
            Vertex* v0 = f->get_halfedge()->get_origin();
            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

            fset[3*i] = v0->get_attributes().tempInt;
            fset[3*i+1] = v1->get_attributes().tempInt;
            fset[3*i+2] = v2->get_attributes().tempInt;
            i++;
        }

        //vtkReader reader;

        //std::string filename = "debug-patch-errormetric";
        //std::stringstream str;
        //str << qf->get_attributes().id;

        //std::string filename1 = filename + std::string("-3D-") + str.str();
        //reader.saveVTKTriPatchWithoutTvalue(filename1, size, newCoordSet, qf->get_attributes().faces.size(), fset);

    }

    // ****************************************
    //            SIMPLE ERROR METRIC
    // ****************************************

    /**
       * \fn bool simple_error_metric(QFace* qf)
       *
       * \brief Returns true if the amount of triangles in a region passes through a certain threshold.
       *
       * \param qf Region of verification.
       *
       */

    bool simple_error_metric(QFace* qf){


        if (qf->get_attributes().faces.size() > _threshold)
            return true;

        return false;

    }

    // ****************************************
    //            TEMPLATE MARK EDGES
    // ****************************************

    //    void templateMarkEdges(QMesh* qm){

    //        QHalfedge* qh;
    //        QFace* qf;
    //        QVertex* qv;
    //        int templateType = 0;

    //        for ( QMesh::QFaceIterator it = qm->faces_begin(); it != qm->faces_end(); it++ ) {

    //            qf = *it;
    //            qh = qf->get_halfedge();
    //            qv = qh->get_origin();

    //            do {

    //                if ( qh->get_edge()->get_attributes().subdivide == true )
    //                    templateType++;

    //                qh = qh->get_next();

    //            } while (qh != qf->get_halfedge() );

    //            qf->get_attributes().orientation =
    //            qf->get_attributes().subdivide = true;

    //            templateType = 0;
    //        }
    //    }

    // ****************************************
    //            TEMPLATE MARK FACES
    // ****************************************

    bool templateMarkFaces(QMesh* qm){

        for ( QMesh::QEdgeIterator it = qm->edges_begin() ; it != qm->edges_end() ; it++ ){

            QEdge* qe = *it;

            qe->get_attributes().checkEdge = false;
            qe->get_attributes().subdivide = false; // verificar

        }

        for ( QMesh::QFaceIterator it = qm->faces_begin() ; it != qm->faces_end() ; it++ ){

            QFace* qf = *it;

            if ( simple_error_metric(qf) ) {

                qf->get_halfedge()->get_edge()->get_attributes().checkEdge = true;
                qf->get_halfedge()->get_next()->get_edge()->get_attributes().checkEdge = true;
                qf->get_halfedge()->get_next()->get_next()->get_edge()->get_attributes().checkEdge = true;
                qf->get_halfedge()->get_prev()->get_edge()->get_attributes().checkEdge = true;

            }

        }

        int count = 0;
        for ( QMesh::QFaceIterator it = qm->faces_begin() ; it != qm->faces_end() ; it++ ){

            QFace* qf = *it;
            templateMarkFace(qf);

            if( qf->get_attributes().orientation != myQFaceAttributes::TEMPLATE_NNNN )
                count++;

        }

        return (count > 0);

    }

    void templateMarkFace(QFace* qf){

        QHalfedge* qh = qf->get_halfedge();

        if (qh->get_edge()->get_attributes().checkEdge == true){

            // Possibilities are: DNNN, DNND, DNDN, DNDD, DDNN, DDND, DDDN, DDDD

            if ( qh->get_next()->get_edge()->get_attributes().checkEdge == true ){

                // Possibilities are: DDNN, DDND, DDDN, DDDD

                if (qh->get_next()->get_next()->get_edge()->get_attributes().checkEdge == true){

                    // Possibilities are: DDDN, DDDD

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true){

                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDDD; //Cross Template;


                    } else {

                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDDN; // Invalid

                    }

                } else {

                    // Possibilities are: DDNN, DDND

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true){
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDND; // Invalid

                    }
                    else {
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDNN; // Left Bottom Subdivision

                    }

                }

            } else {

                // Possibilities; DNNN, DNND, DNDN, DNDD

                if (qh->get_next()->get_next()->get_edge()->get_attributes().checkEdge == true){

                    // Possibilities: DNDN, DNDD

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true)
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNDD; // Invalid
                    else
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNDN; // Left to Right Subdivision

                    // Possibilities: DNND, DNNN

                } else {

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true)
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNND; // Left top Subdivision
                    else
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNNN; // Invalid

                }

            }

        } else {

            // Possibilities: NNNN, NNND, NNDN, NNDD, NDNN, NDND, NDDN, NDDD

            if (qh->get_next()->get_edge()->get_attributes().checkEdge == true){

                // Possibilities: NDNN, NDND, NDDN, NDDD

                if (qh->get_next()->get_next()->get_edge()->get_attributes().checkEdge == true){

                    // Possibilities: NDDN, NDDD

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true)
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDD; // Invalid
                    else
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDN;

                } else {

                    // Possibilities: NDNN, NDND

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true)
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDND;
                    else
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDNN;

                }

            } else {

                // Possibilities: NNNN, NNND, NNDN, NNDD

                if (qh->get_next()->get_next()->get_edge()->get_attributes().checkEdge == true){

                    // Possibilities: NNDN, NNDD

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true)
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDD;
                    else
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDN;

                } else {

                    // Possibilities: NNNN, NNND

                    if (qh->get_prev()->get_edge()->get_attributes().checkEdge == true){
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNND;
                    }
                    else
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNNN;

                }

            }

        }
    }

    // ****************************************************************
    //            ELIMINATE TEMPLATE 3-NEIGHBOURS OR 4-NEIGHBOURS
    // ****************************************************************

    //    QMesh* eliminateTemplate3or4(QMesh* qm){

    //        bool reboot;
    //        QFace* qf;

    //        do {

    //            reboot = false;

    //            for ( QMesh::QFaceIterator it = qm->faces_begin(); it != qm->faces_end(); it++ ) {

    //                qf = *it;

    //                if ( qf->get_attributes().templateType == 4 ||
    //                     qf->get_attributes().templateType == 1) {

    //                    qf->get_attributes().subdivide = true;
    //                    reboot = true;
    //                }

    //            }

    //        } while (reboot);

    //    }

    // ****************************************************************
    //            SETTING QEDGE'S ATTRIBUTE SUBDIVIDE
    // ****************************************************************

    QMesh* settingEdgeSubdivide(QMesh* qm){

        QFace* qf;
        QHalfedge* qh;

        qh = qf->get_halfedge();

        for ( QMesh::QFaceIterator it = qm->faces_begin(); it != qm->faces_end(); it++ ) {

            qf = *it;

            if ( qf->get_attributes().subdivide == true ){

                do {

                    if ( qh->get_mate()->get_face()->get_attributes().subdivide == true )
                        qh->get_edge()->get_attributes().subdivide = true;

                    qh = qh->get_next();

                } while (qh != qf->get_halfedge() );
            }
        }
    }

    // ************************************************************************
    //            DEFINING THE QFACE's ATTRIBUTE BORDER  OF THE REGION D
    // ************************************************************************

    QMesh* definingBBorder(QMesh* qm){

        QEdge* borderEdge;
        QHalfedge* qh;

        // Getting one border edge

        for ( QMesh::QEdgeIterator it = qm->edges_begin(); it != qm->edges_end(); it++ ){

            borderEdge = *it;

            if ( (borderEdge->get_first_halfedge()->get_face()->get_attributes().subdivide == true &&
                  borderEdge->get_first_halfedge()->get_face()->get_attributes().border == false &&
                  borderEdge->get_second_halfedge()->get_face()->get_attributes().subdivide == false &&
                  borderEdge->get_second_halfedge()->get_face()->get_attributes().border == false) )

            {

                qh = borderEdge->get_second_halfedge();

                // Marking all the border faces with face.border = true
                do {

                    qh->get_edge()->get_attributes().subdivide = true;
                    qh->get_face()->get_attributes().border = true;
                    qh = qh->get_prev()->get_mate();

                } while ( qh != borderEdge->get_first_halfedge() );
            }

            if ( (borderEdge->get_second_halfedge()->get_face()->get_attributes().subdivide == true &&
                  borderEdge->get_second_halfedge()->get_face()->get_attributes().border == false &&
                  borderEdge->get_first_halfedge()->get_face()->get_attributes().subdivide == false &&
                  borderEdge->get_first_halfedge()->get_face()->get_attributes().border == false))

            {

                qh = borderEdge->get_first_halfedge();

                // Marking all the border faces with face.border = true
                do {

                    qh->get_edge()->get_attributes().subdivide = true;
                    qh->get_face()->get_attributes().border = true;
                    qh = qh->get_prev()->get_mate();

                } while ( qh != borderEdge->get_second_halfedge() );

            }

        }

    }

    // ****************************************************************
    //                        FIND MIDDLE VERTEX
    // ****************************************************************

    /**
      * \fn std::list< Vertex* >::iterator findMiddleVertex(std::list< Vertex* >& lv)
      *
      * \brief Returns an iterator for the central tri vertex of the lv vertex list.
      *
      * \param lv vertex list made by the geodesical cut, with minimum size >= 3.
      *
      */

    std::list< Vertex* >::iterator findMiddleVertex(std::list< Vertex* >& lv, geodesicMeshing* geo, std::set< pack3<Vertex*,double,double> >* vertices_uv ){

        ASSERT ( lv.size() >= 2 );

        double disttotal = 0, distatual = 0;
        std::list< Vertex* >::iterator it1 = lv.begin();
        std::list< Vertex* >::iterator it2 = it1;
        it2++;

        Vertex* v1 = *it1;
        Vertex* v2 = *it2;

        // Sum of all distances between the pairs of 3D Vertex
        while(it2 != lv.end()) {

            double d = common::Geometric::distance_3d( v1->x(), v1->y(), v1->z(), v2->x(), v2->y(), v2->z() );
            disttotal += d;

            it1++;
            it2++;
        }

        it1 = lv.begin();
        it2 = it1;
        it2++;

        v1 = *it1;
        v2 = *it2;

        // Actually getting the middle vertex
        double disttotal2 = disttotal/2.0;
        double lastedge = common::Geometric::distance_3d( v1->x(), v1->y(), v1->z(), v2->x(), v2->y(), v2->z() );
        while ((distatual+lastedge) < disttotal2) {

            distatual += lastedge;

            it1++;
            it2++;

            lastedge = common::Geometric::distance_3d( v1->x(), v1->y(), v1->z(), v2->x(), v2->y(), v2->z() );
        }


        ASSERT( disttotal2 > distatual );
        ASSERT( (disttotal2 - distatual) < lastedge + C_THRESHOLD );
        double prop = (disttotal2 - distatual)/lastedge;

        if( fabs(prop) < P_THRESHOLD) {
            ASSERT( lv.size() > 2 );
            return it1;
        }
        else if( fabs(prop-1.0) < P_THRESHOLD) {
            ASSERT( lv.size() > 2 );
            return it2;
        }

        Vertex* v = *it1;
        Halfedge* he = v->get_halfedge();

        while ( he->get_next()->get_origin() != *it2 ) {
            he = he->get_prev()->get_mate();
            ASSERT(he != v->get_halfedge());
        }

        ASSERT( he->get_face()->get_attributes().region.size() == 1 );
        ASSERT( he->get_mate()->get_face()->get_attributes().region.size() == 1 );

        QFace* f1 = *(he->get_face()->get_attributes().region.begin());
        QFace* f2 = *(he->get_mate()->get_face()->get_attributes().region.begin());

        Vertex* viz1 = he->get_prev()->get_origin();
        Vertex* viz2 = he->get_mate()->get_prev()->get_origin();


        std::list< geodesics::GP*> pts;
        pts.push_back( new geodesics::GPV(viz1) );
        if( he->get_edge()->get_first_halfedge() == he)
            pts.push_back( new geodesics::GPE(he->get_edge(), prop) );
        else
            pts.push_back( new geodesics::GPE(he->get_edge(), 1.0 - prop) );
        pts.push_back( new geodesics::GPV(viz2) );

        std::list< Vertex* > vtxs;
        geo ->remeshingGeodesics(pts, vtxs);

        for( std::list< geodesics::GP*>::iterator it = pts.begin(); it != pts.end(); ++it)
            delete *it;
        pts.clear();

        ASSERT( vtxs.size() == 1 );

        std::list< Vertex* >::iterator itnew = lv.insert(it2, *(vtxs.begin()));

        Vertex* vnew = *itnew;
        Halfedge *he_new = vnew->get_halfedge();
        do {

            if( he_new->get_next()->get_origin() == *it1 ) {
                he_new->get_face()->get_attributes().region.insert(f2);
                f2->get_attributes().faces.insert( he_new->get_face() );
            } else if( he_new->get_next()->get_origin() == viz1 ) {
                he_new->get_face()->get_attributes().region.insert(f1);
                f1->get_attributes().faces.insert( he_new->get_face() );
            } else if( he_new->get_next()->get_origin() == *it2 ) {
                he_new->get_face()->get_attributes().region.insert(f1);
                f1->get_attributes().faces.insert( he_new->get_face() );
            } else if( he_new->get_next()->get_origin() == viz2 ) {
                he_new->get_face()->get_attributes().region.insert(f2);
                f2->get_attributes().faces.insert( he_new->get_face() );
            } else {
                ASSERT( false );
            }
            ASSERT( he_new->get_face()->get_attributes().region.size() == 1);
            he_new = he_new->get_prev()->get_mate();
        } while( he_new != vnew->get_halfedge());

        if( vertices_uv ) {
            std::set< pack3<Vertex*,double,double> >::iterator it1uv = vertices_uv->find( pack3<Vertex*,double,double>(*it1,0,0) );
            ASSERT( it1uv != vertices_uv->end() );
            std::set< pack3<Vertex*,double,double> >::iterator it2uv = vertices_uv->find( pack3<Vertex*,double,double>(*it2,0,0) );
            ASSERT( it2uv != vertices_uv->end() );

            vertices_uv->insert(pack3<Vertex*,double,double>(vnew,
                                                             (it1uv->second + it2uv->second)/2.0,
                                                             (it1uv->third + it2uv->third)/2.0 ) ); // deveria ser prop ou 1.0-prop
        }



        return itnew;
    }

    // ****************************************************************
    //                         ADD VERTEX
    // ****************************************************************

    /**
      * \fn QVertex* add_vertex(double x, double y, double z, QMesh* qm)
      *
      * \brief Creates a QVertex, inserts in the QMesh and returns a pointer for itself.
      *
      * \param x First Cartesian coordinate of the recently created QVertex.
      *
      * \param y Second Cartesian coordinate of the recently created QVertex.
      *
      * \param z Third Cartesian coordinate of the recently created QVertex.
      *
      * \param qm Pointer to a QMesh whose QVertex will belong.
      */

    QVertex* add_vertex(double x, double y, double z, QMesh* qm) {

        QVertex* vertex = new QVertex(x, y, z, 0);

        qm->add_vertex(vertex);

        return vertex;
    }

    /**
      * \fn void checkNeighbours(QVertex **newVertex, std::list< Vertex* >::iterator vertexIterators[], bool vertexValidation[], QHalfedge *he, QMesh* qm)
      *
      * \brief Creates or sets the vertex for subdivision. It also gives the vertexIterators of the new vertex and a
      *        flag vertexValidation to know if the neighbour is subdivided
      *
      * \param newVertex Array of the new QVertex created.
      *
      * \param vertexIterators Iterators for getting the position of the middle vertex in the given lv list of vertex.
      *
      * \param vertexValidation Booleans to know which neighbour is already subdivided.
      *
      * \param qe Half-edge of the left edge of the quad.
      *
      * \param qm Pointer to a QMesh whose QVertex will belong.
      */

    void checkNeighbours(std::vector<QVertex*> *newVertex, std::list< Vertex* >::iterator vertexIterators[], std::vector<bool> *vertexValidation, QHalfedge *he, QMesh* qm, geodesicMeshing* geo){

        std::vector<Vertex*> vertexTri;
        QHalfedge* heVerification = he;

        /**
           Counter-clockwise order, beginning from left edge
           (Left - Bottom - Right - Top)
          */

        for ( int i = 0 ; i < 4 ; i++ ){

            if (he->get_edge()->get_attributes().subdivide == true){

                newVertex->push_back(he->get_mate()->get_next()->get_origin());
                vertexValidation->push_back(true);

            } else {

                vertexIterators[i] = findMiddleVertex(he->get_edge()->get_attributes().geoline, geo, 0);

                vertexValidation->push_back(false);

                vertexTri.push_back(*vertexIterators[i]);

                newVertex->push_back(add_vertex(vertexTri.back()->x(),
                                                vertexTri.back()->y(),
                                                vertexTri.back()->z(), qm));

                newVertex->back()->get_attributes().setGP(new geodesics::GPV(vertexTri.back()));
            }

            he = he->get_next();

        }

        ASSERT( heVerification == he );

    }

    /**
      * \fn void checkNeighboursTemplate1(QVertex **newVertex, std::list< Vertex* >::iterator vertexIterators[], bool vertexValidation[], QHalfedge *he, QMesh* qm)
      *
      * \brief Creates or sets the vertex for subdivision. It also gives the vertexIterators of the new vertex and a
      *        flag vertexValidation to know if the neighbour is subdivided
      *
      * \param newVertex Array of the new QVertex created.
      *
      * \param vertexIterators Iterators for getting the position of the middle vertex in the given lv list of vertex.
      *
      * \param vertexValidation Booleans to know which neighbour is already subdivided.
      *
      * \param qe Half-edge of the left edge of the quad.
      *
      * \param qm Pointer to a QMesh whose QVertex will belong.
      */

    void checkNeighboursTemplate1(std::vector<QVertex*> &newVertex, std::list< Vertex* >::iterator vertexIterators[], std::vector<bool>& vertexValidation, QHalfedge *he, QMesh* qm, geodesicMeshing* geo){

        std::vector<Vertex*> vertexTri;
        QHalfedge* heVerification = he;

        /**
           Counter-clockwise order, beginning from left edge
           (Left - Bottom - Right - Top)
          */



        for ( int i = 0 ; i < 4 ; i ++ ){

            if ( i % 2 == 1 ){

                if (he->get_edge()->get_attributes().subdivide == true){

                    newVertex.push_back(he->get_mate()->get_next()->get_origin());
                    vertexValidation.push_back(true);

                } else {

                    vertexIterators[i] = findMiddleVertex(he->get_edge()->get_attributes().geoline, geo, 0);

                    vertexValidation.push_back(false);

                    vertexTri.push_back(*vertexIterators[i]);

                    newVertex.push_back(add_vertex(vertexTri.back()->x(),
                                                   vertexTri.back()->y(),
                                                   vertexTri.back()->z(), qm));

                    newVertex.back()->get_attributes().setGP(new geodesics::GPV(vertexTri.back()));
                }
            }

            he = he->get_next();

        }

        ASSERT( heVerification == he );

    }

    /**
      * \fn void checkNeighboursTemplate3(QVertex **newVertex, std::list< Vertex* >::iterator vertexIterators[], bool vertexValidation[], QHalfedge *he, QMesh* qm)
      *
      * \brief Creates or sets the vertex for subdivision. It also gives the vertexIterators of the new vertex and a
      *        flag vertexValidation to know if the neighbour is subdivided
      *
      * \param newVertex Array of the new QVertex created.
      *
      * \param vertexIterators Iterators for getting the position of the middle vertex in the given lv list of vertex.
      *
      * \param vertexValidation Booleans to know which neighbour is already subdivided.
      *
      * \param qe Half-edge of the left edge of the quad.
      *
      * \param qm Pointer to a QMesh whose QVertex will belong.
      */

    void checkNeighboursTemplate3(std::vector<QVertex*> &newVertex, std::list< Vertex* >::iterator vertexIterators[], std::vector<bool> &vertexValidation, QHalfedge *he, QMesh* qm, geodesicMeshing *geo){

        std::vector<Vertex*> vertexTri;
        QHalfedge* heVerification = he;

        // Left Edge

        if (he->get_edge()->get_attributes().subdivide == true){

            newVertex.push_back(he->get_mate()->get_next()->get_origin());
            vertexValidation.push_back(true);

        } else {

            vertexIterators[0] = findMiddleVertex(he->get_edge()->get_attributes().geoline, geo, 0);

            vertexValidation.push_back(false);

            vertexTri.push_back(*vertexIterators[0]);

            newVertex.push_back(add_vertex(vertexTri.back()->x(),
                                           vertexTri.back()->y(),
                                           vertexTri.back()->z(), qm));

            newVertex.back()->get_attributes().setGP(new geodesics::GPV(vertexTri.back()));
        }

        he = he->get_prev();

        // Top Edge

        if (he->get_edge()->get_attributes().subdivide == true){

            newVertex.push_back(he->get_mate()->get_next()->get_origin());
            vertexValidation.push_back(true);

        } else {

            vertexIterators[3] = findMiddleVertex(he->get_edge()->get_attributes().geoline, geo, 0);

            vertexValidation.push_back(false);

            vertexTri.push_back(*vertexIterators[3]);

            newVertex.push_back(add_vertex(vertexTri.back()->x(),
                                           vertexTri.back()->y(),
                                           vertexTri.back()->z(), qm));

            newVertex.back()->get_attributes().setGP(new geodesics::GPV(vertexTri.back()));
        }

        he = he->get_next();

        ASSERT( heVerification == he );
    }

    // ****************************************************************
    //                        FIXING COLOR
    // ****************************************************************

    /**
      * \fn void fixingColor(std::list< Vertex* >& lv, QFace* qf)
      *
      * \brief Sets the color of the new triangles created by the geodesic cut.
      *
      * \param lv List of vertex from the geodesic cut.
      *
      * \param qf QFace (region) which the new triangles will be inserted.
      */

    void fixingColor(std::list< Vertex* >& lv, QFace* qf){

        if( lv.size() == 2 )
        {
            std::list< Vertex* >::iterator iti = lv.begin();
            std::list< Vertex* >::iterator ite = lv.begin();
            ite++;

            Vertex* v = *iti;
            Halfedge* he = v->get_halfedge();

            while ( he->get_next()->get_origin() != *ite ) {
                he = he->get_prev()->get_mate();
                ASSERT(he != v->get_halfedge());
            }

            he->get_face()->get_attributes().region.clear();
            he->get_face()->get_attributes().region.insert(qf);

            return;
        }

        ASSERT( lv.size() > 2 );
        std::list< Vertex* >::iterator it = lv.begin();
        it++;
        std::list< Vertex* >::iterator itprox = it;
        itprox++;

        while ( itprox != lv.end() )
        {
            Vertex* v = *it;
            Halfedge* he = v->get_halfedge();

            do {

                he->get_face()->get_attributes().region.clear();
                he->get_face()->get_attributes().region.insert(qf);

                he = he->get_mate()->get_next();

            } while ( he != v->get_halfedge() );

            it++;
            itprox++;
        }
    }

    // ****************************************************************
    //                 CORRECT ORIENTATION for Edges
    // ****************************************************************

    /**
      * \fn void correctOrientation( QEdge* qe, std::list<Vertex*>& lv )
      *
      * \brief Corrects the geoline orientation based on the QEdge.
      *
      * \param qe QEdge which first Half-edge must be the guide of the geoline's direction.
      *
      * \param lv List of vertex from the geodesic cut.
      */

    void correctOrientation( QEdge* qe, std::list<Vertex*>& lv ) {

        ASSERT( qe->get_first_halfedge() );
        correctOrientation( qe->get_first_halfedge(), lv);
    }

    // ****************************************************************
    //                 CORRECT ORIENTATION for Half-Edges
    // ****************************************************************

    /*
      Corrige a orientação da geoline de acordo com a half-edge.

      Parâmetros:
            qhe       = QHalfedge cuja direção a geoline deve se basear.
            lv        = lista de vértices provenientes do corte da geodésica.
    */

    /**
      * \fn void correctOrientation( QHalfedge* qhe, std::list<Vertex*>& lv )
      *
      * \brief Corrects the geoline orientation based on the QHalfedge.
      *
      * \param qhe QHalfedge which will be the guide of the geoline's direction.
      *
      * \param lv List of vertex from the geodesic cut.
      */

    void correctOrientation( QHalfedge* qhe, std::list<Vertex*>& lv ) {

        //        if( (debug::iteration==4)&&(debug::templatenumber==7)){
        //            int a=10;
        //        }


        Vertex* helv1 = 0;
        Vertex* helv2 = 0;


        Vertex* vi = lv.front();
        Vertex* ve = lv.back();

        QVertex* qvfirst, *qvsecond;
        qvfirst = qhe->get_origin();
        qvsecond = qhe->get_next()->get_origin();
        Vertex* first = dynamic_cast<geodesics::GPV*>(qvfirst->get_attributes().getGP())->get_vertex();
        Vertex* second = dynamic_cast<geodesics::GPV*>(qvsecond->get_attributes().getGP())->get_vertex();


        if( qhe->get_edge()->get_second_halfedge() ) {

            QHalfedge* qhem = qhe->get_mate();

            QVertex* qvmfirst, *qvmsecond;
            qvmfirst = qhem->get_origin();
            qvmsecond = qhem->get_next()->get_origin();
            Vertex* mfirst = dynamic_cast<geodesics::GPV*>(qvmfirst->get_attributes().getGP())->get_vertex();
            Vertex* msecond = dynamic_cast<geodesics::GPV*>(qvmsecond->get_attributes().getGP())->get_vertex();


            if( qvfirst == qvmsecond ) {
                if(qvsecond == qvmfirst ) {
                    if( (vi == first )&&(ve == second) ) {
                        // orientacao correta
                        ASSERT( vi == msecond );
                        ASSERT( ve == mfirst );
                    } else if ((vi == second )&&(ve == first)) {
                        ASSERT( ve == msecond );
                        ASSERT( vi == mfirst );
                        lv.reverse();
                    } else {
                        ASSERT(false);
                    }
                } else {
                    ASSERT( qhe->get_next()->get_mate()->get_next()->get_next()->get_origin() == qvmfirst );
                    if( vi == first ) {
                    } else if( ve == second ) {
                        lv.reverse();
                    } else {
                        ASSERT( false );
                    }
                }
                helv1 = first;
                helv2 = second;
            } else {
                if(qvsecond == qvmfirst ) {
                    ASSERT( qhem->get_next()->get_mate()->get_next()->get_next()->get_origin() == qvfirst );
                    if( vi == second ) {
                        lv.reverse();
                        helv1 = mfirst;
                        helv2 = msecond;
                    } else if( ve == second ) {
                        helv1 = mfirst;
                        helv2 = msecond;
                    } else if( vi == first ) {
                        QHalfedge* helv = qhem->get_next()->get_mate()->get_next();
                        helv1 = dynamic_cast<geodesics::GPV*>(helv->get_origin()->get_attributes().getGP())->get_vertex();
                        helv2 = dynamic_cast<geodesics::GPV*>(helv->get_next()->get_origin()->get_attributes().getGP())->get_vertex();
                    } else if( ve == first ) {
                        lv.reverse();
                        QHalfedge* helv = qhem->get_next()->get_mate()->get_next();
                        helv1 = dynamic_cast<geodesics::GPV*>(helv->get_origin()->get_attributes().getGP())->get_vertex();
                        helv2 = dynamic_cast<geodesics::GPV*>(helv->get_next()->get_origin()->get_attributes().getGP())->get_vertex();
                    } else {
                        ASSERT( false );
                    }
                } else {
                    ASSERT( false );
                }
            }
        } else {
            if( (vi == first )&&(ve == second) ) {
            } else if ((vi == second )&&(ve == first)) {
                lv.reverse();
            } else {
                ASSERT(false);
            }
            helv1 = first;
            helv2 = second;
        }

        ASSERT( (lv.front() == helv1)||(lv.front() == helv2) );
        ASSERT( (lv.back() == helv1)||(lv.back() == helv2) );
        ASSERT( lv.front() != lv.back() );
        ASSERT( helv1 != helv2 );
    }

    // ****************************************************************
    //                   CORRECT LIST
    // ****************************************************************

    /**
      * \fn void correctList(QMesh* qm, std::list< QFace* >* _list_of_qfaces)
      *
      * \brief Eliminates invalid templates.
      *
      * \param qm Pointer to a quadrilateral mesh.
      *
      * \param _list_of_qfaces First QFaces* of the first main quad.
      */

    void correctList(QMesh *qm){

        std::multiset< pack3<double, QFace*, myQFaceAttributes::OrientationCut> > setFaces;

        for( QMesh::QFaceIterator it =qm->faces_begin(); it != qm->faces_end(); it++ ){

            QFace* qf = *it;

            if (qf->get_attributes().getTemplateType() == 1)
                setFaces.insert( pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qf,qf->get_attributes().orientation) );
            else if (qf->get_attributes().getTemplateType() == 5)
                setFaces.insert( pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qf,qf->get_attributes().orientation) );

        }

        while (!setFaces.empty()){

            std::multiset< pack3<double, QFace*, myQFaceAttributes::OrientationCut> >::iterator it;
            it = setFaces.begin();
            pack3<double, QFace*, myQFaceAttributes::OrientationCut> p = *it;
            setFaces.erase(it);

            QFace* qf = p.second;

            if (qf->get_attributes().orientation == p.third){

                if (qf->get_attributes().getTemplateType() == 5){

                    qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDDD;

                    QHalfedge* qh;
                    qh = qf->get_halfedge();

                    for (int i = 0; i < 4; i++){

                        if (qh->get_edge()->get_attributes().checkEdge == false){

                            qh->get_edge()->get_attributes().checkEdge = true;
                            templateMarkFace(qh->get_mate()->get_face());

                            if (qh->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                                setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qh->get_mate()->get_face(),qh->get_mate()->get_face()->get_attributes().orientation));
                            else if (qh->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                                setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qh->get_mate()->get_face(),qh->get_mate()->get_face()->get_attributes().orientation));
                        }

                        qh = qh->get_next();

                    }

                } else if (qf->get_attributes().getTemplateType() == 1){

                    QHalfedge* qhn = 0;

                    switch(qf->get_attributes().orientation){

                    case myQFaceAttributes::TEMPLATE_DNNN:
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDNN;
                        qhn = qf->get_halfedge()->get_next();
                        break;

                    case myQFaceAttributes::TEMPLATE_NDNN:
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDN;
                        qhn = qf->get_halfedge()->get_next()->get_next();
                        break;

                    case myQFaceAttributes::TEMPLATE_NNDN:
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDD;
                        qhn = qf->get_halfedge()->get_prev();
                        break;

                    case myQFaceAttributes::TEMPLATE_NNND:
                        qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNND;
                        qhn = qf->get_halfedge();
                        break;
                    }

                    ASSERT( qhn );

                    qhn->get_edge()->get_attributes().checkEdge = true;
                    templateMarkFace(qhn->get_mate()->get_face());

                    if (qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                        setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qhn->get_mate()->get_face(),qhn->get_mate()->get_face()->get_attributes().orientation));

                    else if (qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                        setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qhn->get_mate()->get_face(),qhn->get_mate()->get_face()->get_attributes().orientation));

                }
            }
        }
    }


    void correctList2(QMesh *qm){

        std::multiset< pack3<double, QFace*, myQFaceAttributes::OrientationCut> > setFaces;

        for( QMesh::QFaceIterator it =qm->faces_begin(); it != qm->faces_end(); it++ ){

            QFace* qf = *it;

            if (qf->get_attributes().getTemplateType() == 1)
                setFaces.insert( pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qf,qf->get_attributes().orientation) );
            else if (qf->get_attributes().getTemplateType() == 5)
                setFaces.insert( pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qf,qf->get_attributes().orientation) );

        }

        while (!setFaces.empty()){

            std::multiset< pack3<double, QFace*, myQFaceAttributes::OrientationCut> >::iterator it;
            it = setFaces.begin();
            pack3<double, QFace*, myQFaceAttributes::OrientationCut> p = *it;
            setFaces.erase(it);

            QFace* qf = p.second;

            if (qf->get_attributes().orientation == p.third){

                if (qf->get_attributes().getTemplateType() == 5){

                    qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDDD;

                    QHalfedge* qh;
                    qh = qf->get_halfedge();

                    for (int i = 0; i < 4; i++){

                        if (qh->get_edge()->get_attributes().checkEdge == false){

                            qh->get_edge()->get_attributes().checkEdge = true;
                            templateMarkFace(qh->get_mate()->get_face());

                            if (qh->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                                setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qh->get_mate()->get_face(),qh->get_mate()->get_face()->get_attributes().orientation));
                            else if (qh->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                                setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qh->get_mate()->get_face(),qh->get_mate()->get_face()->get_attributes().orientation));
                        }

                        qh = qh->get_next();

                    }

                } else if (qf->get_attributes().getTemplateType() == 1){

                    QHalfedge* qhn = 0;
                    QHalfedge* qhp = 0;

                    switch(qf->get_attributes().orientation){

                    case myQFaceAttributes::TEMPLATE_DNNN:
                        qhn = qf->get_halfedge()->get_next();
                        qhp = qf->get_halfedge()->get_prev();
                        break;

                    case myQFaceAttributes::TEMPLATE_NDNN:
                        qhn = qf->get_halfedge()->get_next()->get_next();
                        qhp = qf->get_halfedge();
                        break;

                    case myQFaceAttributes::TEMPLATE_NNDN:
                        qhn = qf->get_halfedge()->get_prev();
                        qhp = qf->get_halfedge()->get_next();
                        break;

                    case myQFaceAttributes::TEMPLATE_NNND:
                        qhn = qf->get_halfedge();
                        qhp = qf->get_halfedge()->get_next()->get_next();
                        break;
                    }

                    ASSERT( qhn );
                    ASSERT( qhp );

                    if((qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 1)||(qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 5)) {

                        switch(qf->get_attributes().orientation){

                        case myQFaceAttributes::TEMPLATE_DNNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDNN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NDNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNDN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDD;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNND:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNND;
                            break;
                        }

                        qhn->get_edge()->get_attributes().checkEdge = true;
                        templateMarkFace(qhn->get_mate()->get_face());

                        if (qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qhn->get_mate()->get_face(),qhn->get_mate()->get_face()->get_attributes().orientation));

                        else if (qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qhn->get_mate()->get_face(),qhn->get_mate()->get_face()->get_attributes().orientation));

                    } else {

                        switch(qf->get_attributes().orientation){

                        case myQFaceAttributes::TEMPLATE_DNNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNND;
                            break;

                        case myQFaceAttributes::TEMPLATE_NDNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDNN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNDN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNND:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDD;
                            break;
                        }

                        qhp->get_edge()->get_attributes().checkEdge = true;
                        templateMarkFace(qhp->get_mate()->get_face());

                        if (qhp->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qhp->get_mate()->get_face(),qhp->get_mate()->get_face()->get_attributes().orientation));

                        else if (qhp->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qhp->get_mate()->get_face(),qhp->get_mate()->get_face()->get_attributes().orientation));
                    }
                }
            }
        }
    }


    void correctList3(QMesh *qm){

        std::multiset< pack3<double, QFace*, myQFaceAttributes::OrientationCut> > setFaces;

        for( QMesh::QFaceIterator it =qm->faces_begin(); it != qm->faces_end(); it++ ){

            QFace* qf = *it;

            if (qf->get_attributes().getTemplateType() == 1)
                setFaces.insert( pack3<double, QFace*, myQFaceAttributes::OrientationCut>(1,qf,qf->get_attributes().orientation) );
            else if (qf->get_attributes().getTemplateType() == 5)
                setFaces.insert( pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qf,qf->get_attributes().orientation) );

        }

        while (!setFaces.empty()){

            std::multiset< pack3<double, QFace*, myQFaceAttributes::OrientationCut> >::iterator it;
            it = setFaces.begin();
            pack3<double, QFace*, myQFaceAttributes::OrientationCut> p = *it;
            setFaces.erase(it);

            QFace* qf = p.second;

            if (qf->get_attributes().orientation == p.third){

                if (qf->get_attributes().getTemplateType() == 5){

                    qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDDD;

                    QHalfedge* qh;
                    qh = qf->get_halfedge();

                    for (int i = 0; i < 4; i++){

                        if (qh->get_edge()->get_attributes().checkEdge == false){

                            qh->get_edge()->get_attributes().checkEdge = true;
                            templateMarkFace(qh->get_mate()->get_face());

                            if (qh->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                                setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0.4,qh->get_mate()->get_face(),qh->get_mate()->get_face()->get_attributes().orientation));
                            else if (qh->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                                setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qh->get_mate()->get_face(),qh->get_mate()->get_face()->get_attributes().orientation));
                        }

                        qh = qh->get_next();

                    }

                } else if (qf->get_attributes().getTemplateType() == 1){

                    QHalfedge* qhn = 0;
                    QHalfedge* qhp = 0;

                    switch(qf->get_attributes().orientation){

                    case myQFaceAttributes::TEMPLATE_DNNN:
                        qhn = qf->get_halfedge()->get_next();
                        qhp = qf->get_halfedge()->get_prev();
                        break;

                    case myQFaceAttributes::TEMPLATE_NDNN:
                        qhn = qf->get_halfedge()->get_next()->get_next();
                        qhp = qf->get_halfedge();
                        break;

                    case myQFaceAttributes::TEMPLATE_NNDN:
                        qhn = qf->get_halfedge()->get_prev();
                        qhp = qf->get_halfedge()->get_next();
                        break;

                    case myQFaceAttributes::TEMPLATE_NNND:
                        qhn = qf->get_halfedge();
                        qhp = qf->get_halfedge()->get_next()->get_next();
                        break;
                    }

                    ASSERT( qhn );
                    ASSERT( qhp );

                    if((qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 1)||(qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 5)) {

                        switch(qf->get_attributes().orientation){

                        case myQFaceAttributes::TEMPLATE_DNNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDNN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NDNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNDN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDD;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNND:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNND;
                            break;
                        }

                        qhn->get_edge()->get_attributes().checkEdge = true;
                        templateMarkFace(qhn->get_mate()->get_face());

                        if (qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0.5,qhn->get_mate()->get_face(),qhn->get_mate()->get_face()->get_attributes().orientation));

                        else if (qhn->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qhn->get_mate()->get_face(),qhn->get_mate()->get_face()->get_attributes().orientation));

                    } else {

                        switch(qf->get_attributes().orientation){

                        case myQFaceAttributes::TEMPLATE_DNNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DNND;
                            break;

                        case myQFaceAttributes::TEMPLATE_NDNN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_DDNN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNDN:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NDDN;
                            break;

                        case myQFaceAttributes::TEMPLATE_NNND:
                            qf->get_attributes().orientation = myQFaceAttributes::TEMPLATE_NNDD;
                            break;
                        }

                        qhp->get_edge()->get_attributes().checkEdge = true;
                        templateMarkFace(qhp->get_mate()->get_face());

                        if (qhp->get_mate()->get_face()->get_attributes().getTemplateType() == 1)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0.5,qhp->get_mate()->get_face(),qhp->get_mate()->get_face()->get_attributes().orientation));

                        else if (qhp->get_mate()->get_face()->get_attributes().getTemplateType() == 5)
                            setFaces.insert(pack3<double, QFace*, myQFaceAttributes::OrientationCut>(0,qhp->get_mate()->get_face(),qhp->get_mate()->get_face()->get_attributes().orientation));
                    }
                }
            }
        }
    }

    // ****************************************************************
    //                   SUBDIVIDE
    // ****************************************************************

    /**
      * \fn void subdivide(TMesh* tm, QMesh* qm, geodesicMeshing* geo, std::list< QFace* > _list_of_qfaces)
      *
      * \brief Verifies a list of QFaces to see if there are any .
      *
      * \param faces_set Set of faces (attribute from QFace).
      *
      * \param he Base half-edge for the loop.
      */

    void subdivide(TMesh* tm, QMesh* qm, geodesicMeshing* geo ){

        debug::templatenumber = 0;

        for( QMesh::QFaceIterator it = qm->faces_begin(); it != qm->faces_end(); it++ ){

            QFace* qf = *it;

            switch (qf->get_attributes().getTemplateType()) {

            case 0:
                ASSERT(qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_NNNN);
                break;
            case 2:
                if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_NDND)
                    makingTestsTemplate1(tm, qm, qf->get_halfedge(),geo);
                else if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_DNDN)
                    makingTestsTemplate1(tm, qm, qf->get_halfedge()->get_next(), geo);
                else {
                    ASSERT(false);
                }
                break;
            case 3:
                if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_DDNN)
                    makingTestsTemplate3(tm,qm,qf->get_halfedge()->get_next(),geo);
                else if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_DNND)
                    makingTestsTemplate3(tm,qm,qf->get_halfedge(),geo);
                else if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_NDDN)
                    makingTestsTemplate3(tm,qm,qf->get_halfedge()->get_next()->get_next(),geo);
                else if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_NNDD)
                    makingTestsTemplate3(tm,qm,qf->get_halfedge()->get_prev(),geo);
                else {
                    ASSERT(false);
                }
                break;
            case 4:
                if (qf->get_attributes().orientation == myQFaceAttributes::TEMPLATE_DDDD)
                    makingTestsTemplate2(tm, qm, qf->get_halfedge(), geo);
                else {
                    ASSERT(false);
                }
                break;
            default:
                ASSERT(false);
            }

            debug::templatenumber++;

            if( debug::templatenumber==1530 && debug::iteration==6 )
                vtkWriter::saveVTKTri("debug-it6",tm,qm);

        }

    }

    // ****************************************************************
    //                   VERTEX STAR FACE INSERTION
    // ****************************************************************

    /**
      * \fn void vertexStarFaceInsertion( std::set< Face* > &faces_set, Halfedge* he )
      *
      * \brief Loop which inserts the faces of the vertex star.
      *
      * \param faces_set Set of faces (attribute from QFace).
      *
      * \param he Base half-edge for the loop.
      */

    void vertexStarFaceInsertion( std::set< Face* > &faces_set, Halfedge* he ){

        Vertex* v = he->get_origin();

        do {

            faces_set.insert( he->get_face() );

            he = he->get_prev()->get_mate();
            ASSERT( he->get_origin() == v );
        } while( he!= v->get_halfedge());

    }

    // ****************************************************************
    //                          ADD NEW FACES
    // ****************************************************************

    /**
      * \fn void addNewFaces( std::set< Face* > &faces_set, std::list< Vertex* >& geoline )
      *
      * \brief Goes through the geoline inserting all the vertex star's faces. This must occour
      *        in each geodesic cut, because new faces will be inserted in the mesh.
      *
      * \param faces_set Set of faces (attribute from QFace).
      *
      * \param geoline List of vertex from the geodesic cut.
      */

    void addNewFaces( std::set< Face* > &faces_set, std::list< Vertex* >& geoline ) {

        if( geoline.size() == 2 )
            return;

        ASSERT(geoline.size()>2);
        std::list<Vertex*>::iterator it1 = geoline.begin();
        std::list<Vertex*>::iterator it2 = it1;
        it2++;

        it1++; it2++;
        while( it2 != geoline.end() ) {

            Halfedge* he = (*it1)->get_halfedge();

            vertexStarFaceInsertion(faces_set, he);

            it1++; it2++;
        }
    }

    // ****************************************************************
    //                  ADD NEW FACES (Middle Vertex)
    // ****************************************************************

    /**
      * \fn void addNewFaces( std::set< Face* > &faces_set, QVertex *qv )
      *
      * \brief Given the middle vertex, the faces of it's star are inserted.
      *
      * \param faces_set Set of faces (attribute from QFace).
      *
      * \param qv Middle vertex of the geoline (lv).
      */

    /*
      Dado o vértice qv central, a estrela do mesmo é percorrida e as faces
      que correspondem ao mesmo são inseridas.

      Parâmetros:
            faces_set = conjunto de faces de uma QFace.
            qv        = vértice central.
    */

    void addNewFaces( std::set< Face* > &faces_set, QVertex *qv ) {

        Vertex* v = dynamic_cast<geodesics::GPV*>(qv->get_attributes().getGP())->get_vertex();

        Halfedge* he = v->get_halfedge();

        vertexStarFaceInsertion(faces_set, he);

    }

    // ****************************************************************
    //                        MAKING TESTS
    // ****************************************************************

    /**
      * \fn void makingTestsTemplate2(TMesh* tm, QMesh* qm, geodesicMeshing* geo)
      *
      * \brief "Main" file of a test case using the template 2 ( D - D - D - D ).
      *
      * \param tm Pointer to the triangular mesh.
      *
      * \param qm Pointer to the quadrilateral mesh.
      *
      * \param geo Pointer to the geodesic object.
      */

    void makingTestsTemplate2(TMesh* tm, QMesh* qm, QHalfedge* qh, geodesicMeshing* geo){

        // Template 2

        QFace* qf;
        std::vector<QVertex*> middleVertex;
        QVertex* vertexQuad;
        QOperators qop(qm);
        std::vector<bool> vertexValidation;
        std::list< Vertex* >::iterator vertexIterators[4];

        Vertex* middleTriVertex;
        std::list< Vertex* > lvTri;


        qf = qh->get_face();

        //faceQFacePointerVerification(tm, qm);

        middleVertex.clear();
        vertexValidation.clear();

        checkNeighbours(&middleVertex, vertexIterators, &vertexValidation, qh, qm, geo);

        // Get geodesic from middleVertex[0] to middleVertex[2]
        Vertex* vi = dynamic_cast<geodesics::GPV* >( middleVertex.at(0)->get_attributes().getGP())->get_vertex();
        Vertex* ve = dynamic_cast<geodesics::GPV* >( middleVertex.at(2)->get_attributes().getGP())->get_vertex();
        lvTri.clear();

        std::set< pack3<Vertex*,double,double> > vertices_uv;
        floater(qf,vertices_uv);
        //parametrization( qf, vertices_uv );

        getGeodesics(qf, vi,ve, lvTri, vertices_uv, geo);

        lvTri.push_front(vi);
        lvTri.push_back(ve);

        fixingColor(lvTri, qf);
        addNewFaces( qf->get_attributes().faces, lvTri);

        // Creating QVertex with the middle Vertex coordinates
        std::list< Vertex* >::iterator itMiddle;

        itMiddle = findMiddleVertex(lvTri, geo, &vertices_uv);
        middleTriVertex = *itMiddle;


        vertexQuad = add_vertex(middleTriVertex->x(),
                                middleTriVertex->y(),
                                middleTriVertex->z(), qm);

        vertexQuad->get_attributes().setGP(new geodesics::GPV(middleTriVertex));

        // Creating the geolines of the new Edges
        std::list< Vertex* > lvMiddleLeft;
        lvMiddleLeft.splice(lvMiddleLeft.begin(),lvTri,lvTri.begin(),itMiddle);

        lvMiddleLeft.push_back(*itMiddle);

        std::list< Vertex* > lvMiddleRight;
        lvMiddleRight.splice(lvMiddleRight.begin(),lvTri,itMiddle,lvTri.end());

        std::list< Vertex* > lvMiddleTop;

        //geo->getGeodesic(qf,middleTriVertex,dynamic_cast<geodesics::GPV* >( middleVertex[3]->get_attributes().getGP())->get_vertex(),lvMiddleTop);
        getGeodesics(qf, middleTriVertex,dynamic_cast<geodesics::GPV* >( middleVertex[3]->get_attributes().getGP())->get_vertex(),lvMiddleTop, vertices_uv, geo);

        lvMiddleTop.push_front(middleTriVertex);
        lvMiddleTop.push_back(dynamic_cast<geodesics::GPV* >( middleVertex[3]->get_attributes().getGP())->get_vertex());
        fixingColor(lvMiddleTop, qf);
        addNewFaces( qf->get_attributes().faces, lvMiddleTop);

        std::list< Vertex* > lvMiddleBottom;

        //geo->getGeodesic(qf,middleTriVertex,dynamic_cast<geodesics::GPV* >( middleVertex[1]->get_attributes().getGP())->get_vertex(),lvMiddleBottom);
        getGeodesics(qf, middleTriVertex,dynamic_cast<geodesics::GPV* >( middleVertex[1]->get_attributes().getGP())->get_vertex(),lvMiddleBottom, vertices_uv, geo);

        lvMiddleBottom.push_front(middleTriVertex);
        lvMiddleBottom.push_back(dynamic_cast<geodesics::GPV* >( middleVertex[1]->get_attributes().getGP())->get_vertex());
        fixingColor(lvMiddleBottom, qf);
        addNewFaces( qf->get_attributes().faces, lvMiddleBottom);

        // Splitting the existing geolines of the old Edges

        std::list< Vertex* > lvCCDirection[8];
        int lvcont = 0;

        for ( int i = 0 ; i < 4 ; i++ ){

            if ( vertexValidation[i] == false ){

                std::list< Vertex* >&lv = qh->get_edge()->get_attributes().geoline;
                ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_first_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.front() );
                ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_second_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.back() );

                lvCCDirection[lvcont].splice(lvCCDirection[lvcont].begin(),lv, lv.begin(),vertexIterators[i]);

                lvCCDirection[lvcont].push_back(*vertexIterators[i]);

                lvCCDirection[lvcont+1].splice(lvCCDirection[lvcont+1].begin(),lv, vertexIterators[i], lv.end());

                ASSERT( lv.empty() );


                if( qh->get_edge()->get_second_halfedge() == qh )
                {
                    lvCCDirection[lvcont].swap( lvCCDirection[lvcont+1] );
                    lvCCDirection[lvcont].reverse();
                    lvCCDirection[lvcont+1].reverse();
                }

                qh = qh->get_next();

            } else {

                ASSERT( qh->get_mate()->get_next()->get_mate()->get_next()->get_next()->get_origin() == qh->get_origin() );
                lvCCDirection[lvcont] = qh->get_mate()->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline;
                lvCCDirection[lvcont+1] = qh->get_mate()->get_edge()->get_attributes().geoline;

                qh = qh->get_next();
            }

            lvcont += 2;

        }

        ASSERT ( lvcont == 8 );

        ASSERT ( qh == qf->get_halfedge() );

        // Store the fifth QVertex for the split method
        middleVertex.push_back(vertexQuad);

        // 2D after the geodesic cut

        //        double *uv_ = new double[3*vertices_uv.size()];

        //        int l = 0;

        //        int u = 0;

        //        for (std::set< pack3<Vertex*,double,double> >::iterator it = vertices_uv.begin() ;
        //             it != vertices_uv.end() ; it++){

        //            Vertex* v = it->first;
        //            v->get_attributes().tempInt = u;
        //            u++;

        //            uv_[3*l] = it->second;
        //            uv_[3*l+1] = it->third;
        //            uv_[3*l+2] = 0.0;
        //            l++;

        //        }

        //        l = 0;
        //        unsigned* fset = new unsigned[ 3*qf->get_attributes().faces.size()];
        //        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {
        //            Face* f = *it;
        //            Vertex* v0 = f->get_halfedge()->get_origin();
        //            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
        //            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

        //            fset[3*l] = v0->get_attributes().tempInt;
        //            fset[3*l+1] = v1->get_attributes().tempInt;
        //            fset[3*l+2] = v2->get_attributes().tempInt;
        //            l++;
        //        }

        //vtkReader reader;

        //vtkWriter::saveVTKTriPatchWithoutTvalue("Teste", vertices_uv.size(), uv_, qf->get_attributes().faces.size(), fset);

        qop.split_template2(qh,&middleVertex);

        addNewFaces( qf->get_attributes().faces, lvMiddleLeft);
        addNewFaces( qf->get_attributes().faces, lvMiddleRight);
        addNewFaces( qf->get_attributes().faces, lvMiddleBottom);
        addNewFaces( qf->get_attributes().faces, lvMiddleTop);
        addNewFaces( qf->get_attributes().faces, middleVertex[4]);

        std::set< Face* >  faces_set = qf->get_attributes().faces;
        qf->get_attributes().faces.clear();

        // Cleaning the link between face and triangles inside
        for ( std::set< Face* >::iterator it = faces_set.begin(); it !=faces_set.end(); it++ ){

            Face* f;
            f = *it;
            f->get_attributes().region.clear();
        }

        // Setting the geolines to the new Edges
        QFace* qfNeighbour;

        int j = 0;

        for ( int i = 0 ; i < 8 ; i += 2 ){


            correctOrientation( qh->get_edge(), lvCCDirection[i] );
            qh->get_edge()->get_attributes().geoline = lvCCDirection[i];
            correctOrientation( qh, lvCCDirection[i] );
            qfNeighbour = qh->get_mate()->get_face();
            edgeMark( lvCCDirection[i], qh->get_face(), qfNeighbour );

            qh = qh->get_next()->get_mate()->get_next();
            correctOrientation( qh->get_edge(), lvCCDirection[i+1] );
            qh->get_edge()->get_attributes().geoline = lvCCDirection[i+1];
            correctOrientation( qh, lvCCDirection[i+1] );
            if(vertexValidation[j])
                qfNeighbour = qh->get_mate()->get_face();
            edgeMark( lvCCDirection[i+1], qh->get_face(), qfNeighbour);

            qh = qh->get_next();

            j++;
        }

        ASSERT ( j == 4 );

        if( (debug::iteration==6)&&(debug::templatenumber==1530) )
            int a= 0;

        correctOrientation( qh->get_next()->get_edge(), lvMiddleLeft );
        qh->get_next()->get_edge()->get_attributes().geoline = lvMiddleLeft;
        //correctOrientation( qh->get_next(), lvMiddleLeft );
        edgeMark(lvMiddleLeft,qh->get_next()->get_edge());

        correctOrientation( qh->get_next()->get_mate()->get_prev()->get_edge(), lvMiddleBottom );
        qh->get_next()->get_mate()->get_prev()->get_edge()->get_attributes().geoline = lvMiddleBottom;
        //correctOrientation( qh->get_next()->get_mate()->get_prev(), lvMiddleBottom );
        edgeMark(lvMiddleBottom,qh->get_next()->get_mate()->get_prev()->get_edge());

        correctOrientation( qh->get_next()->get_next()->get_edge(), lvMiddleTop );
        qh->get_next()->get_next()->get_edge()->get_attributes().geoline = lvMiddleTop;
        //correctOrientation( qh->get_next()->get_next(), lvMiddleTop );
        edgeMark(lvMiddleTop,qh->get_next()->get_next()->get_edge());

        correctOrientation( qh->get_next()->get_next()->get_mate()->get_next()->get_edge(), lvMiddleRight );
        qh->get_next()->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline = lvMiddleRight;
        //correctOrientation( qh->get_next()->get_next()->get_mate()->get_next(), lvMiddleRight );
        edgeMark(lvMiddleRight,qh->get_next()->get_next()->get_mate()->get_next()->get_edge());

        int i = 0;
        for ( std::set< Face* >::iterator it = faces_set.begin(); it !=faces_set.end(); it++ ){

            Face* f;
            f = *it;
            if( f->get_attributes().region.empty() ) {
                flood_region(f);
                i++;
            }
        }
        //ASSERT(i==4);

        ASSERT( faces_set.size() == (qh->get_face()->get_attributes().faces.size() +
                                     qh->get_next()->get_mate()->get_face()->get_attributes().faces.size() +
                                     qh->get_next()->get_next()->get_mate()->get_face()->get_attributes().faces.size() +
                                     qh->get_next()->get_mate()->get_prev()->get_mate()->get_face()->get_attributes().faces.size() ));

        for(i=0; i<4; i++) {
            qh->get_edge()->get_attributes().subdivide = !vertexValidation[i];
            qh = qh->get_next()->get_mate()->get_next();
            qh->get_edge()->get_attributes().subdivide = !vertexValidation[i];
            qh = qh->get_next();
        }


    }

    // ****************************************************************
    //                    MAKING TESTS TEMPLATE 3
    // ****************************************************************

    /**
      * \fn void makingTestsTemplate3(TMesh* tm, QMesh* qm, geodesicMeshing* geo)
      *
      * \brief "Main" file of a test case using the template 3 ( D - D Twisted ).
      *
      * \param tm Triangular mesh
      *
      * \param qm Quadrilateral mesh
      *
      * \param geo Pointer to the geodesic object.
      */

    void makingTestsTemplate3(TMesh* tm, QMesh* qm, QHalfedge* qh, geodesicMeshing *geo){

        // Template 3

        QFace* qf;
        std::vector<QVertex*> middleVertex;
        QVertex* vertexQuad;
        QOperators qop(qm);
        std::vector<bool> vertexValidation;
        std::list< Vertex* >::iterator vertexIterators[4];
        //      std::list< QHalfedge* > _list_of_qhalfedges;

        Vertex* middleTriVertex;
        std::list< Vertex* > lvTri;

        //faceQFacePointerVerification(tm, qm);

        //        // Retrieving a star from the vertex

        //        qf = *(qm->faces_begin());
        //        QHalfedge* qhBegin = qf->get_halfedge();

        //        do {

        //            _list_of_qhalfedges.push_back(qhBegin);
        //            qhBegin = qhBegin->get_prev()->get_mate();

        //        } while ( qhBegin != qf->get_halfedge() );

        qf = qh->get_face();

        middleVertex.clear();
        vertexValidation.clear();

        checkNeighboursTemplate3(middleVertex, vertexIterators, vertexValidation, qh, qm, geo);

        Vertex* vi = dynamic_cast<geodesics::GPV* >( middleVertex.at(0)->get_attributes().getGP())->get_vertex();
        Vertex* ve = dynamic_cast<geodesics::GPV* >( qh->get_next()->get_next()->get_origin()->get_attributes().getGP())->get_vertex();
        lvTri.clear();

        std::set< pack3<Vertex*,double,double> > vertices_uv;
        floater(qf,vertices_uv);
        getGeodesics(qf, vi, ve, lvTri, vertices_uv, geo);

        lvTri.push_front(vi);
        lvTri.push_back(ve);

        fixingColor(lvTri, qf);
        addNewFaces( qf->get_attributes().faces, lvTri);

        // Creating QVertex with the middle Vertex coordinates
        std::list< Vertex* >::iterator itMiddle;

        itMiddle = findMiddleVertex(lvTri, geo, &vertices_uv);
        middleTriVertex = *itMiddle;

        vertexQuad = add_vertex(middleTriVertex->x(),
                                middleTriVertex->y(),
                                middleTriVertex->z(), qm);

        vertexQuad->get_attributes().setGP(new geodesics::GPV(middleTriVertex));

        // Creating the geolines of the new Edges
        std::list< Vertex* > lvMiddleLeft;
        lvMiddleLeft.splice(lvMiddleLeft.begin(),lvTri,lvTri.begin(),itMiddle);
        lvMiddleLeft.push_back(*itMiddle);

        std::list< Vertex* > lvMiddleRight;
        lvMiddleRight.splice(lvMiddleRight.begin(),lvTri,itMiddle,lvTri.end());

        std::list< Vertex* > lvMiddleTop;
        getGeodesics(qf, middleTriVertex,dynamic_cast<geodesics::GPV* >( middleVertex[1]->get_attributes().getGP())->get_vertex(),lvMiddleTop, vertices_uv, geo);

        lvMiddleTop.push_front(middleTriVertex);
        lvMiddleTop.push_back(dynamic_cast<geodesics::GPV* >( middleVertex[1]->get_attributes().getGP())->get_vertex());
        fixingColor(lvMiddleTop, qf);
        addNewFaces( qf->get_attributes().faces, lvMiddleTop);

        // Splitting the existing geolines of the old Edges

        std::list< Vertex* > lvLeftTop;
        std::list< Vertex* > lvLeftBottom;

        if ( vertexValidation[0] == false ){

            std::list< Vertex* >&lv = qh->get_edge()->get_attributes().geoline;
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_first_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.front() );
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_second_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.back() );

            lvLeftTop.splice(lvLeftTop.begin(), lv, lv.begin(), vertexIterators[0]);
            lvLeftTop.push_back(*vertexIterators[0]);
            lvLeftBottom.splice(lvLeftBottom.begin(), lv, vertexIterators[0], lv.end());

            ASSERT( lv.empty() );

            if( qh->get_edge()->get_second_halfedge() == qh )
            {
                lvLeftTop.swap( lvLeftBottom );
                lvLeftTop.reverse();
                lvLeftBottom.reverse();
            }

        } else {

            ASSERT( qh->get_mate()->get_next()->get_mate()->get_next()->get_next()->get_origin() == qh->get_origin() );
            lvLeftTop = qh->get_mate()->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline;
            lvLeftBottom = qh->get_mate()->get_edge()->get_attributes().geoline;
        }


        qh = qh->get_prev();

        std::list< Vertex* > lvTopRight;
        std::list< Vertex* > lvTopLeft;

        if ( vertexValidation[1] == false ){

            std::list< Vertex* >&lv = qh->get_edge()->get_attributes().geoline;
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_first_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.front() );
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_second_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.back() );

            lvTopRight.splice(lvTopRight.begin(), lv, lv.begin(), vertexIterators[3]);
            lvTopRight.push_back(*vertexIterators[3]);
            lvTopLeft.splice(lvTopLeft.begin(), lv, vertexIterators[3], lv.end());

            ASSERT( lv.empty() );

            if( qh->get_edge()->get_second_halfedge() == qh )
            {
                lvTopRight.swap( lvTopLeft );
                lvTopRight.reverse();
                lvTopLeft.reverse();
            }

        } else {

            ASSERT( qh->get_mate()->get_next()->get_mate()->get_next()->get_next()->get_origin() == qh->get_origin() );
            lvTopRight = qh->get_mate()->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline;
            lvTopLeft = qh->get_mate()->get_edge()->get_attributes().geoline;

        }

        qh = qh->get_next();

        //ASSERT ( qh == qf->get_halfedge() );

        middleVertex.push_back(vertexQuad);

        // 2D after the geodesic cut

        //        double *uv_ = new double[3*vertices_uv.size()];

        //        int l = 0;

        //        int u = 0;

        //        for (std::set< pack3<Vertex*,double,double> >::iterator it = vertices_uv.begin() ;
        //             it != vertices_uv.end() ; it++){

        //            Vertex* v = it->first;
        //            v->get_attributes().tempInt = u;
        //            u++;

        //            uv_[3*l] = it->second;
        //            uv_[3*l+1] = it->third;
        //            uv_[3*l+2] = 0.0;
        //            l++;

        //        }

        //        l = 0;

        //        unsigned* fset = new unsigned[ 3*qf->get_attributes().faces.size()];
        //        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {
        //            Face* f = *it;
        //            Vertex* v0 = f->get_halfedge()->get_origin();
        //            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
        //            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

        //            fset[3*l] = v0->get_attributes().tempInt;
        //            fset[3*l+1] = v1->get_attributes().tempInt;
        //            fset[3*l+2] = v2->get_attributes().tempInt;
        //            l++;
        //        }

        //vtkReader reader;

        //reader.saveVTKTriPatchWithoutTvalue("Abracadabra", vertices_uv.size(), uv_, qf->get_attributes().faces.size(), fset);

        qop.split_template3(qh,&middleVertex);

        addNewFaces( qf->get_attributes().faces, lvMiddleLeft);
        addNewFaces( qf->get_attributes().faces, lvMiddleRight);
        addNewFaces( qf->get_attributes().faces, lvMiddleTop);
        addNewFaces( qf->get_attributes().faces, vertexQuad);

        std::set< Face* >  faces_set = qf->get_attributes().faces;
        qf->get_attributes().faces.clear();

        // Cleaning the link between face and triangles inside
        for ( std::set< Face* >::iterator it = faces_set.begin(); it !=faces_set.end(); it++ ){

            Face* f;
            f = *it;
            f->get_attributes().region.clear();
        }

        ASSERT( qh == qf->get_halfedge() );

        // Setting the geolines to the new Edges
        QFace* qfNeighbour;

        // lvLeftTop
        correctOrientation( qh->get_edge(), lvLeftTop );
        qh->get_edge()->get_attributes().geoline = lvLeftTop;
        correctOrientation( qh, lvLeftTop );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark( lvLeftTop, qh->get_face(), qfNeighbour );

        // lvLeftBottom
        qh = qh->get_next()->get_mate()->get_next();
        correctOrientation( qh->get_edge(), lvLeftBottom );
        qh->get_edge()->get_attributes().geoline = lvLeftBottom;
        correctOrientation( qh, lvLeftBottom );
        if(vertexValidation[0])
            qfNeighbour = qh->get_mate()->get_face();
        edgeMark( lvLeftBottom, qh->get_face(), qfNeighbour);

        // lvBottom
        qh = qh->get_next();
        std::list< Vertex* > lvBottom = qh->get_edge()->get_attributes().geoline;
        correctOrientation( qh->get_edge(), lvBottom );
        qh->get_edge()->get_attributes().geoline = lvBottom;
        correctOrientation( qh, lvBottom );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvBottom, qh->get_face(), qfNeighbour);

        // lvRight
        qh = qh->get_next()->get_mate()->get_next();
        std::list< Vertex* > lvRight = qh->get_edge()->get_attributes().geoline;
        correctOrientation( qh->get_edge(), lvRight );
        qh->get_edge()->get_attributes().geoline = lvRight;
        correctOrientation( qh, lvRight );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvRight, qh->get_face(), qfNeighbour);

        // lvTopRight
        qh = qh->get_next();
        correctOrientation( qh->get_edge(), lvTopRight );
        qh->get_edge()->get_attributes().geoline = lvTopRight;
        correctOrientation( qh, lvTopRight );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark( lvTopRight, qh->get_face(), qfNeighbour );

        // lvTopLeft
        qh = qh->get_next()->get_mate()->get_next();
        correctOrientation( qh->get_edge(), lvTopLeft );
        qh->get_edge()->get_attributes().geoline = lvTopLeft;
        correctOrientation( qh, lvTopLeft );
        if(vertexValidation[1])
            qfNeighbour = qh->get_mate()->get_face();
        edgeMark( lvTopLeft, qh->get_face(), qfNeighbour);

        qh = qh->get_next()->get_next();

        // lvMiddleLeft
        correctOrientation( qh->get_edge(), lvMiddleLeft );
        qh->get_edge()->get_attributes().geoline = lvMiddleLeft;
        edgeMark(lvMiddleLeft,qh->get_edge());

        // lvMiddleRight
        correctOrientation( qh->get_next()->get_mate()->get_next()->get_edge(), lvMiddleRight );
        qh->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline = lvMiddleRight;
        edgeMark(lvMiddleRight,qh->get_next()->get_mate()->get_next()->get_edge());

        // lvMiddleTop
        correctOrientation( qh->get_next()->get_edge(), lvMiddleTop );
        qh->get_next()->get_edge()->get_attributes().geoline = lvMiddleTop;
        edgeMark(lvMiddleTop,qh->get_next()->get_edge());

        qh = qh->get_prev();

        int i = 0;
        for ( std::set< Face* >::iterator it = faces_set.begin(); it !=faces_set.end(); it++ ){

            Face* f;
            f = *it;
            if( f->get_attributes().region.empty() ) {
                flood_region(f);
                i++;
            }
        }

        ASSERT( faces_set.size() == (qh->get_face()->get_attributes().faces.size() +
                                     qh->get_next()->get_mate()->get_face()->get_attributes().faces.size() +
                                     qh->get_next()->get_next()->get_mate()->get_face()->get_attributes().faces.size() ));


        qh->get_edge()->get_attributes().subdivide = !vertexValidation[0];
        qh = qh->get_next()->get_mate()->get_next();
        qh->get_edge()->get_attributes().subdivide = !vertexValidation[0];
        qh = qh->get_next();

        qh = qh->get_next()->get_mate()->get_next();
        qh = qh->get_next();

        qh->get_edge()->get_attributes().subdivide = !vertexValidation[1];
        qh = qh->get_next()->get_mate()->get_next();
        qh->get_edge()->get_attributes().subdivide = !vertexValidation[1];
        qh = qh->get_next();



    }
    // ****************************************************************
    //                 MAKING TESTS WITH TEMPLATE 1
    // ****************************************************************

    /**
      * \fn void makingTestsTemplate1(QMesh* qm, geodesicMeshing* geo)
      *
      * \brief "Main" file of a test case using the template 1 ( D - D Straight ).
      *
      * \param tm Pointer to the triangular mesh.
      *
      * \param qm Pointer to the quadrilateral mesh.
      *
      * \param geo Pointer to the geodesic object.
      */

    void makingTestsTemplate1(TMesh* tm, QMesh* qm, QHalfedge* qh, geodesicMeshing* geo){

        // Template 1

        QFace* qf;
        std::vector<QVertex*> middleVertex;
        QOperators qop(qm);
        std::vector<bool> vertexValidation;
        std::list< Vertex* >::iterator vertexIterators[4];

        //      std::list< QHalfedge* > _list_of_qhalfedges;
        std::list< Vertex* > lvTri;

        //        // Retrieving a strip from the quad

        //        qf = *(qm->faces_begin());
        //        QHalfedge* qhBegin = qf->get_halfedge();

        //        do {

        //            _list_of_qhalfedges.push_back(qhBegin);
        //            qhBegin = qhBegin->get_next()->get_mate()->get_next();

        //        } while ( qhBegin != qf->get_halfedge() );

        // faceQFacePointerVerification(tm, qm);

        qf = qh->get_face();

        middleVertex.clear();
        vertexValidation.clear();

        checkNeighboursTemplate1(middleVertex, vertexIterators, vertexValidation, qh, qm, geo);

        // Get geodesic from middleVertex[1] to middleVertex[3]
        Vertex* vi = dynamic_cast<geodesics::GPV* >( middleVertex.at(0)->get_attributes().getGP())->get_vertex();
        Vertex* ve = dynamic_cast<geodesics::GPV* >( middleVertex.at(1)->get_attributes().getGP())->get_vertex();
        lvTri.clear();

        std::set< pack3<Vertex*,double,double> > vertices_uv;
        floater(qf,vertices_uv);

        getGeodesics(qf, vi, ve, lvTri, vertices_uv, geo);

        lvTri.push_front(vi);
        lvTri.push_back(ve);

        fixingColor(lvTri, qf);
        addNewFaces( qf->get_attributes().faces, lvTri);

        // Splitting the existing geolines of the old Edges

        qh = qh->get_next();

        // If Bottom Edge isn't subdivided ..

        std::list< Vertex* > lvBottomLeft;
        std::list< Vertex* > lvBottomRight;

        if ( vertexValidation[0] == false ){

            std::list< Vertex* >&lv = qh->get_edge()->get_attributes().geoline;
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_first_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.front() );
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_second_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.back() );

            lvBottomLeft.splice(lvBottomLeft.begin(),lv, lv.begin(),vertexIterators[1]);

            lvBottomLeft.push_back(*vertexIterators[1]);

            lvBottomRight.splice(lvBottomRight.begin(),lv, vertexIterators[1], lv.end());

            ASSERT( lv.empty() );

            if( qh->get_edge()->get_second_halfedge() == qh )
            {
                lvBottomLeft.swap( lvBottomRight );
                lvBottomLeft.reverse();
                lvBottomRight.reverse();
            }

        } else {

            ASSERT( qh->get_mate()->get_next()->get_mate()->get_next()->get_next()->get_origin() == qh->get_origin() );
            lvBottomLeft = qh->get_mate()->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline;
            lvBottomRight = qh->get_mate()->get_edge()->get_attributes().geoline;

        }

        qh = qh->get_next()->get_next();

        // ************************************************

        // If top edge isn't subdivided ..

        std::list< Vertex* > lvTopLeft;
        std::list< Vertex* > lvTopRight;

        if ( vertexValidation[1] == false ){

            std::list< Vertex* >& lv = qh->get_edge()->get_attributes().geoline;
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_first_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.front() );
            ASSERT( dynamic_cast<geodesics::GPV* >(qh->get_edge()->get_second_halfedge()->get_origin()->get_attributes().getGP())->get_vertex() == lv.back() );

            lvTopRight.splice(lvTopRight.begin(),lv,lv.begin(),vertexIterators[3]);

            lvTopRight.push_back(*vertexIterators[3]);

            lvTopLeft.splice(lvTopLeft.begin(),lv,vertexIterators[3],lv.end());

            ASSERT( lv.empty() );

            if( qh->get_edge()->get_second_halfedge() == qh )
            {
                lvTopRight.swap( lvTopLeft );
                lvTopRight.reverse();
                lvTopLeft.reverse();
            }

        } else {

            ASSERT( qh->get_mate()->get_next()->get_mate()->get_next()->get_next()->get_origin() == qh->get_origin() );
            lvTopRight = qh->get_mate()->get_next()->get_mate()->get_next()->get_edge()->get_attributes().geoline;
            lvTopLeft = qh->get_mate()->get_edge()->get_attributes().geoline;
        }

        qh = qh->get_next();

        // 2D after the geodesic cut

        //        double *uv_ = new double[3*vertices_uv.size()];

        //        int l = 0;

        //        int u = 0;

        //        for (std::set< pack3<Vertex*,double,double> >::iterator it = vertices_uv.begin() ;
        //             it != vertices_uv.end() ; it++){

        //            Vertex* v = it->first;
        //            v->get_attributes().tempInt = u;
        //            u++;

        //            uv_[3*l] = it->second;
        //            uv_[3*l+1] = it->third;
        //            uv_[3*l+2] = 0.0;
        //            l++;

        //        }

        //        l = 0;

        //        unsigned* fset = new unsigned[ 3*qf->get_attributes().faces.size()];
        //        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {
        //            Face* f = *it;
        //            Vertex* v0 = f->get_halfedge()->get_origin();
        //            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
        //            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

        //            fset[3*l] = v0->get_attributes().tempInt;
        //            fset[3*l+1] = v1->get_attributes().tempInt;
        //            fset[3*l+2] = v2->get_attributes().tempInt;
        //            l++;
        //        }

        // vtkReader reader;

        //reader.saveVTKTriPatchWithoutTvalue("Teste", vertices_uv.size(), uv_, qf->get_attributes().faces.size(), fset);

        // Splitting the face
        qop.split_template1(qh,&middleVertex);

        addNewFaces( qf->get_attributes().faces, lvTri);

        std::set< Face* >  faces_set = qf->get_attributes().faces;
        qf->get_attributes().faces.clear();

        // Cleaning the link between face and triangles inside
        for ( std::set< Face* >::iterator it = faces_set.begin(); it !=faces_set.end(); it++ ){

            Face* f;
            f = *it;
            f->get_attributes().region.clear();
        }

        // Setting the geolines to the new Edges
        QFace* qfNeighbour;

        // lvLeft

        std::list< Vertex* > lvLeft = qh->get_edge()->get_attributes().geoline;
        correctOrientation( qh->get_edge(), lvLeft);
        qh->get_edge()->get_attributes().geoline = lvLeft;
        correctOrientation( qh, lvLeft );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvLeft, qh->get_face(), qfNeighbour);

        // Verification of Bottom Edge

        // lvBottom

        qh = qh->get_next();

        correctOrientation( qh->get_edge(), lvBottomLeft );
        qh->get_edge()->get_attributes().geoline = lvBottomLeft;
        correctOrientation( qh, lvBottomLeft );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvBottomLeft, qh->get_face(), qfNeighbour);

        qh = qh->get_next()->get_mate()->get_next();
        correctOrientation( qh->get_edge(), lvBottomRight );
        qh->get_edge()->get_attributes().geoline = lvBottomRight;
        correctOrientation( qh, lvBottomRight );
        if(vertexValidation[0])
            qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvBottomRight,qh->get_face(), qfNeighbour);

        // lvRight

        qh = qh->get_next();

        std::list< Vertex* > lvRight = qh->get_edge()->get_attributes().geoline;
        correctOrientation( qh->get_edge(), lvRight);
        qh->get_edge()->get_attributes().geoline = lvRight;
        correctOrientation( qh, lvRight );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvRight, qh->get_face(), qfNeighbour);

        // Verification of Top Edge

        // lvTop
        qh = qh->get_next();
        correctOrientation( qh->get_edge(), lvTopRight );
        qh->get_edge()->get_attributes().geoline = lvTopRight;
        correctOrientation( qh, lvTopRight );
        qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvTopRight,qh->get_face(), qfNeighbour);

        qh = qh->get_next()->get_mate()->get_next();
        correctOrientation( qh->get_edge(), lvTopLeft );
        qh->get_edge()->get_attributes().geoline = lvTopLeft;
        correctOrientation( qh, lvTopLeft );
        if(vertexValidation[1])
            qfNeighbour = qh->get_mate()->get_face();
        edgeMark(lvTopLeft,qh->get_face(), qfNeighbour);

        //***************************

        correctOrientation( qh->get_prev()->get_edge(), lvTri );
        qh->get_prev()->get_edge()->get_attributes().geoline = lvTri;
        //correctOrientation( qh->get_next(), lvMiddleLeft );
        edgeMark(lvTri,qh->get_prev()->get_edge());

        qh = qh->get_next();

        int i = 0;
        for ( std::set< Face* >::iterator it = faces_set.begin(); it !=faces_set.end(); it++ ){

            Face* f;
            f = *it;
            if( f->get_attributes().region.empty() ) {
                flood_region(f);
                i++;
            }
        }

        ASSERT( faces_set.size() == (qh->get_face()->get_attributes().faces.size() +
                                     qh->get_next()->get_next()->get_mate()->get_face()->get_attributes().faces.size() ));


        qh = qh->get_next();
        qh->get_edge()->get_attributes().subdivide = !vertexValidation[0];
        qh = qh->get_next()->get_mate()->get_next();
        qh->get_edge()->get_attributes().subdivide = !vertexValidation[0];

        qh = qh->get_next()->get_next();
        qh->get_edge()->get_attributes().subdivide = !vertexValidation[1];
        qh = qh->get_next()->get_mate()->get_next();
        qh->get_edge()->get_attributes().subdivide = !vertexValidation[1];

        qh = qh->get_next();

    }

    void parametrization( QFace* qf, std::set< pack3<Vertex*,double,double> >& vertices_uv ) {

        ASSERT( qf );
        ASSERT( ! qf->get_attributes().faces.empty() );
        //ASSERT( vertices_uv.empty() );

        std::set< Vertex* > vertices;

        for( std::set< Face* >::iterator it = qf->get_attributes().faces.begin(); it != qf->get_attributes().faces.end(); it++ ) {

            Face* f = *it;
            Vertex* v0 = f->get_halfedge()->get_origin();
            Vertex* v1 = f->get_halfedge()->get_next()->get_origin();
            Vertex* v2 = f->get_halfedge()->get_prev()->get_origin();

            vertices.insert( v0 );
            vertices.insert( v1 );
            vertices.insert( v2 );
        }

        cgalMethods::parametrization(vertices, qf->get_attributes().faces, vertices_uv);
        //cm.path(path1);
    }

};

#endif // TOOLS_H
